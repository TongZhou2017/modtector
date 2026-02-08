use chrono;
use clap::Args;
use quick_xml::events::Event;
use quick_xml::Reader;
use std::error::Error;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Read, Write};
use std::path::Path;
use std::time::Instant;

/// Validate convert command parameters
fn validate_convert_args(args: &ConvertArgs) -> Result<(), Box<dyn Error>> {
    // Check if using dual input mode (mutation + stop)
    let using_dual_input = args.input_mutation.is_some() || args.input_stop.is_some();
    
    if using_dual_input {
        // Dual input mode: at least one of input_mutation or input_stop must be provided
        // Files are optional - if missing, corresponding columns will be empty
        if let Some(ref path) = args.input_mutation {
            if !path.trim().is_empty() && !Path::new(path).exists() {
                return Err(format!("Error: Mutation input file does not exist: {}", path).into());
            }
        }
        if let Some(ref path) = args.input_stop {
            if !path.trim().is_empty() && !Path::new(path).exists() {
                return Err(format!("Error: Stop input file does not exist: {}", path).into());
            }
        }
    } else {
        // Single input mode: validate main input file
        if args.input.trim().is_empty() {
            return Err("Error: Input file path cannot be empty (or use --input-mutation/--input-stop)".into());
        }
        if !Path::new(&args.input).exists() {
            return Err(format!("Error: Input file does not exist: {}", args.input).into());
        }
    }

    // Validate output file
    if args.output.trim().is_empty() {
        return Err("Error: Output file path cannot be empty".into());
    }
    if !args.output.ends_with(".csv") {
        return Err(format!(
            "Error: Output file path must end with .csv: {}",
            args.output
        )
        .into());
    }

    // Validate strand
    match args.strand.as_str() {
        "+" | "-" => {}
        _ => {
            return Err(format!(
                "Error: Invalid strand value: {}. Must be '+' or '-'",
                args.strand
            )
            .into());
        }
    }

    Ok(())
}

#[derive(Args, Debug)]
pub struct ConvertArgs {
    /// Input file (bamreadcount txt format or rf-rctools TSV format)
    /// When using dual input mode (--input-mutation/--input-stop), this can be empty
    #[arg(short = 'i', long = "input")]
    pub input: String,
    /// Output pileup CSV file
    #[arg(short = 'o', long = "output")]
    pub output: String,
    /// Input format: 'bamreadcount', 'rf-rctools', 'rf-norm', 'rf-norm-xml', 'shapemapper2', 'shapemapper-profile', 'samtools-mpileup', 'icSHAPE-rt', or 'bedgraph' (default: auto-detect)
    #[arg(short = 'f', long = "format")]
    pub format: Option<String>,
    /// For rf-rctools format: if true, input was generated with rf-count -m (mutation mode)
    /// When true, field[1] is mutation_count; when false, field[1] is RT_stops
    #[arg(long = "rf-count-mutations", default_value = "false")]
    pub rf_count_mutations: bool,
    /// Strand orientation ('+' for forward, '-' for reverse)
    #[arg(short = 's', long = "strand", default_value = "+")]
    pub strand: String,
    /// Log file path (optional)
    #[arg(short = 'l', long = "log")]
    pub log: Option<String>,
    /// For shapemapper2 format: reference FASTA file path (required for shapemapper2 format)
    #[arg(long = "ref-fasta")]
    pub ref_fasta: Option<String>,
    /// For shapemapper2 format: chromosome/contig name (default: inferred from reference or "unknown")
    #[arg(long = "chr-name")]
    pub chr_name: Option<String>,
    /// Optional: Mutation input file (rf-rctools format with -m flag)
    /// When provided with --input-stop, both files will be merged into output
    /// If file is missing, mutation columns will be set to 0
    #[arg(long = "input-mutation")]
    pub input_mutation: Option<String>,
    /// Optional: Stop input file (rf-rctools format without -m flag)
    /// When provided with --input-mutation, both files will be merged into output
    /// If file is missing, stop columns will be set to 0
    #[arg(long = "input-stop")]
    pub input_stop: Option<String>,
    /// For icSHAPE-rt format: filter by strand ('+' or '-' or 'None' to process all strands)
    /// Default: '+' (only process positive strand)
    #[arg(long = "filter-strand", default_value = "+")]
    pub filter_strand: Option<String>,
    /// For rf-norm-xml format: signal type for reactivity column mapping.
    /// 'stop' = Rouskin RT-stop → write to reactivity_stop, reactivity_mutation=NaN;
    /// 'mut' = Zubradt mutation → write to reactivity_mutation, reactivity_stop=NaN (default).
    #[arg(long = "rf-norm-signal", default_value = "mut")]
    pub rf_norm_signal: String,
}

/// Parse a line from bamreadcount format.
/// Returns: (chr, position, ref_base, depth, base_a, base_c, base_g, base_t, ins, del)
fn parse_bamreadcount_line(
    line: &str,
) -> Result<(String, u64, char, u64, u64, u64, u64, u64, u64, u64), Box<dyn Error>> {
    let fields: Vec<&str> = line.split('\t').collect();
    if fields.len() < 4 {
        return Err("Invalid bam-readcount format: insufficient fields".into());
    }

    let chr = fields[0].to_string();
    let position: u64 = fields[1].parse()?;
    let reference_base = fields[2]
        .chars()
        .next()
        .ok_or("Invalid reference base")?;
    let depth: u64 = fields[3].parse()?;

    // Parse base counts from remaining fields (base:count:...)
    let mut base_a = 0u64;
    let mut base_c = 0u64;
    let mut base_g = 0u64;
    let mut base_t = 0u64;
    let mut ins = 0u64;
    let mut del = 0u64;

    // Process each base field (starting from field 4)
    for i in 4..fields.len() {
        let base_field = fields[i];
        let parts: Vec<&str> = base_field.split(':').collect();
        if parts.len() < 2 {
            continue;
        }

        let base = parts[0];
        let count: u64 = parts[1].parse().unwrap_or(0);

        match base {
            "A" => base_a = count,
            "C" => base_c = count,
            "G" => base_g = count,
            "T" | "U" => base_t = count,
            "+" => ins = count, // Insertion
            "-" => del = count, // Deletion
            _ => {}
        }
    }

    Ok((
        chr,
        position,
        reference_base,
        depth,
        base_a,
        base_c,
        base_g,
        base_t,
        ins,
        del,
    ))
}

/// Parse a line from rf-rctools view format.
/// When using rf-count without -m: TSV format is (base, RT_stops, coverage)
/// When using rf-count with -m: TSV format is (base, mutation_count, coverage)
/// Returns: (ref_base, mutation_count, rt_stops, coverage)
/// Note: When -m is used, mutation_count is from field[1], rt_stops = 0
///       When -m is not used, rt_stops is from field[1], mutation_count = 0
fn parse_rf_rctools_line(
    line: &str,
    _chr: &str,
    _position: u64,
    count_mutations: bool,
) -> Result<(char, u64, u64, u64), Box<dyn Error>> {
    let fields: Vec<&str> = line.split('\t').collect();
    if fields.len() < 3 {
        return Err("Invalid rf-rctools format: insufficient fields".into());
    }

    let ref_base = fields[0]
        .chars()
        .next()
        .ok_or("Invalid reference base")?;
    let coverage: u64 = fields[2].parse()?;
    
    let (mutation_count, rt_stops) = if count_mutations {
        // When -m is used, field[1] is mutation count
        let mutation_count: u64 = fields[1].parse()?;
        (mutation_count, 0)
    } else {
        // When -m is not used, field[1] is RT stops
        let rt_stops: u64 = fields[1].parse()?;
        (0, rt_stops)
    };

    Ok((ref_base, mutation_count, rt_stops, coverage))
}

/// Calculate mutation count (non-reference bases)
fn calculate_mutation_count(
    ref_base: char,
    base_a: u64,
    base_c: u64,
    base_g: u64,
    base_t: u64,
) -> u64 {
    let ref_base_upper = ref_base.to_ascii_uppercase();
    match ref_base_upper {
        'A' => base_c + base_g + base_t,
        'C' => base_a + base_g + base_t,
        'G' => base_a + base_c + base_t,
        'T' | 'U' => base_a + base_c + base_g,
        _ => 0,
    }
}

/// Parse a line from samtools mpileup format.
/// Format: chr\tposition\tref_base\tstrand\tdepth\tA\tC\tG\tT\tN
/// Returns: (chr, position, ref_base, strand, depth, base_a, base_c, base_g, base_t, base_n)
fn parse_samtools_mpileup_line(
    line: &str,
) -> Result<(String, u64, char, String, u64, u64, u64, u64, u64, u64), Box<dyn Error>> {
    let fields: Vec<&str> = line.split('\t').collect();
    if fields.len() < 10 {
        return Err("Invalid samtools mpileup format: expected at least 10 columns (chr, position, ref_base, strand, depth, A, C, G, T, N)".into());
    }

    let chr = fields[0].to_string();
    let position: u64 = fields[1].parse()?;
    let ref_base = fields[2]
        .chars()
        .next()
        .ok_or("Invalid reference base")?;
    let strand = fields[3].to_string();
    let depth: u64 = fields[4].parse()?;
    let base_a: u64 = fields[5].parse().unwrap_or(0);
    let base_c: u64 = fields[6].parse().unwrap_or(0);
    let base_g: u64 = fields[7].parse().unwrap_or(0);
    let base_t: u64 = fields[8].parse().unwrap_or(0);
    let base_n: u64 = fields[9].parse().unwrap_or(0);

    Ok((
        chr,
        position,
        ref_base,
        strand,
        depth,
        base_a,
        base_c,
        base_g,
        base_t,
        base_n,
    ))
}

/// Parse a line from shapemapper2 mutations.txt format.
/// Returns: (mutation_count, depth, base_a, base_c, base_g, base_t, ins, del)
/// shapemapper2 format columns (35 total):
/// 0: A-, 1: T-, 2: G-, 3: C- (deletions)
/// 4: -A, 5: -T, 6: -G, 7: -C, 8: -N (insertions)
/// 9-20: AT, AG, AC, TA, TG, TC, GA, GT, GC, CA, CT, CG (dinucleotide mutations = single nucleotide substitutions)
/// 21: multinuc_deletion, 22: multinuc_insertion
/// 23-27: A_multinuc_mismatch, C_multinuc_mismatch, G_multinuc_mismatch, T_multinuc_mismatch, N_multinuc_mismatch
/// 28: complex_deletion, 29: complex_insertion
/// 30: read_depth, 31: effective_depth, 32: off_target_mapped_depth, 33: low_mapq_mapped_depth, 34: mapped_depth
/// 
/// Note: Dinucleotide mutations are actually single nucleotide substitutions:
/// - AT = A->T, AG = A->G, AC = A->C
/// - TA = T->A, TG = T->G, TC = T->C
/// - GA = G->A, GT = G->T, GC = G->C
/// - CA = C->A, CT = C->T, CG = C->G
/// 
/// This function returns base counts and total mutation count, but the actual mutation_count
/// should be calculated based on the reference base (non-reference bases only).
fn parse_shapemapper2_line(
    line: &str,
) -> Result<(u64, u64, u64, u64, u64, u64, u64, u64), Box<dyn Error>> {
    let fields: Vec<&str> = line.split('\t').collect();
    if fields.len() < 35 {
        return Err(format!(
            "Invalid shapemapper2 format: expected at least 35 columns, got {}",
            fields.len()
        )
        .into());
    }

    // Parse numeric fields
    let a_del: u64 = fields[0].parse().unwrap_or(0);  // A-
    let t_del: u64 = fields[1].parse().unwrap_or(0);   // T-
    let g_del: u64 = fields[2].parse().unwrap_or(0); // G-
    let c_del: u64 = fields[3].parse().unwrap_or(0);  // C-
    
    let a_ins: u64 = fields[4].parse().unwrap_or(0);  // -A
    let t_ins: u64 = fields[5].parse().unwrap_or(0);  // -T
    let g_ins: u64 = fields[6].parse().unwrap_or(0);  // -G
    let c_ins: u64 = fields[7].parse().unwrap_or(0);  // -C
    
    // Dinucleotide mutations (AT, AG, AC, TA, TG, TC, GA, GT, GC, CA, CT, CG)
    let at_mut: u64 = fields[9].parse().unwrap_or(0);
    let ag_mut: u64 = fields[10].parse().unwrap_or(0);
    let ac_mut: u64 = fields[11].parse().unwrap_or(0);
    let ta_mut: u64 = fields[12].parse().unwrap_or(0);
    let tg_mut: u64 = fields[13].parse().unwrap_or(0);
    let tc_mut: u64 = fields[14].parse().unwrap_or(0);
    let ga_mut: u64 = fields[15].parse().unwrap_or(0);
    let gt_mut: u64 = fields[16].parse().unwrap_or(0);
    let gc_mut: u64 = fields[17].parse().unwrap_or(0);
    let ca_mut: u64 = fields[18].parse().unwrap_or(0);
    let ct_mut: u64 = fields[19].parse().unwrap_or(0);
    let cg_mut: u64 = fields[20].parse().unwrap_or(0);
    
    // Multinucleotide mismatches
    let a_multinuc: u64 = fields[23].parse().unwrap_or(0);
    let c_multinuc: u64 = fields[24].parse().unwrap_or(0);
    let g_multinuc: u64 = fields[25].parse().unwrap_or(0);
    let t_multinuc: u64 = fields[26].parse().unwrap_or(0);
    
    // Calculate mutation counts per base
    // For A: mutations to A are from T->A, G->A, C->A, and A multinuc
    let base_a_mutations = ta_mut + ga_mut + ca_mut + a_multinuc;
    // For C: mutations to C are from A->C, T->C, G->C, and C multinuc
    let base_c_mutations = ac_mut + tc_mut + gc_mut + c_multinuc;
    // For G: mutations to G are from A->G, T->G, C->G, and G multinuc
    let base_g_mutations = ag_mut + tg_mut + cg_mut + g_multinuc;
    // For T: mutations to T are from A->T, G->T, C->T, and T multinuc
    let base_t_mutations = at_mut + gt_mut + ct_mut + t_multinuc;
    
    // Total mutation count (sum of all mutations)
    let mutation_count = base_a_mutations + base_c_mutations + base_g_mutations + base_t_mutations;
    
    // Depth (use effective_depth if available, otherwise read_depth)
    let depth: u64 = if fields.len() > 31 {
        fields[31].parse().unwrap_or_else(|_| fields[30].parse().unwrap_or(0))
    } else if fields.len() > 30 {
        fields[30].parse().unwrap_or(0)
    } else {
        0
    };
    
    // Insertions (sum of all insertions)
    let ins = a_ins + t_ins + g_ins + c_ins;
    
    // Deletions (sum of all deletions)
    let del = a_del + t_del + g_del + c_del;
    
    // Base counts: we don't have direct base counts in shapemapper2 format
    // We approximate by using depth and mutation counts
    // This is an approximation - actual base counts would require more information
    let base_a = base_a_mutations;
    let base_c = base_c_mutations;
    let base_g = base_g_mutations;
    let base_t = base_t_mutations;
    
    Ok((mutation_count, depth, base_a, base_c, base_g, base_t, ins, del))
}

/// Load reference sequence from FASTA file
fn load_reference_sequence(fasta_path: &str, chr_name: Option<&str>) -> Result<Vec<char>, Box<dyn Error>> {
    use std::collections::HashMap;
    
    let file = File::open(fasta_path)?;
    let reader = BufReader::new(file);
    let mut sequences: HashMap<String, String> = HashMap::new();
    let mut current_seq_name = String::new();
    let mut current_sequence = String::new();
    
    for line in reader.lines() {
        let line = line?;
        if line.starts_with('>') {
            // Save previous sequence
            if !current_seq_name.is_empty() {
                sequences.insert(current_seq_name.clone(), current_sequence.clone());
            }
            // Start new sequence
            current_seq_name = line[1..].split_whitespace().next().unwrap_or("").to_string();
            current_sequence.clear();
        } else {
            current_sequence.push_str(&line.trim());
        }
    }
    
    // Save last sequence
    if !current_seq_name.is_empty() {
        sequences.insert(current_seq_name.clone(), current_sequence.clone());
    }
    
    // Get the requested sequence
    let seq = if let Some(chr) = chr_name {
        sequences.get(chr)
            .ok_or_else(|| format!("Chromosome/contig '{}' not found in reference file", chr))?
    } else {
        // Use first sequence if no chr_name specified
        sequences.values().next()
            .ok_or("No sequences found in reference file")?
    };
    
    Ok(seq.chars().collect())
}

/// Auto-detect input format by examining the first few lines or file extension
fn detect_input_format(input_path: &str) -> Result<String, Box<dyn Error>> {
    // Check file extension first
    if input_path.ends_with(".xml") {
        // Try to verify it's rf-norm XML by checking the first few bytes
        let file = File::open(input_path)?;
        let mut reader = BufReader::new(file);
        let mut buffer = [0u8; 200];
        if let Ok(n) = reader.read(&mut buffer) {
            let content = String::from_utf8_lossy(&buffer[..n]);
            if content.contains("rf-norm") || content.contains("<data") && content.contains("<transcript") {
                return Ok("rf-norm-xml".to_string());
            }
        }
    }
    
    let file = File::open(input_path)?;
    let reader = BufReader::new(file);
    
    let mut lines = reader.lines();
    
    // Read first line
    if let Some(Ok(first_line)) = lines.next() {
        let trimmed = first_line.trim();
        
        // Check if it's icSHAPE-rt format (starts with @ColNum or @ChrID)
        if trimmed.starts_with("@ColNum") || trimmed.starts_with("@ChrID") {
            return Ok("icSHAPE-rt".to_string());
        }
        
        // Check if it's rf-norm format (CSV with header: ChrID,Strand,Position,Base,reactivity_stop,reactivity_mutation)
        if trimmed.contains(',') {
            let fields: Vec<&str> = trimmed.split(',').collect();
            if fields.len() >= 6 && 
               fields[0].trim().eq_ignore_ascii_case("ChrID") &&
               fields[1].trim().eq_ignore_ascii_case("Strand") &&
               fields[4].trim().contains("reactivity") {
                return Ok("rf-norm".to_string());
            }
        }
        
        // Check if it's shapemapper profile format (tab-separated with header containing Nucleotide, Sequence, Norm_profile)
        if trimmed.contains('\t') {
            let fields: Vec<&str> = trimmed.split('\t').collect();
            let fields_lower: Vec<String> = fields.iter().map(|f| f.trim().to_lowercase()).collect();
            // shapemapper profile format has columns: Nucleotide, Sequence, ..., Norm_profile
            if fields_lower.contains(&"nucleotide".to_string()) && 
               fields_lower.contains(&"sequence".to_string()) &&
               fields_lower.contains(&"norm_profile".to_string()) {
                return Ok("shapemapper-profile".to_string());
            }
            // shapemapper2 format has 35 columns with specific header pattern
            if fields.len() >= 30 && 
               (fields[0].trim() == "A-" || fields[0].contains("A-")) &&
               (fields.contains(&"read_depth") || fields.contains(&"effective_depth") || fields.contains(&"mapped_depth")) {
                return Ok("shapemapper2".to_string());
            }
        }
        
        // Check if it's a transcript name (single word, no tabs)
        if !trimmed.contains('\t') && !trimmed.contains(',') {
            // Likely rf-rctools format (first line is transcript name)
            return Ok("rf-rctools".to_string());
        }
        
        // Check if it's samtools mpileup format (tab-separated with header: chr\tposition\tref_base\tstrand\tdepth\tA\tC\tG\tT\tN)
        let fields: Vec<&str> = trimmed.split('\t').collect();
        if fields.len() >= 10 {
            // Check if it matches m1A mpileup header pattern
            if fields[0].trim().eq_ignore_ascii_case("chr") &&
               fields[1].trim().eq_ignore_ascii_case("position") &&
               fields[2].trim().eq_ignore_ascii_case("ref_base") &&
               fields[3].trim().eq_ignore_ascii_case("strand") &&
               fields[4].trim().eq_ignore_ascii_case("depth") {
                return Ok("samtools-mpileup".to_string());
            }
            // Check if it's a data line (position is numeric, strand is + or -)
            if fields.len() >= 10 && fields[1].parse::<u64>().is_ok() && 
               (fields[3].trim() == "+" || fields[3].trim() == "-") {
                return Ok("samtools-mpileup".to_string());
            }
        }
        
        // Check if it's bedgraph format
        // Old format (5 fields): chr\tstart\tend\tcoverage\tstrand
        // New format (6 fields): chr\tstart\tend\tstop_count\tdepth\tstrand
        if fields.len() >= 5 {
            // Check if start and end are numeric, and strand is + or -
            let has_valid_coords = fields[1].parse::<u64>().is_ok() && 
                                  fields[2].parse::<u64>().is_ok();
            let has_valid_strand = fields.len() >= 5 && 
                                  (fields[4].trim() == "+" || fields[4].trim() == "-");
            
            if has_valid_coords && has_valid_strand {
                // Check if it's new format (6 fields with stop_count and depth)
                if fields.len() >= 6 && 
                   fields[3].parse::<u64>().is_ok() && 
                   fields[4].parse::<u64>().is_ok() &&
                   (fields[5].trim() == "+" || fields[5].trim() == "-") {
                    return Ok("bedgraph".to_string());  // New format with stop_count and depth
                }
                // Check if it's old format (5 fields with coverage)
                else if fields.len() >= 5 && 
                        fields[3].parse::<f64>().is_ok() &&
                        (fields[4].trim() == "+" || fields[4].trim() == "-") {
                    return Ok("bedgraph".to_string());  // Old format with coverage only
                }
            }
        }
        
        // Check if it's bamreadcount format (tab-separated with position in second field)
        if fields.len() >= 4 {
            // Try to parse second field as position
            if fields[1].parse::<u64>().is_ok() {
                return Ok("bamreadcount".to_string());
            }
        }
    }
    
    // Default to bamreadcount if can't determine
    Ok("bamreadcount".to_string())
}

/// Convert input format (bamreadcount or rf-rctools) to modtector pileup CSV format
pub fn convert_bamreadcount_to_pileup(
    args: &ConvertArgs,
    logger: &mut crate::Logger,
) -> Result<(), Box<dyn Error>> {
    // Validate arguments
    validate_convert_args(args)?;

    // Determine input format
    let format = match &args.format {
        Some(f) => f.clone(),
        None => {
            logger.log("Auto-detecting input format...")?;
            detect_input_format(&args.input)?
        }
    };

    let start_time = Instant::now();

    // Record environment information and parameters
    logger.log("=== ModDetector Convert Function Log ===")?;
    logger.log(&format!("Software Version: v{}", crate::VERSION))?;
    logger.log(&format!(
        "Runtime: {}",
        chrono::Utc::now().format("%Y-%m-%d %H:%M:%S")
    ))?;
    
    // Check if using dual input mode
    let using_dual_input = args.input_mutation.is_some() || args.input_stop.is_some();
    
    if using_dual_input {
        logger.log("Using dual input mode (mutation + stop)")?;
        if let Some(ref path) = args.input_mutation {
            logger.log(&format!("Mutation input file: {}", path))?;
        }
        if let Some(ref path) = args.input_stop {
            logger.log(&format!("Stop input file: {}", path))?;
        }
    } else {
        logger.log(&format!("Input File: {}", args.input))?;
    }
    
    logger.log(&format!("Output File: {}", args.output))?;
    logger.log(&format!("Input Format: {}", format))?;
    logger.log(&format!("Strand: {}", args.strand))?;
    logger.log("Starting conversion...")?;

    // Open output file
    let output_file = File::create(&args.output)?;
    let mut writer = BufWriter::new(output_file);

    // Process based on format
    if using_dual_input {
        // Dual input mode: merge mutation and stop files
        // Write pileup header
        writeln!(
            writer,
            "ChrID,pipe_truncation_Strand,pipe_truncation_ChrPos,rf_mutation_Base,rf_mutation_Count,pipe_truncation_count,depth,rf_mutation_ins,rf_mutation_del,base_A,base_C,base_G,base_T"
        )?;
        
        merge_rf_rctools_to_pileup(
            args.input_mutation.as_ref(),
            args.input_stop.as_ref(),
            &mut writer,
            args,
            logger,
        )?;
    } else if format == "rf-norm-xml" {
        // Use streaming processing to avoid loading entire file into memory
        logger.log("Streaming input file (not loading into memory)...")?;

        // Open input file with large buffer for better NFS performance
        let file = File::open(&args.input)?;
        let _reader = BufReader::with_capacity(8 * 1024 * 1024, file); // 8MB buffer
        // For XML format, write norm header directly
        writeln!(
            writer,
            "ChrID,Strand,Position,Base,reactivity_stop,reactivity_mutation"
        )?;
        convert_rf_norm_xml_to_norm(&args.input, &mut writer, args, logger)?;
    } else if format == "shapemapper-profile" {
        // Use streaming processing to avoid loading entire file into memory
        logger.log("Streaming input file (not loading into memory)...")?;

        // Open input file with large buffer for better NFS performance
        let file = File::open(&args.input)?;
        let reader = BufReader::with_capacity(8 * 1024 * 1024, file); // 8MB buffer
        // For shapemapper profile format, write norm header directly
        writeln!(
            writer,
            "ChrID,Strand,Position,Base,reactivity_stop,reactivity_mutation"
        )?;
        convert_shapemapper_profile_to_norm(reader, &mut writer, args, logger)?;
    } else {
        // Use streaming processing to avoid loading entire file into memory
        logger.log("Streaming input file (not loading into memory)...")?;

        // Open input file with large buffer for better NFS performance
        let file = File::open(&args.input)?;
        let reader = BufReader::with_capacity(8 * 1024 * 1024, file); // 8MB buffer
        
        // Write pileup header for other formats
        writeln!(
            writer,
            "ChrID,pipe_truncation_Strand,pipe_truncation_ChrPos,rf_mutation_Base,rf_mutation_Count,pipe_truncation_count,depth,rf_mutation_ins,rf_mutation_del,base_A,base_C,base_G,base_T"
        )?;
        
        if format == "rf-rctools" {
            convert_rf_rctools_to_pileup(reader, &mut writer, args, logger)?;
        } else if format == "rf-norm" {
            convert_rf_norm_to_norm(reader, &mut writer, args, logger)?;
        } else if format == "shapemapper2" {
            convert_shapemapper2_to_pileup(reader, &mut writer, args, logger)?;
        } else if format == "samtools-mpileup" {
            convert_samtools_mpileup_to_pileup(reader, &mut writer, args, logger)?;
        } else if format == "icSHAPE-rt" {
            convert_icshape_rt_to_pileup(reader, &mut writer, args, logger)?;
        } else if format == "bedgraph" {
            convert_bedgraph_to_pileup(reader, &mut writer, args, logger)?;
        } else {
            convert_bamreadcount_to_pileup_internal(reader, &mut writer, args, logger)?;
        }
    }

    // Flush output
    writer.flush()?;

    let elapsed = start_time.elapsed();
    logger.log(&format!(
        "Conversion completed in {:.2} seconds",
        elapsed.as_secs_f64()
    ))?;
    logger.log(&format!("Output file: {}", args.output))?;

    Ok(())
}

/// Convert bamreadcount format to pileup (internal function)
fn convert_bamreadcount_to_pileup_internal(
    reader: BufReader<File>,
    writer: &mut BufWriter<File>,
    args: &ConvertArgs,
    logger: &mut crate::Logger,
) -> Result<(), Box<dyn Error>> {
    // Process file line by line
    let mut line_count = 0u64;
    let mut processed_count = 0u64;
    let mut error_count = 0u64;

    for line in reader.lines() {
        let line = line?;
        line_count += 1;

        // Progress reporting every 1 million lines
        if line_count % 1_000_000 == 0 {
            logger.log_and_progress(&format!(
                "  Processed {} lines ({} valid, {} errors)...",
                line_count, processed_count, error_count
            ))?;
        }

        if line.trim().is_empty() {
            continue;
        }

        match parse_bamreadcount_line(&line) {
            Ok((chr, position, ref_base, depth, base_a, base_c, base_g, base_t, ins, del)) => {
                // Calculate mutation count
                let mutation_count =
                    calculate_mutation_count(ref_base, base_a, base_c, base_g, base_t);

                // Write CSV line
                writeln!(
                    writer,
                    "{},{},{},{},{},{},{},{},{},{},{},{},{}",
                    chr,
                    args.strand,
                    position,
                    ref_base,
                    mutation_count,
                    0, // pipe_truncation_count (not available in bam-readcount)
                    depth,
                    ins,
                    del,
                    base_a,
                    base_c,
                    base_g,
                    base_t
                )?;

                processed_count += 1;
            }
            Err(e) => {
                error_count += 1;
                if error_count <= 10 {
                    // Log first 10 errors
                    logger.log(&format!("Warning: Failed to parse line {}: {}", line_count, e))?;
                }
                // Continue processing despite errors
            }
        }
    }

    logger.finish_progress()?;
    logger.log(&format!(
        "Processing completed: total {} lines, {} valid, {} errors",
        line_count, processed_count, error_count
    ))?;

    Ok(())
}

/// Merge two rf-rctools files (mutation and stop) into a single pileup
/// If a file is missing, corresponding columns will be set to 0
fn merge_rf_rctools_to_pileup(
    mutation_path: Option<&String>,
    stop_path: Option<&String>,
    writer: &mut BufWriter<File>,
    args: &ConvertArgs,
    logger: &mut crate::Logger,
) -> Result<(), Box<dyn Error>> {
    use std::collections::HashMap;
    
    logger.log("Merging mutation and stop rf-rctools files...")?;
    
    // Helper function to read a file and return a HashMap of (chr, position) -> (ref_base, mutation_count, rt_stops, coverage)
    fn read_rf_file(
        path: Option<&String>,
        is_mutation: bool,
        logger: &mut crate::Logger,
    ) -> Result<HashMap<(String, u64), (char, u64, u64, u64)>, Box<dyn Error>> {
        let mut data = HashMap::new();
        
        if let Some(file_path) = path {
            if file_path.trim().is_empty() || !Path::new(file_path).exists() {
                logger.log(&format!("Warning: {} file not found or empty, will use zero values", 
                    if is_mutation { "mutation" } else { "stop" }))?;
                return Ok(data);
            }
            
            logger.log(&format!("Reading {} file: {}", 
                if is_mutation { "mutation" } else { "stop" }, file_path))?;
            
            let file = File::open(file_path)?;
            let reader = BufReader::new(file);
            let mut line_iter = reader.lines();
            
            // First line is transcript name
            let mut chr = match line_iter.next() {
                Some(Ok(line)) => line.trim().to_string(),
                Some(Err(e)) => return Err(format!("Failed to read transcript name: {}", e).into()),
                None => return Err("Input file is empty".into()),
            };
            
            let mut position = 1u64;
            let mut line_count = 0u64;
            
            for line in line_iter {
                let line = line?;
                line_count += 1;
                
                if line_count % 100_000 == 0 {
                    logger.log_and_progress(&format!("  Processed {} lines from {} file...", 
                        line_count, if is_mutation { "mutation" } else { "stop" }))?;
                }
                
                let trimmed = line.trim();
                if trimmed.is_empty() {
                    continue;
                }
                
                // Check if this is a new gene ID
                let fields: Vec<&str> = trimmed.split('\t').collect();
                if fields.len() == 1 && !trimmed.chars().next().map_or(false, |c| matches!(c.to_ascii_uppercase(), 'A' | 'C' | 'G' | 'T' | 'U')) {
                    let new_chr = trimmed.to_string();
                    if new_chr != chr {
                        chr = new_chr;
                        position = 1;
                    }
                    continue;
                }
                
                match parse_rf_rctools_line(&line, &chr, position, is_mutation) {
                    Ok((ref_base, mutation_count, rt_stops, coverage)) => {
                        data.insert((chr.clone(), position), (ref_base, mutation_count, rt_stops, coverage));
                        position += 1;
                    }
                    Err(_) => {
                        position += 1; // Continue on error
                    }
                }
            }
            
            logger.log(&format!("Read {} positions from {} file", data.len(), 
                if is_mutation { "mutation" } else { "stop" }))?;
        }
        
        Ok(data)
    }
    
    // Read both files
    let mutation_data = read_rf_file(mutation_path, true, logger)?;
    let stop_data = read_rf_file(stop_path, false, logger)?;
    
    // Get all unique (chr, position) keys
    let mut all_keys: Vec<(String, u64)> = mutation_data.keys().cloned().collect();
    for key in stop_data.keys() {
        if !all_keys.contains(key) {
            all_keys.push(key.clone());
        }
    }
    
    // Sort by chr and position
    all_keys.sort_by(|a, b| {
        a.0.cmp(&b.0).then_with(|| a.1.cmp(&b.1))
    });
    
    logger.log(&format!("Merging {} unique positions...", all_keys.len()))?;
    
    // Merge and write
    let mut processed_count = 0u64;
    let mut current_chr = String::new();
    
    for (chr, position) in all_keys {
        if chr != current_chr {
            current_chr = chr.clone();
            logger.log(&format!("Processing transcript/ChrID: {}", chr))?;
        }
        
        // Get data from both files, defaulting to 0 if missing
        let (mut ref_base, mut mutation_count, mut rt_stops, mut coverage) = 
            mutation_data.get(&(chr.clone(), position))
                .copied()
                .unwrap_or(('N', 0, 0, 0));
        
        // Merge stop data
        if let Some((stop_ref_base, _, stop_rt_stops, stop_coverage)) = 
            stop_data.get(&(chr.clone(), position)) {
            // Use ref_base from mutation if available, otherwise from stop
            if ref_base == 'N' {
                ref_base = *stop_ref_base;
            }
            rt_stops = *stop_rt_stops;
            // Use maximum coverage
            coverage = coverage.max(*stop_coverage);
        }
        
        // If we still don't have a ref_base, try to get it from stop data
        if ref_base == 'N' {
            if let Some((stop_ref_base, _, _, _)) = stop_data.get(&(chr.clone(), position)) {
                ref_base = *stop_ref_base;
            }
        }
        
        // Set base counts based on ref_base and coverage (approximation)
        let (base_a, base_c, base_g, base_t) = match ref_base.to_ascii_uppercase() {
            'A' => (coverage, 0, 0, 0),
            'C' => (0, coverage, 0, 0),
            'G' => (0, 0, coverage, 0),
            'T' | 'U' => (0, 0, 0, coverage),
            _ => (0, 0, 0, 0),
        };
        
        // Write CSV line
        writeln!(
            writer,
            "{},{},{},{},{},{},{},{},{},{},{},{},{}",
            chr,
            args.strand,
            position,
            ref_base,
            mutation_count, // rf_mutation_Count
            rt_stops, // pipe_truncation_count
            coverage, // depth
            0, // rf_mutation_ins (not available)
            0, // rf_mutation_del (not available)
            base_a,
            base_c,
            base_g,
            base_t
        )?;
        
        processed_count += 1;
        
        if processed_count % 100_000 == 0 {
            logger.log_and_progress(&format!("  Merged {} positions...", processed_count))?;
        }
    }
    
    logger.finish_progress()?;
    logger.log(&format!("Merge completed: {} positions processed", processed_count))?;
    
    Ok(())
}

/// Convert rf-rctools format to pileup (internal function)
fn convert_rf_rctools_to_pileup(
    reader: BufReader<File>,
    writer: &mut BufWriter<File>,
    args: &ConvertArgs,
    logger: &mut crate::Logger,
) -> Result<(), Box<dyn Error>> {
    // First line is transcript name (ChrID)
    let mut line_iter = reader.lines();
    let mut chr = match line_iter.next() {
        Some(Ok(line)) => line.trim().to_string(),
        Some(Err(e)) => return Err(format!("Failed to read transcript name: {}", e).into()),
        None => return Err("Input file is empty".into()),
    };

    logger.log(&format!("Transcript/ChrID: {}", chr))?;

    // Process data lines
    let mut line_count = 0u64;
    let mut processed_count = 0u64;
    let mut error_count = 0u64;
    let mut position = 1u64; // Positions start from 1

    for line in line_iter {
        let line = line?;
        line_count += 1;

        // Progress reporting every 100k lines
        if line_count % 100_000 == 0 {
            logger.log_and_progress(&format!(
                "  Processed {} lines ({} valid, {} errors)...",
                line_count, processed_count, error_count
            ))?;
        }

        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }

        // Check if this line is a new gene ID (not tab-separated, looks like a gene ID)
        // Gene IDs typically contain underscores and don't start with a nucleotide
        let fields: Vec<&str> = trimmed.split('\t').collect();
        if fields.len() == 1 && !trimmed.chars().next().map_or(false, |c| matches!(c.to_ascii_uppercase(), 'A' | 'C' | 'G' | 'T' | 'U')) {
            // This looks like a new gene ID line
            let new_chr = trimmed.to_string();
            if new_chr != chr {
                logger.log(&format!("New transcript/ChrID detected: {} (previous: {})", new_chr, chr))?;
                chr = new_chr;
                position = 1; // Reset position for new gene
            }
            continue;
        }

        match parse_rf_rctools_line(&line, &chr, position, args.rf_count_mutations) {
            Ok((ref_base, mutation_count, rt_stops, coverage)) => {
                // For rf-rctools format:
                // - When rf_count_mutations=true: field[1] is mutation_count, rt_stops=0
                // - When rf_count_mutations=false: field[1] is rt_stops, mutation_count=0
                // - coverage -> depth
                // - No base counts, so we set base counts based on ref_base and coverage
                //   (This is an approximation - we don't have actual base counts)
                let (base_a, base_c, base_g, base_t) = match ref_base.to_ascii_uppercase() {
                    'A' => (coverage, 0, 0, 0),
                    'C' => (0, coverage, 0, 0),
                    'G' => (0, 0, coverage, 0),
                    'T' | 'U' => (0, 0, 0, coverage),
                    _ => (0, 0, 0, 0),
                };

                // Write CSV line
                writeln!(
                    writer,
                    "{},{},{},{},{},{},{},{},{},{},{},{},{}",
                    chr,
                    args.strand,
                    position,
                    ref_base,
                    mutation_count, // rf_mutation_Count
                    rt_stops, // pipe_truncation_count
                    coverage, // depth
                    0, // rf_mutation_ins (not available)
                    0, // rf_mutation_del (not available)
                    base_a,
                    base_c,
                    base_g,
                    base_t
                )?;

                processed_count += 1;
                position += 1;
            }
            Err(e) => {
                error_count += 1;
                if error_count <= 10 {
                    // Log first 10 errors
                    logger.log(&format!("Warning: Failed to parse line {}: {}", line_count, e))?;
                }
                // Continue processing despite errors
                position += 1; // Still increment position
            }
        }
    }

    logger.finish_progress()?;
    logger.log(&format!(
        "Processing completed: total {} lines, {} valid, {} errors",
        line_count, processed_count, error_count
    ))?;

    Ok(())
}


/// Convert rf-norm output format to ModDetector norm format
/// Input format: ChrID,Strand,Position,Base,reactivity_stop,reactivity_mutation
/// Output format: ChrID,Strand,Position,Base,reactivity_stop,reactivity_mutation (same format, but validated and standardized)
fn convert_rf_norm_to_norm(
    reader: BufReader<File>,
    writer: &mut BufWriter<File>,
    args: &ConvertArgs,
    logger: &mut crate::Logger,
) -> Result<(), Box<dyn Error>> {
    let mut line_iter = reader.lines();
    let mut line_count = 0u64;
    let mut processed_count = 0u64;
    let mut error_count = 0u64;

    // Read and validate header
    let header = match line_iter.next() {
        Some(Ok(h)) => {
            line_count += 1;
            h
        }
        Some(Err(e)) => return Err(format!("Error reading header: {}", e).into()),
        None => return Err("Error: Input file is empty".into()),
    };

    // Validate header format
    let header_fields: Vec<&str> = header.split(',').map(|s| s.trim()).collect();
    if header_fields.len() < 6 {
        return Err(format!(
            "Error: Invalid rf-norm format header. Expected at least 6 columns, got {}",
            header_fields.len()
        )
        .into());
    }

    // Expected header: ChrID,Strand,Position,Base,reactivity_stop,reactivity_mutation
    if !header_fields[0].eq_ignore_ascii_case("ChrID")
        || !header_fields[1].eq_ignore_ascii_case("Strand")
        || !header_fields[2].eq_ignore_ascii_case("Position")
        || !header_fields[3].eq_ignore_ascii_case("Base")
    {
        return Err("Error: Invalid rf-norm format header. Expected: ChrID,Strand,Position,Base,reactivity_stop,reactivity_mutation".into());
    }

    // Write header (same format)
    writeln!(writer, "{}", header)?;

    // Process data lines
    for line in line_iter {
        let line = line?;
        line_count += 1;

        // Progress reporting every 100,000 lines
        if line_count % 100_000 == 0 {
            logger.log(&format!("  Processed {} lines...", line_count))?;
        }

        if line.trim().is_empty() {
            continue;
        }

        // Parse and validate line
        let fields: Vec<&str> = line.split(',').map(|s| s.trim()).collect();
        if fields.len() < 6 {
            error_count += 1;
            if error_count <= 10 {
                logger.log(&format!(
                    "Warning: Skipping line {} (insufficient fields: {})",
                    line_count,
                    fields.len()
                ))?;
            }
            continue;
        }

        // Validate and write line
        // Fields: ChrID, Strand, Position, Base, reactivity_stop, reactivity_mutation
        let chr = fields[0];
        let strand = fields[1];
        let position = fields[2];
        let base = fields[3];
        let reactivity_stop = fields[4];
        let reactivity_mutation = fields[5];

        // Validate strand
        if strand != "+" && strand != "-" {
            error_count += 1;
            if error_count <= 10 {
                logger.log(&format!(
                    "Warning: Invalid strand '{}' at line {}, using original value",
                    strand, line_count
                ))?;
            }
            // Continue with original strand, let downstream handle it
        }

        // Validate position is numeric
        if position.parse::<u64>().is_err() {
            error_count += 1;
            if error_count <= 10 {
                logger.log(&format!(
                    "Warning: Invalid position '{}' at line {}",
                    position, line_count
                ))?;
            }
            continue;
        }

        // Validate reactivity values are numeric or NaN (just check, don't modify)
        let _stop_val: Result<f64, _> = if reactivity_stop.eq_ignore_ascii_case("NaN") {
            Ok(f64::NAN)
        } else {
            reactivity_stop.parse()
        };
        let _mut_val: Result<f64, _> = if reactivity_mutation.eq_ignore_ascii_case("NaN") {
            Ok(f64::NAN)
        } else {
            reactivity_mutation.parse()
        };

        // Write line (same format, just validated)
        writeln!(
            writer,
            "{},{},{},{},{},{}",
            chr, strand, position, base, reactivity_stop, reactivity_mutation
        )?;

        processed_count += 1;
    }

    logger.log(&format!(
        "Processing completed: total {} lines, {} valid, {} errors",
        line_count, processed_count, error_count
    ))?;

    Ok(())
}

/// Convert rf-norm XML output format to ModDetector norm format
/// Input format: XML file from rf-norm with <data><transcript><sequence> and <reactivity>
/// Output format: ChrID,Strand,Position,Base,reactivity_stop,reactivity_mutation
fn convert_rf_norm_xml_to_norm(
    input_path: &str,
    writer: &mut BufWriter<File>,
    args: &ConvertArgs,
    logger: &mut crate::Logger,
) -> Result<(), Box<dyn Error>> {
    logger.log("Parsing XML file...")?;
    
    let mut file = File::open(input_path)?;
    let mut xml_content = String::new();
    file.read_to_string(&mut xml_content)?;
    
    let mut reader = Reader::from_str(&xml_content);
    reader.config_mut().trim_text(true);
    
    let mut transcript_id = String::new();
    let mut sequence = String::new();
    let mut reactivity_text = String::new();
    let mut in_transcript = false;
    let mut in_sequence = false;
    let mut in_reactivity = false;
    
    // Parse XML
    loop {
        match reader.read_event() {
            Ok(Event::Start(e)) => {
                match e.name().as_ref() {
                    b"transcript" => {
                        in_transcript = true;
                        // Get transcript id from attribute
                        for attr in e.attributes() {
                            if let Ok(attr) = attr {
                                if attr.key.as_ref() == b"id" {
                                    transcript_id = String::from_utf8_lossy(&attr.value).to_string();
                                }
                            }
                        }
                    }
                    b"sequence" => {
                        in_sequence = true;
                        sequence.clear();
                    }
                    b"reactivity" => {
                        in_reactivity = true;
                        reactivity_text.clear();
                    }
                    _ => {}
                }
            }
            Ok(Event::Text(e)) => {
                let text = e.unescape().unwrap_or_default();
                if in_sequence {
                    sequence.push_str(&text);
                } else if in_reactivity {
                    reactivity_text.push_str(&text);
                }
            }
            Ok(Event::End(e)) => {
                match e.name().as_ref() {
                    b"transcript" => {
                        in_transcript = false;
                    }
                    b"sequence" => {
                        in_sequence = false;
                    }
                    b"reactivity" => {
                        in_reactivity = false;
                    }
                    _ => {}
                }
            }
            Ok(Event::Eof) => break,
            Err(e) => return Err(format!("XML parsing error: {}", e).into()),
            _ => {}
        }
    }
    
    if transcript_id.is_empty() || sequence.is_empty() || reactivity_text.is_empty() {
        return Err("Error: Could not find transcript, sequence, or reactivity in XML file".into());
    }
    
    // Clean sequence and reactivity
    sequence = sequence.replace('\n', "").replace('\t', "").replace(' ', "");
    reactivity_text = reactivity_text.replace('\n', "").replace('\t', "");
    
    // Parse reactivity values
    let reactivity_values: Vec<&str> = reactivity_text
        .split(',')
        .map(|s| s.trim())
        .filter(|s| !s.is_empty())
        .collect();
    
    if sequence.len() != reactivity_values.len() {
        return Err(format!(
            "Error: Sequence length ({}) does not match reactivity values count ({})",
            sequence.len(),
            reactivity_values.len()
        )
        .into());
    }
    
    // Determine strand (default to '+')
    let strand = &args.strand;
    
    // Write data lines: --rf-norm-signal stop → reactivity_stop column; mut → reactivity_mutation column
    let is_stop = args.rf_norm_signal.eq_ignore_ascii_case("stop");
    let mut processed_count = 0u64;
    for (pos, (base, react_val)) in sequence.chars().zip(reactivity_values.iter()).enumerate() {
        let position = (pos + 1) as u64;
        if is_stop {
            writeln!(writer, "{},{},{},{},{},NaN", transcript_id, strand, position, base, react_val)?;
        } else {
            writeln!(writer, "{},{},{},{},NaN,{}", transcript_id, strand, position, base, react_val)?;
        }
        processed_count += 1;
    }
    
    logger.log(&format!(
        "Processing completed: {} positions processed",
        processed_count
    ))?;
    
    Ok(())
}

/// Load all reference sequences from FASTA file with their names and lengths
fn load_all_reference_sequences(fasta_path: &str) -> Result<Vec<(String, Vec<char>)>, Box<dyn Error>> {
    use std::collections::HashMap;
    
    let file = File::open(fasta_path)?;
    let reader = BufReader::new(file);
    let mut sequences: Vec<(String, Vec<char>)> = Vec::new();
    let mut current_seq_name = String::new();
    let mut current_sequence = String::new();
    
    for line in reader.lines() {
        let line = line?;
        if line.starts_with('>') {
            // Save previous sequence
            if !current_seq_name.is_empty() {
                sequences.push((current_seq_name.clone(), current_sequence.chars().collect()));
            }
            // Start new sequence
            current_seq_name = line[1..].split_whitespace().next().unwrap_or("").to_string();
            current_sequence.clear();
        } else {
            current_sequence.push_str(&line.trim());
        }
    }
    
    // Save last sequence
    if !current_seq_name.is_empty() {
        sequences.push((current_seq_name.clone(), current_sequence.chars().collect()));
    }
    
    Ok(sequences)
}

/// Convert shapemapper2 mutations.txt format to pileup (internal function)
/// Supports multi-sequence reference files by processing all sequences in order
fn convert_shapemapper2_to_pileup(
    reader: BufReader<File>,
    writer: &mut BufWriter<File>,
    args: &ConvertArgs,
    logger: &mut crate::Logger,
) -> Result<(), Box<dyn Error>> {
    // Check if reference FASTA is provided
    let ref_fasta = args.ref_fasta.as_ref()
        .ok_or("Error: --ref-fasta is required for shapemapper2 format")?;
    
    // Load all reference sequences
    logger.log("Loading reference sequences...")?;
    let ref_sequences = load_all_reference_sequences(ref_fasta)?;
    
    if ref_sequences.is_empty() {
        return Err("No sequences found in reference file".into());
    }
    
    logger.log(&format!("Loaded {} reference sequence(s)", ref_sequences.len()))?;
    for (name, seq) in &ref_sequences {
        logger.log(&format!("  {}: {} bases", name, seq.len()))?;
    }
    
    // If chr_name is specified, use only that sequence; otherwise use all sequences
    let sequences_to_process: Vec<(String, Vec<char>)> = if let Some(chr_name) = args.chr_name.as_ref() {
        ref_sequences.into_iter()
            .filter(|(name, _)| name == chr_name)
            .collect()
    } else {
        ref_sequences
    };
    
    if sequences_to_process.is_empty() {
        return Err(format!(
            "No matching sequences found. Requested: {:?}, Available: {}",
            args.chr_name,
            if args.chr_name.is_some() { "see above" } else { "all" }
        ).into());
    }
    
    logger.log(&format!("Processing {} sequence(s)...", sequences_to_process.len()))?;
    
    // Process file line by line
    let mut line_iter = reader.lines();
    let mut line_count = 0u64;
    let mut processed_count = 0u64;
    let mut error_count = 0u64;
    
    // Skip header line
    if let Some(Ok(_header)) = line_iter.next() {
        line_count += 1;
        logger.log("Skipping header line...")?;
    }
    
    // Track current sequence and position
    let mut current_seq_idx = 0usize;
    let mut current_seq_offset = 0u64; // Cumulative offset for multi-sequence processing
    
    // Calculate total expected length
    let total_expected_length: u64 = sequences_to_process.iter()
        .map(|(_, seq)| seq.len() as u64)
        .sum();
    
    // Process data lines
    let mut data_lines_processed = 0u64;
    for line in line_iter {
        let line = line?;
        line_count += 1;
        
        // Progress reporting every 100k lines
        if line_count % 100_000 == 0 {
            logger.log_and_progress(&format!(
                "  Processed {} lines ({} valid, {} errors)...",
                line_count, processed_count, error_count
            ))?;
        }
        
        if line.trim().is_empty() {
            continue;
        }
        
        data_lines_processed += 1;
        
        // Check if we need to move to next sequence
        while current_seq_idx < sequences_to_process.len() {
            let (chr_name, ref_sequence) = &sequences_to_process[current_seq_idx];
            let seq_length = ref_sequence.len() as u64;
            let position_in_seq = data_lines_processed - current_seq_offset;
            
            if position_in_seq <= seq_length {
                // This line belongs to current sequence
                let position = position_in_seq; // 1-based position within sequence
                let ref_base = ref_sequence[(position - 1) as usize];
                
                match parse_shapemapper2_line(&line) {
                    Ok((_total_mutations, depth, base_a_mut, base_c_mut, base_g_mut, base_t_mut, ins, del)) => {
                        // Calculate mutation count based on reference base (non-reference bases only)
                        // This matches the logic used by ModDetector and other methods
                        let mutation_count = match ref_base.to_ascii_uppercase() {
                            'A' => base_c_mut + base_g_mut + base_t_mut,
                            'C' => base_a_mut + base_g_mut + base_t_mut,
                            'G' => base_a_mut + base_c_mut + base_t_mut,
                            'T' | 'U' => base_a_mut + base_c_mut + base_g_mut,
                            _ => 0,
                        };
                        
                        // Calculate actual base counts from depth and mutation counts
                        // shapemapper2 format only provides mutation counts, not actual base counts
                        // We need to infer base counts: ref_base_count = depth - mutations - deletions
                        // For non-reference bases, use mutation counts directly
                        let (base_a, base_c, base_g, base_t) = match ref_base.to_ascii_uppercase() {
                            'A' => {
                                // Reference is A, so base_A = depth - mutations - deletions
                                let base_a_count = depth.saturating_sub(base_c_mut + base_g_mut + base_t_mut + del);
                                (base_a_count, base_c_mut, base_g_mut, base_t_mut)
                            },
                            'C' => {
                                let base_c_count = depth.saturating_sub(base_a_mut + base_g_mut + base_t_mut + del);
                                (base_a_mut, base_c_count, base_g_mut, base_t_mut)
                            },
                            'G' => {
                                let base_g_count = depth.saturating_sub(base_a_mut + base_c_mut + base_t_mut + del);
                                (base_a_mut, base_c_mut, base_g_count, base_t_mut)
                            },
                            'T' | 'U' => {
                                let base_t_count = depth.saturating_sub(base_a_mut + base_c_mut + base_g_mut + del);
                                (base_a_mut, base_c_mut, base_g_mut, base_t_count)
                            },
                            _ => (base_a_mut, base_c_mut, base_g_mut, base_t_mut),
                        };
                        
                        // Write CSV line
                        writeln!(
                            writer,
                            "{},{},{},{},{},{},{},{},{},{},{},{},{}",
                            chr_name,
                            args.strand,
                            position,
                            ref_base,
                            mutation_count,
                            0, // pipe_truncation_count (not available in shapemapper2)
                            depth,
                            ins,
                            del,
                            base_a,
                            base_c,
                            base_g,
                            base_t
                        )?;
                        
                        processed_count += 1;
                    }
                    Err(e) => {
                        error_count += 1;
                        if error_count <= 10 {
                            logger.log(&format!("Warning: Failed to parse line {}: {}", line_count, e))?;
                        }
                        // Continue processing despite errors
                    }
                }
                break; // Processed this line, move to next
            } else {
                // Move to next sequence
                current_seq_offset += seq_length;
                current_seq_idx += 1;
                if current_seq_idx < sequences_to_process.len() {
                    logger.log(&format!("Switching to sequence: {}", sequences_to_process[current_seq_idx].0))?;
                }
            }
        }
        
        // If we've processed all sequences but there are more lines, log warning
        if current_seq_idx >= sequences_to_process.len() {
            error_count += 1;
            if error_count <= 10 {
                logger.log(&format!(
                    "Warning: Line {} exceeds total reference sequence length, skipping",
                    line_count
                ))?;
            }
        }
    }
    
    // Check if data is incomplete and warn (but don't fill missing positions)
    if data_lines_processed < total_expected_length {
        logger.log(&format!(
            "Warning: mutations.txt has {} data lines, but reference has {} positions. Missing {} positions will not be included in output.",
            data_lines_processed, total_expected_length, total_expected_length - data_lines_processed
        ))?;
        
        // Log which sequences are incomplete
        let mut remaining_positions = total_expected_length - data_lines_processed;
        let mut current_check_idx = current_seq_idx;
        let mut current_check_offset = current_seq_offset;
        
        while current_check_idx < sequences_to_process.len() && remaining_positions > 0 {
            let (chr_name, ref_sequence) = &sequences_to_process[current_check_idx];
            let seq_length = ref_sequence.len() as u64;
            let current_position_in_seq = if current_check_idx == current_seq_idx {
                data_lines_processed - current_check_offset
            } else {
                0
            };
            
            if current_position_in_seq < seq_length {
                let missing_in_seq = seq_length - current_position_in_seq;
                logger.log(&format!(
                    "  Sequence {}: missing {} positions (from position {} to {})",
                    chr_name, missing_in_seq, current_position_in_seq + 1, seq_length
                ))?;
                remaining_positions -= missing_in_seq.min(remaining_positions);
            }
            
            current_check_offset += seq_length;
            current_check_idx += 1;
        }
        
        logger.log("Note: Missing positions are not filled with zeros. Only positions with actual data are included in the output.")?;
    }
    
    logger.finish_progress()?;
    logger.log(&format!(
        "Processing completed: total {} lines, {} valid, {} errors",
        line_count, processed_count, error_count
    ))?;
    
    Ok(())
}

/// Convert shapemapper profile format to ModDetector norm format
/// Input format: Tab-separated file with columns: Nucleotide, Sequence, ..., Norm_profile, ...
/// Output format: ChrID,Strand,Position,Base,reactivity_stop,reactivity_mutation
/// Supports multi-sequence profiles: if --ref-fasta is provided, sequences are inferred from reference order
fn convert_shapemapper_profile_to_norm(
    reader: BufReader<File>,
    writer: &mut BufWriter<File>,
    args: &ConvertArgs,
    logger: &mut crate::Logger,
) -> Result<(), Box<dyn Error>> {
    logger.log("Converting shapemapper profile to norm format...")?;
    
    // Load reference sequences if provided (for multi-sequence support)
    let ref_sequences: Vec<(String, Vec<char>)> = if let Some(ref_fasta) = args.ref_fasta.as_ref() {
        logger.log("Loading reference sequences for multi-sequence support...")?;
        let sequences = load_all_reference_sequences(ref_fasta)?;
        logger.log(&format!("Loaded {} reference sequence(s)", sequences.len()))?;
        for (name, seq) in &sequences {
            logger.log(&format!("  {}: {} bases", name, seq.len()))?;
        }
        sequences
    } else {
        Vec::new()
    };
    
    let mut line_iter = reader.lines();
    let mut line_count = 0u64;
    let mut processed_count = 0u64;
    let mut error_count = 0u64;

    // Read and parse header
    let header = match line_iter.next() {
        Some(Ok(h)) => {
            line_count += 1;
            h
        }
        Some(Err(e)) => return Err(format!("Error reading header: {}", e).into()),
        None => return Err("Error: Input file is empty".into()),
    };

    // Parse header to find column indices
    let header_fields: Vec<&str> = header.split('\t').map(|s| s.trim()).collect();
    let mut header_map: std::collections::HashMap<String, usize> = std::collections::HashMap::new();
    for (i, field) in header_fields.iter().enumerate() {
        header_map.insert(field.to_lowercase(), i);
    }

    // Find required columns
    let nucleotide_idx = header_map
        .get("nucleotide")
        .copied()
        .ok_or("Error: Header missing 'Nucleotide' column")?;
    
    // Check if Sequence_name column exists (added during merge, inserted after Nucleotide)
    let sequence_name_idx = header_map.get("sequence_name").copied();
    let has_sequence_name_col = sequence_name_idx.is_some();
    
    // Adjust indices: if Sequence_name was added after Nucleotide, subsequent columns shift by 1
    // But we use header_map which is case-insensitive and finds the actual column position
    // So we don't need to adjust - header_map already accounts for the actual column positions
    
    let sequence_idx = header_map
        .get("sequence")
        .copied()
        .ok_or("Error: Header missing 'Sequence' column")?;
    let norm_profile_idx = header_map
        .get("norm_profile")
        .copied()
        .ok_or("Error: Header missing 'Norm_profile' column")?;
    
    // Determine sequence names: prefer Sequence_name column if available, otherwise use reference sequences
    let use_sequence_name_col = has_sequence_name_col;
    let sequence_names: Vec<String> = if !ref_sequences.is_empty() {
        ref_sequences.iter().map(|(name, _)| name.clone()).collect()
    } else {
        vec![args.chr_name.as_deref().unwrap_or("unknown").to_string()]
    };
    
    let strand = &args.strand;

    logger.log(&format!(
        "Found columns: Nucleotide={}, Sequence={}, Norm_profile={}",
        nucleotide_idx, sequence_idx, norm_profile_idx
    ))?;
    if use_sequence_name_col {
        logger.log(&format!("Found Sequence_name column at index {} - will use it to identify sequences", sequence_name_idx.unwrap()))?;
    } else {
        logger.log(&format!("No Sequence_name column - will infer from reference sequences: {:?}", sequence_names))?;
    }
    
    // Track current sequence and position
    let mut current_seq_idx = 0usize;
    let mut last_nucleotide = 0u64;
    let mut current_seq_name: Option<String> = None;

    // Process data lines
    for line in line_iter {
        let line = line?;
        line_count += 1;

        // Progress reporting every 100,000 lines
        if line_count % 100_000 == 0 {
            logger.log(&format!("  Processed {} lines...", line_count))?;
        }

        if line.trim().is_empty() {
            continue;
        }

        // Parse tab-separated fields
        let fields: Vec<&str> = line.split('\t').collect();
        
        // Check if we have enough fields (include sequence_name_idx if present)
        let max_idx = nucleotide_idx.max(sequence_idx).max(norm_profile_idx)
            .max(sequence_name_idx.unwrap_or(0));
        if fields.len() <= max_idx {
            error_count += 1;
            if error_count <= 10 {
                logger.log(&format!(
                    "Warning: Skipping line {} (insufficient fields: {}, need at least {})",
                    line_count,
                    fields.len(),
                    max_idx + 1
                ))?;
            }
            continue;
        }

        // Extract values
        // Note: If Sequence_name column was added during merge, nucleotide_idx and other indices may need adjustment
        // But since Sequence_name is inserted after Nucleotide, nucleotide_idx remains the same
        let nucleotide_str = fields[nucleotide_idx].trim();
        let sequence = fields[sequence_idx].trim();
        let norm_profile_str = fields[norm_profile_idx].trim();
        
        // Extract sequence name from Sequence_name column if available
        let line_seq_name = if let Some(seq_name_idx) = sequence_name_idx {
            if fields.len() > seq_name_idx {
                Some(fields[seq_name_idx].trim().to_string())
            } else {
                error_count += 1;
                if error_count <= 10 {
                    logger.log(&format!(
                        "Warning: Sequence_name column index {} exceeds field count {} at line {}",
                        seq_name_idx, fields.len(), line_count
                    ))?;
                }
                None
            }
        } else {
            None
        };

        // Parse position (nucleotide number)
        let nucleotide_pos = match nucleotide_str.parse::<u64>() {
            Ok(p) => p,
            Err(_) => {
                error_count += 1;
                if error_count <= 10 {
                    logger.log(&format!(
                        "Warning: Invalid nucleotide position '{}' at line {}",
                        nucleotide_str, line_count
                    ))?;
                }
                continue;
            }
        };

        // Determine current sequence name
        let chr_name = if use_sequence_name_col {
            // Use Sequence_name column if available
            if let Some(ref seq_name) = line_seq_name {
                // Check if sequence changed
                if current_seq_name.as_ref() != Some(seq_name) {
                    if current_seq_name.is_some() {
                        logger.log(&format!("Switching to sequence: {}", seq_name))?;
                    } else {
                        logger.log(&format!("Starting sequence: {}", seq_name))?;
                    }
                    current_seq_name = Some(seq_name.clone());
                }
                seq_name
            } else {
                // Fallback if Sequence_name column is missing
                if current_seq_name.is_none() {
                    current_seq_name = Some(sequence_names[0].clone());
                }
                current_seq_name.as_ref().unwrap()
            }
        } else {
            // Use sequence boundary detection if Sequence_name column is not available
            // Detect sequence boundary: if nucleotide position resets to 1, we're starting a new sequence
            if !ref_sequences.is_empty() && nucleotide_pos == 1 && last_nucleotide > 1 {
                // Check if we've completed the current sequence
                if current_seq_idx < ref_sequences.len() {
                    let (seq_name, seq_data) = &ref_sequences[current_seq_idx];
                    let expected_length = seq_data.len() as u64;
                    if last_nucleotide == expected_length {
                        logger.log(&format!("Completed sequence {} ({} positions)", seq_name, expected_length))?;
                        current_seq_idx += 1;
                        if current_seq_idx < ref_sequences.len() {
                            logger.log(&format!("Starting sequence {}", ref_sequences[current_seq_idx].0))?;
                        }
                    }
                }
            }
            
            if current_seq_idx < sequence_names.len() {
                &sequence_names[current_seq_idx]
            } else {
                // Fallback: use last sequence name if we've exceeded expected sequences
                sequence_names.last().unwrap()
            }
        };
        
        last_nucleotide = nucleotide_pos;

        // Get base from Sequence column (convert to uppercase)
        let base = if sequence.is_empty() {
            'N'
        } else {
            sequence.chars().next().unwrap_or('N').to_ascii_uppercase()
        };

        // Parse Norm_profile value (handle 'nan' string)
        let mutation_reactivity_str = if norm_profile_str.eq_ignore_ascii_case("nan") || norm_profile_str.is_empty() {
            "NaN".to_string()
        } else {
            match norm_profile_str.parse::<f64>() {
                Ok(val) => {
                    // Check if it's NaN
                    if val.is_nan() {
                        "NaN".to_string()
                    } else {
                        // Use the original string to preserve formatting
                        norm_profile_str.to_string()
                    }
                }
                Err(_) => {
                    error_count += 1;
                    if error_count <= 10 {
                        logger.log(&format!(
                            "Warning: Invalid Norm_profile value '{}' at line {}, using NaN",
                            norm_profile_str, line_count
                        ))?;
                    }
                    "NaN".to_string()
                }
            }
        };

        // For stop reactivity, set to NaN (shapemapper doesn't provide stop signal)
        let stop_reactivity = "NaN";

        // Write row
        writeln!(
            writer,
            "{},{},{},{},{},{}",
            chr_name, strand, nucleotide_pos, base, stop_reactivity, mutation_reactivity_str
        )?;

        processed_count += 1;
    }

    logger.log(&format!(
        "Processing completed: total {} lines, {} valid, {} errors",
        line_count, processed_count, error_count
    ))?;

    Ok(())
}

/// Convert samtools mpileup format to pileup (internal function)
/// Input format: chr\tposition\tref_base\tstrand\tdepth\tA\tC\tG\tT\tN
/// Output format: ModDetector pileup CSV
fn convert_samtools_mpileup_to_pileup(
    reader: BufReader<File>,
    writer: &mut BufWriter<File>,
    args: &ConvertArgs,
    logger: &mut crate::Logger,
) -> Result<(), Box<dyn Error>> {
    logger.log("Converting samtools mpileup format to pileup...")?;
    
    let mut line_iter = reader.lines();
    let mut line_count = 0u64;
    let mut processed_count = 0u64;
    let mut error_count = 0u64;
    
    // Skip header line if present
    let mut first_line_processed = false;
    
    // Process data lines
    for line in line_iter {
        let line = line?;
        line_count += 1;
        
        // Check if first line is header
        if !first_line_processed {
            first_line_processed = true;
            if line.trim().starts_with("chr\t") && line.contains("position\tref_base\tstrand") {
                logger.log("Skipping header line...")?;
                continue;
            }
        }
        
        // Progress reporting every 100k lines
        if line_count % 100_000 == 0 {
            logger.log_and_progress(&format!(
                "  Processed {} lines ({} valid, {} errors)...",
                line_count, processed_count, error_count
            ))?;
        }
        
        if line.trim().is_empty() {
            continue;
        }
        
        match parse_samtools_mpileup_line(&line) {
            Ok((chr, position, ref_base, strand, depth, base_a, base_c, base_g, base_t, _base_n)) => {
                // Calculate mutation count
                let mutation_count =
                    calculate_mutation_count(ref_base, base_a, base_c, base_g, base_t);
                
                // Use strand from input file, or fall back to args.strand
                let output_strand: &str = if strand == "+" || strand == "-" {
                    &strand
                } else {
                    &args.strand
                };
                
                // Write CSV line
                writeln!(
                    writer,
                    "{},{},{},{},{},{},{},{},{},{},{},{},{}",
                    chr,
                    output_strand,
                    position,
                    ref_base,
                    mutation_count,
                    0, // pipe_truncation_count (not available in samtools mpileup)
                    depth,
                    0, // rf_mutation_ins (not available)
                    0, // rf_mutation_del (not available)
                    base_a,
                    base_c,
                    base_g,
                    base_t
                )?;
                
                processed_count += 1;
            }
            Err(e) => {
                error_count += 1;
                if error_count <= 10 {
                    logger.log(&format!("Warning: Failed to parse line {}: {}", line_count, e))?;
                }
                // Continue processing despite errors
            }
        }
    }
    
    logger.finish_progress()?;
    logger.log(&format!(
        "Processing completed: total {} lines, {} valid, {} errors",
        line_count, processed_count, error_count
    ))?;
    
    Ok(())
}

/// Convert bedgraph format to pileup (internal function)
/// Supports two formats:
///   - Old format (5 fields): chr\tstart\tend\tcoverage\tstrand
///     For 5' end counting (bedtools genomecov -5), coverage represents depth (total reads at that position)
///     RT stop count = coverage[i] - coverage[i+1] (for consecutive positions on the same strand)
///   - New format (6 fields): chr\tstart\tend\tstop_count\tdepth\tstrand
///     Directly provides RT stop count and depth
/// Note: start is 0-based, end is exclusive (0-based), we convert start+1 to 1-based position
fn convert_bedgraph_to_pileup(
    reader: BufReader<File>,
    writer: &mut BufWriter<File>,
    args: &ConvertArgs,
    logger: &mut crate::Logger,
) -> Result<(), Box<dyn Error>> {
    logger.log("Converting bedgraph format to pileup...")?;
    
    // Check if reference FASTA is provided (needed for base information)
    let ref_fasta = args.ref_fasta.as_ref();
    let ref_bases = if let Some(ref_path) = ref_fasta {
        logger.log("Loading reference sequences...")?;
        let sequences = load_all_reference_sequences(ref_path)?;
        if sequences.is_empty() {
            return Err("No sequences found in reference file".into());
        }
        use std::collections::HashMap;
        let mut bases: HashMap<String, Vec<char>> = HashMap::new();
        for (name, seq) in sequences {
            bases.insert(name, seq);
        }
        logger.log(&format!("Loaded {} reference sequence(s)", bases.len()))?;
        Some(bases)
    } else {
        logger.log("Warning: No reference FASTA provided, base information will be 'N'")?;
        None
    };
    
    // Get filter_strand parameter (if provided)
    let filter_strand = args.filter_strand.as_ref()
        .map(|s| s.trim())
        .filter(|s| !s.is_empty() && s != &"None");
    
    if let Some(strand) = filter_strand {
        logger.log(&format!("Filtering by strand: {}", strand))?;
    } else {
        logger.log("Processing all strands")?;
    }
    
    // First pass: read all bedgraph entries into memory, grouped by (chr, strand)
    logger.log("Reading bedgraph entries...")?;
    use std::collections::HashMap;
    #[derive(Debug)]
    struct BedgraphEntry {
        chr_id: String,
        strand: String,
        pos: u64,  // 1-based position
        coverage: u64,  // depth at this position
        stop_count: Option<u64>,  // RT stop count (if available in new format)
    }
    
    let mut entries_by_key: HashMap<(String, String), Vec<BedgraphEntry>> = HashMap::new();
    let mut line_count = 0u64;
    let mut error_count = 0u64;
    
    for line in reader.lines() {
        let line = line?;
        line_count += 1;
        
        if line_count % 100_000 == 0 {
            logger.log_and_progress(&format!("  Read {} lines...", line_count))?;
        }
        
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        
        // Parse bedgraph format
        // Old format (5 fields): chr\tstart\tend\tcoverage\tstrand
        // New format (6 fields): chr\tstart\tend\tstop_count\tdepth\tstrand
        let fields: Vec<&str> = trimmed.split('\t').collect();
        
        // Parse bedgraph format
        // New format (6 fields): chr\tstart\tend\tstop_count\tdepth\tstrand
        // Old format (5 fields): chr\tstart\tend\tcoverage\tstrand
        if fields.len() < 5 {
            error_count += 1;
            if error_count <= 10 {
                logger.log(&format!("Warning: Skipping line {} (insufficient fields): {}", line_count, trimmed.chars().take(50).collect::<String>()))?;
            }
            continue;
        }
        
        let start_result = fields[1].parse::<u64>();
        let _end_result = fields[2].parse::<u64>();
        
        let (stop_count_result, depth_result, strand_idx) = if fields.len() >= 6 {
            // New format: chr\tstart\tend\tstop_count\tdepth\tstrand
            let stop_count_result = fields[3].parse::<u64>();
            let depth_result = fields[4].parse::<u64>();
            let strand_idx = 5usize;
            (stop_count_result, depth_result, strand_idx)
        } else {
            // Old format: chr\tstart\tend\tcoverage\tstrand
            let stop_count_result: Result<u64, std::num::ParseIntError> = Ok(0u64);  // stop_count will be calculated later
            let coverage_result = fields[3].parse::<f64>();
            // Convert ParseFloatError to ParseIntError for type consistency
            let depth_result = coverage_result.map(|f| f as u64)
                .map_err(|_| "0".parse::<u64>().unwrap_err());  // Dummy error conversion
            let strand_idx = 4usize;
            (stop_count_result, depth_result, strand_idx)
        };
        
        if let (Ok(start), Ok(_end), Ok(stop_count), Ok(depth), _) = (start_result, _end_result, stop_count_result, depth_result, strand_idx) {
            let chr_id = fields[0].to_string();
            let strand = fields[strand_idx].trim().to_string();
            
            // Apply strand filter if specified
            if let Some(filter) = &filter_strand {
                if strand != *filter {
                    continue;
                }
            }
            
            // Convert 0-based start to 1-based position
            let pos = start + 1;
            
            let key = (chr_id.clone(), strand.clone());
            entries_by_key.entry(key).or_insert_with(Vec::new).push(BedgraphEntry {
                chr_id,
                strand,
                pos,
                coverage: depth,  // Store depth in coverage field
                stop_count: if fields.len() >= 6 { Some(stop_count) } else { None },  // Store stop_count if available
            });
        } else {
            error_count += 1;
            if error_count <= 10 {
                logger.log(&format!("Warning: Failed to parse line {} (line: {})", line_count, trimmed.chars().take(50).collect::<String>()))?;
            }
        }
    }
    
    logger.finish_progress()?;
    logger.log(&format!("Read {} lines, {} errors", line_count, error_count))?;
    logger.log(&format!("Grouped into {} (chr, strand) combinations", entries_by_key.len()))?;
    
    // Second pass: sort entries by position and calculate RT stop count
    logger.log("Calculating RT stop counts and writing pileup...")?;
    let mut processed_count = 0u64;
    
    for (key, mut entries) in entries_by_key {
        // Sort by position
        entries.sort_by_key(|e| e.pos);
        
        // Process each entry
        for (idx, entry) in entries.iter().enumerate() {
            let depth = entry.coverage;
            
            // Get RT stop count
            let stop_count = if let Some(stop) = entry.stop_count {
                // New format: stop_count is directly provided
                stop
            } else {
                // Old format: calculate RT stop count from adjacent positions
                if idx + 1 < entries.len() {
                    let next_entry = &entries[idx + 1];
                    if next_entry.pos == entry.pos + 1 {
                        // Consecutive position: RT stop = current coverage - next coverage
                        if depth > next_entry.coverage {
                            depth - next_entry.coverage
                        } else {
                            // If next coverage is larger, it means no stops at current position
                            0
                        }
                    } else {
                        // Non-consecutive: all coverage at this position are stops
                        depth
                    }
                } else {
                    // Last position: all coverage at this position are stops
                    depth
                }
            };
            
            // Get reference base
            let ref_base = if let Some(ref bases_map) = &ref_bases {
                if let Some(seq) = bases_map.get(&entry.chr_id) {
                    if entry.pos > 0 && entry.pos <= seq.len() as u64 {
                        seq[(entry.pos - 1) as usize].to_ascii_uppercase()
                    } else {
                        'N'
                    }
                } else {
                    'N'
                }
            } else {
                'N'
            };
            
            // Write pileup line
            writeln!(
                writer,
                "{},{},{},{},{},{},{},{},{},{},{},{},{}",
                entry.chr_id,
                entry.strand,
                entry.pos,
                ref_base,
                0, // mutation_count
                stop_count, // pipe_truncation_count
                depth,
                0, // ins
                0, // del
                0, // base_A
                0, // base_C
                0, // base_G
                0  // base_T
            )?;
            
            processed_count += 1;
        }
    }
    
    logger.log(&format!(
        "Processing completed: {} entries written",
        processed_count
    ))?;
    
    Ok(())
}

/// Convert icSHAPE-pipe RT format to pileup (internal function)
/// Format: chr_id\tstrand\tposition\tRT_count\tBD_count
/// Comments start with @
fn convert_icshape_rt_to_pileup(
    reader: BufReader<File>,
    writer: &mut BufWriter<File>,
    args: &ConvertArgs,
    logger: &mut crate::Logger,
) -> Result<(), Box<dyn Error>> {
    // Check if reference FASTA is provided
    let ref_fasta = args.ref_fasta.as_ref()
        .ok_or("Error: --ref-fasta is required for icSHAPE-rt format")?;
    
    // Load all reference sequences
    logger.log("Loading reference sequences...")?;
    let ref_sequences = load_all_reference_sequences(ref_fasta)?;
    
    if ref_sequences.is_empty() {
        return Err("No sequences found in reference file".into());
    }
    
    // Create a HashMap for quick lookup
    use std::collections::HashMap;
    let mut ref_bases: HashMap<String, Vec<char>> = HashMap::new();
    for (name, seq) in ref_sequences {
        ref_bases.insert(name, seq);
    }
    
    logger.log(&format!("Loaded {} reference sequence(s)", ref_bases.len()))?;
    for (name, seq) in &ref_bases {
        logger.log(&format!("  {}: {} bases", name, seq.len()))?;
    }
    
    // Get filter_strand parameter
    let filter_strand = args.filter_strand.as_ref()
        .map(|s| s.trim())
        .filter(|s| !s.is_empty() && s != &"None");
    
    if let Some(strand) = filter_strand {
        logger.log(&format!("Filtering by strand: {}", strand))?;
    } else {
        logger.log("Processing all strands")?;
    }
    
    // Process file line by line
    let mut line_count = 0u64;
    let mut processed_count = 0u64;
    let mut skipped_count = 0u64;
    let mut error_count = 0u64;
    
    for line in reader.lines() {
        let line = line?;
        line_count += 1;
        
        // Progress reporting every 100k lines
        if line_count % 100_000 == 0 {
            logger.log_and_progress(&format!(
                "  Processed {} lines ({} valid, {} skipped, {} errors)...",
                line_count, processed_count, skipped_count, error_count
            ))?;
        }
        
        let trimmed = line.trim();
        
        // Skip comment lines (starting with @) and empty lines
        if trimmed.is_empty() || trimmed.starts_with('@') {
            continue;
        }
        
        // Parse RT file format: chr_id\tstrand\tposition\tRT_count\tBD_count
        let fields: Vec<&str> = trimmed.split('\t').collect();
        if fields.len() < 5 {
            error_count += 1;
            if error_count <= 10 {
                logger.log(&format!("Warning: Skipping line {} (insufficient fields): {}", line_count, trimmed.chars().take(50).collect::<String>()))?;
            }
            continue;
        }
        
        let chr_id = fields[0].to_string();
        let strand = fields[1].to_string();
        
        // Parse numeric fields
        let pos = match fields[2].parse::<u64>() {
            Ok(p) => p,
            Err(e) => {
                error_count += 1;
                if error_count <= 10 {
                    logger.log(&format!("Warning: Failed to parse position on line {}: {} (line: {})", line_count, e, trimmed.chars().take(50).collect::<String>()))?;
                }
                continue;
            }
        };
        
        let rt_count_f = match fields[3].parse::<f64>() {
            Ok(c) => c,
            Err(e) => {
                error_count += 1;
                if error_count <= 10 {
                    logger.log(&format!("Warning: Failed to parse RT_count on line {}: {} (line: {})", line_count, e, trimmed.chars().take(50).collect::<String>()))?;
                }
                continue;
            }
        };
        
        let bd_count_f = match fields[4].parse::<f64>() {
            Ok(c) => c,
            Err(e) => {
                error_count += 1;
                if error_count <= 10 {
                    logger.log(&format!("Warning: Failed to parse BD_count on line {}: {} (line: {})", line_count, e, trimmed.chars().take(50).collect::<String>()))?;
                }
                continue;
            }
        };
        
        // Apply strand filter if specified
        if let Some(filter) = filter_strand {
            if strand != filter {
                skipped_count += 1;
                continue;
            }
        }
        
        let rt_count = rt_count_f as u64;
        let bd_count = bd_count_f as u64;
        
        // Get reference base
        let ref_base = if let Some(seq) = ref_bases.get(&chr_id) {
            if pos > 0 && pos <= seq.len() as u64 {
                seq[(pos - 1) as usize].to_ascii_uppercase()
            } else {
                'N'
            }
        } else {
            'N'
        };
        
        // depth = RT_count + BD_count (total depth)
        let depth = if bd_count > 0 {
            rt_count + bd_count
        } else {
            rt_count
        };
        
        // Write pileup line
        // For stop signal, mutation_count=0, pipe_truncation_count=RT_count
        writeln!(
            writer,
            "{},{},{},{},{},{},{},{},{},{},{},{},{}",
            chr_id,
            strand,
            pos,
            ref_base,
            0, // mutation_count
            rt_count, // pipe_truncation_count
            depth,
            0, // ins
            0, // del
            0, // base_A
            0, // base_C
            0, // base_G
            0  // base_T
        )?;
        
        processed_count += 1;
    }
    
    logger.finish_progress()?;
    logger.log(&format!(
        "Processing completed: total {} lines, {} valid, {} skipped, {} errors",
        line_count, processed_count, skipped_count, error_count
    ))?;
    
    Ok(())
}
