use std::collections::HashMap;
use std::error::Error;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use rayon::prelude::*;

/// Gene region information from GTF file
#[derive(Debug, Clone)]
pub struct GeneRegion {
    pub seqname: String,
    pub start: u64,
    pub end: u64,
    pub strand: String,
    pub gene_id: String,
    pub gene_name: Option<String>,
}

/// Parse GTF file and extract gene regions matching the target gene name/ID
pub fn parse_gtf_for_gene(
    gtf_file: &str,
    target_gene: &str,
) -> Result<Vec<GeneRegion>, Box<dyn Error>> {
    let file = File::open(gtf_file)?;
    let reader = BufReader::new(file);
    let mut gene_regions: HashMap<String, GeneRegion> = HashMap::new();

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 9 {
            continue;
        }

        let feature = fields[2];
        // Only process gene and transcript features
        if feature != "gene" && feature != "transcript" && feature != "exon" {
            continue;
        }

        let seqname = fields[0].to_string();
        let start: u64 = fields[3].parse().unwrap_or(0);
        let end: u64 = fields[4].parse().unwrap_or(0);
        let strand = fields[6].to_string();
        let attributes = fields[8];

        // Parse attributes to find gene_id and gene_name
        let mut gene_id: Option<String> = None;
        let mut gene_name: Option<String> = None;

        for attr in attributes.split(';') {
            let attr = attr.trim();
            if attr.starts_with("gene_id \"") {
                if let Some(start_quote) = attr.find('"') {
                    if let Some(end_quote) = attr[start_quote + 1..].find('"') {
                        gene_id = Some(attr[start_quote + 1..start_quote + 1 + end_quote].to_string());
                    }
                }
            } else if attr.starts_with("gene_name \"") {
                if let Some(start_quote) = attr.find('"') {
                    if let Some(end_quote) = attr[start_quote + 1..].find('"') {
                        gene_name = Some(attr[start_quote + 1..start_quote + 1 + end_quote].to_string());
                    }
                }
            }
        }

        // Check if this gene matches our target
        let matches = if let Some(ref gid) = gene_id {
            gid.contains(target_gene) || gid == target_gene
        } else {
            false
        } || if let Some(ref gname) = gene_name {
            gname.contains(target_gene) || gname == target_gene
        } else {
            false
        };

        if matches {
            if let Some(ref gid) = gene_id {
                // For gene features, store the full gene region
                if feature == "gene" {
                    gene_regions.insert(
                        gid.clone(),
                        GeneRegion {
                            seqname: seqname.clone(),
                            start,
                            end,
                            strand: strand.clone(),
                            gene_id: gid.clone(),
                            gene_name: gene_name.clone(),
                        },
                    );
                } else if feature == "transcript" || feature == "exon" {
                    // For transcript/exon features, update the gene region to include all exons
                    let entry = gene_regions.entry(gid.clone()).or_insert_with(|| {
                        GeneRegion {
                            seqname: seqname.clone(),
                            start,
                            end,
                            strand: strand.clone(),
                            gene_id: gid.clone(),
                            gene_name: gene_name.clone(),
                        }
                    });
                    // Expand the region to include this exon
                    entry.start = entry.start.min(start);
                    entry.end = entry.end.max(end);
                }
            }
        }
    }

    let regions: Vec<GeneRegion> = gene_regions.into_values().collect();
    
    if regions.is_empty() {
        return Err(format!(
            "No matching gene found: '{}'. Please check the gene_id or gene_name fields in the GTF file.",
            target_gene
        ).into());
    }

    Ok(regions)
}

/// Extract count results for a single gene region and calculate average depth
pub fn extract_single_gene_region(
    count_csv: &str,
    region: &GeneRegion,
    output_csv: &str,
    use_relative_position: bool,
) -> Result<(usize, f64), Box<dyn Error>> {
    let file = File::open(count_csv)?;
    let reader = BufReader::new(file);
    let mut output_file = BufWriter::new(File::create(output_csv)?);

    let mut line_count = 0;
    let mut extracted_count = 0;
    let mut total_depth = 0u64;
    let mut depth_count = 0u64;
    let mut header_written = false;

    for line in reader.lines() {
        let line = line?;
        line_count += 1;

        // Write header
        if line_count == 1 {
            writeln!(output_file, "{}", line)?;
            header_written = true;
            continue;
        }

        // Parse CSV line
        let fields: Vec<&str> = line.split(',').collect();
        if fields.len() < 7 {
            continue;
        }

        let chr_id = fields[0].trim_matches('"');
        let strand = fields[1].trim();
        let pos_str = fields[2].trim();
        
        // Check if this position is in the gene region
        if chr_id == region.seqname && strand == region.strand {
            if let Ok(pos) = pos_str.parse::<u64>() {
                if pos >= region.start && pos <= region.end {
                    // Parse depth (column 7, index 6)
                    if let Ok(depth) = fields[6].trim().parse::<u64>() {
                        total_depth += depth;
                        depth_count += 1;
                    }
                    
                    // Convert to relative position if requested
                    if use_relative_position {
                        let relative_pos = pos - region.start + 1; // 1-based relative position
                        // Reconstruct line with relative position
                        let mut new_fields: Vec<String> = fields.iter().map(|s| s.to_string()).collect();
                        new_fields[2] = relative_pos.to_string();
                        writeln!(output_file, "{}", new_fields.join(","))?;
                    } else {
                        writeln!(output_file, "{}", line)?;
                    }
                    extracted_count += 1;
                }
            }
        }
    }

    if !header_written {
        return Err("Count CSV file format error: missing header".into());
    }

    let avg_depth = if depth_count > 0 {
        total_depth as f64 / depth_count as f64
    } else {
        0.0
    };

    Ok((extracted_count, avg_depth))
}

/// Extract count results for multiple gene regions, output separate files for each
pub fn extract_gene_regions_separate(
    count_csv: &str,
    gene_regions: &[GeneRegion],
    output_prefix: &str,
    threads: usize,
    use_relative_position: bool,
) -> Result<Vec<(String, usize, f64)>, Box<dyn Error>> {
    // Use parallel processing if multiple threads requested
    if threads > 1 && gene_regions.len() > 1 {
        let results: Result<Vec<_>, Box<dyn Error + Send + Sync>> = gene_regions
            .par_iter()
            .enumerate()
            .map(|(idx, region)| -> Result<Option<(String, usize, f64)>, Box<dyn Error + Send + Sync>> {
                let gene_display_name = region.gene_name.as_ref()
                    .unwrap_or(&region.gene_id)
                    .replace('/', "_")
                    .replace(' ', "_")
                    .replace(':', "_")
                    .replace('|', "_")
                    .replace('\\', "_");
                
                let temp_output = format!("{}_temp_{}.csv", output_prefix, idx);
                
                match extract_single_gene_region(count_csv, region, &temp_output, use_relative_position) {
                    Ok((count, avg_depth)) => {
                        if count > 0 {
                            let final_output = format!("{}_{}_depth{:.1}.csv", 
                                output_prefix, 
                                gene_display_name,
                                avg_depth
                            );
                            
                            std::fs::rename(&temp_output, &final_output)
                                .map_err(|e| Box::new(std::io::Error::from(e)) as Box<dyn Error + Send + Sync>)?;
                            Ok(Some((final_output, count, avg_depth)))
                        } else {
                            let _ = std::fs::remove_file(&temp_output);
                            Ok(None)
                        }
                    }
                    Err(e) => Err(format!("{}", e).into()),
                }
            })
            .collect::<Result<Vec<_>, Box<dyn Error + Send + Sync>>>();
        
        let mut final_results = Vec::new();
        let collected_results = results.map_err(|e| format!("{}", e))?;
        for (idx, result) in collected_results.into_iter().enumerate() {
            if let Some((filename, count, avg_depth)) = result {
                let gene_display_name = gene_regions[idx].gene_name.as_ref()
                    .unwrap_or(&gene_regions[idx].gene_id)
                    .replace('/', "_")
                    .replace(' ', "_")
                    .replace(':', "_")
                    .replace('|', "_")
                    .replace('\\', "_");
                
                println!(
                    "  ✓ Gene {}: extracted {} rows, average depth {:.2}, output file: {}",
                    gene_display_name, count, avg_depth, filename
                );
                final_results.push((filename, count, avg_depth));
            } else {
                let gene_display_name = gene_regions[idx].gene_name.as_ref()
                    .unwrap_or(&gene_regions[idx].gene_id)
                    .replace('/', "_")
                    .replace(' ', "_")
                    .replace(':', "_")
                    .replace('|', "_")
                    .replace('\\', "_");
                println!(
                    "  ✗ Gene {}: no data found, skipped",
                    gene_display_name
                );
            }
        }
        Ok(final_results)
    } else {
        // Sequential processing for single thread or single region
        let mut results = Vec::new();

        for (idx, region) in gene_regions.iter().enumerate() {
            // Generate output filename with gene name and average depth placeholder
            // We'll update it after calculating the depth
            let gene_display_name = region.gene_name.as_ref()
                .unwrap_or(&region.gene_id)
                .replace('/', "_")
                .replace(' ', "_")
                .replace(':', "_")
                .replace('|', "_")
                .replace('\\', "_");
            
            let temp_output = format!("{}_temp_{}.csv", output_prefix, idx);
            
            // Extract data for this region
            let (count, avg_depth) = extract_single_gene_region(count_csv, region, &temp_output, use_relative_position)?;
            
            if count > 0 {
                // Generate final filename with average depth
                // Format: {prefix}_{gene_name}_depth{avg_depth}.csv
                let final_output = format!("{}_{}_depth{:.1}.csv", 
                    output_prefix, 
                    gene_display_name,
                    avg_depth
                );
                
                // Rename temp file to final filename
                std::fs::rename(&temp_output, &final_output)?;
                
                results.push((final_output.clone(), count, avg_depth));
                
                println!(
                    "  ✓ Gene {}: extracted {} rows, average depth {:.2}, output file: {}",
                    gene_display_name, count, avg_depth, final_output
                );
            } else {
                // Remove empty temp file
                let _ = std::fs::remove_file(&temp_output);
                println!(
                    "  ✗ Gene {}: no data found, skipped",
                    gene_display_name
                );
            }
        }

        Ok(results)
    }
}

/// Main extraction function - extracts each gene region to separate files
pub fn extract_gene_from_count(
    count_csv: &str,
    gtf_file: &str,
    target_gene: &str,
    output_prefix: &str,
    threads: usize,
    use_relative_position: bool,
) -> Result<(), Box<dyn Error>> {
    println!("Parsing GTF file: {}", gtf_file);
    let gene_regions = parse_gtf_for_gene(gtf_file, target_gene)?;
    
    println!("Found {} matching gene regions:", gene_regions.len());
    for region in &gene_regions {
        println!(
            "  Gene ID: {}, chromosome: {}, position: {}-{}, strand: {}",
            region.gene_id,
            region.seqname,
            region.start,
            region.end,
            region.strand
        );
        if let Some(ref name) = region.gene_name {
            println!("    Gene name: {}", name);
        }
    }

    println!("Extracting data from count results: {}", count_csv);
    if threads > 1 {
        println!("Using {} threads for parallel processing...", threads);
    }
    println!("Generating separate output files for each gene region...");
    if use_relative_position {
        println!("Using relative positions (1-based, relative to gene start position)");
    }
    
    let results = extract_gene_regions_separate(count_csv, &gene_regions, output_prefix, threads, use_relative_position)?;
    
    println!("\nExtraction completed! Generated {} files:", results.len());
    for (filename, count, avg_depth) in &results {
        println!("  {}: {} rows, average depth {:.2}", filename, count, avg_depth);
    }
    
    Ok(())
}

/// Parse bam-readcount format line and extract base counts
/// Format: chr position reference_base depth base:count:...
/// Returns: (chr, position, reference_base, depth, base_a, base_c, base_g, base_t, ins, del)
fn parse_bamreadcount_line(line: &str) -> Result<(String, u64, char, u64, u64, u64, u64, u64, u64, u64), Box<dyn Error>> {
    let fields: Vec<&str> = line.split('\t').collect();
    if fields.len() < 4 {
        return Err("Invalid bam-readcount format: insufficient fields".into());
    }

    let chr = fields[0].to_string();
    let position: u64 = fields[1].parse()?;
    let reference_base = fields[2].chars().next().ok_or("Invalid reference base")?;
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

    Ok((chr, position, reference_base, depth, base_a, base_c, base_g, base_t, ins, del))
}

/// Extract gene regions from bam-readcount format file
pub fn extract_gene_from_bamreadcount(
    bamreadcount_file: &str,
    gtf_file: &str,
    target_gene: &str,
    output_prefix: &str,
    threads: usize,
    use_relative_position: bool,
) -> Result<(), Box<dyn Error>> {
    println!("Parsing GTF file: {}", gtf_file);
    let gene_regions = parse_gtf_for_gene(gtf_file, target_gene)?;
    
    println!("Found {} matching gene regions:", gene_regions.len());
    for region in &gene_regions {
        println!(
            "  Gene ID: {}, chromosome: {}, position: {}-{}, strand: {}",
            region.gene_id,
            region.seqname,
            region.start,
            region.end,
            region.strand
        );
        if let Some(ref name) = region.gene_name {
            println!("    Gene name: {}", name);
        }
    }

    println!("Extracting data from bam-readcount results: {}", bamreadcount_file);
    if threads > 1 {
        println!("Using {} threads for parallel processing...", threads);
    }
    println!("Generating separate output files for each gene region...");
    if use_relative_position {
        println!("Using relative positions (1-based, relative to gene start position)");
    }
    
    let results = extract_gene_regions_from_bamreadcount(bamreadcount_file, &gene_regions, output_prefix, threads, use_relative_position)?;
    
    println!("\nExtraction completed! Generated {} files:", results.len());
    for (filename, count, avg_depth) in &results {
        println!("  {}: {} rows, average depth {:.2}", filename, count, avg_depth);
    }
    
    Ok(())
}

/// Extract gene regions from bam-readcount format, output separate files for each
/// Uses streaming processing to avoid loading entire file into memory
fn extract_gene_regions_from_bamreadcount(
    bamreadcount_file: &str,
    gene_regions: &[GeneRegion],
    output_prefix: &str,
    _threads: usize,
    use_relative_position: bool,
) -> Result<Vec<(String, usize, f64)>, Box<dyn Error>> {
    // Use streaming processing instead of loading entire file into memory
    // This is critical for large files (e.g., 35GB) on NFS filesystems
    println!("Streaming bam-readcount file (not loading into memory)...");
    
    // Open input file with large buffer for better NFS performance
    let file = File::open(bamreadcount_file)?;
    let reader = BufReader::with_capacity(8 * 1024 * 1024, file); // 8MB buffer
    
    // Pre-open all output files and initialize statistics
    let mut output_files: Vec<Option<BufWriter<File>>> = Vec::new();
    let mut temp_outputs: Vec<String> = Vec::new();
    let mut extracted_counts: Vec<usize> = Vec::new();
    let mut total_depths: Vec<u64> = Vec::new();
    let mut depth_counts: Vec<u64> = Vec::new();
    
    for (idx, _region) in gene_regions.iter().enumerate() {
        let temp_output = format!("{}_temp_{}.csv", output_prefix, idx);
        let output_file = BufWriter::new(File::create(&temp_output)?);
        temp_outputs.push(temp_output);
        output_files.push(Some(output_file));
        extracted_counts.push(0);
        total_depths.push(0);
        depth_counts.push(0);
    }
    
    // Write headers to all output files
    for output_file in &mut output_files {
        if let Some(ref mut file) = output_file {
            writeln!(file, "ChrID,pipe_truncation_Strand,pipe_truncation_ChrPos,rf_mutation_Base,rf_mutation_Count,pipe_truncation_count,depth,rf_mutation_ins,rf_mutation_del,base_A,base_C,base_G,base_T")?;
        }
    }
    
    // Single pass: stream file line by line, filter and write to appropriate outputs
    let mut line_count = 0u64;
    for line in reader.lines() {
        let line = line?;
        line_count += 1;
        
        // Progress reporting every 10 million lines
        if line_count % 10_000_000 == 0 {
            println!("  Processed {} lines...", line_count);
        }
        
        if line.trim().is_empty() {
            continue;
        }

        match parse_bamreadcount_line(&line) {
            Ok((chr, position, ref_base, depth, base_a, base_c, base_g, base_t, ins, del)) => {
                // Check which regions this position belongs to
                for (idx, region) in gene_regions.iter().enumerate() {
                    if chr == region.seqname && position >= region.start && position <= region.end {
                        // Calculate mutation count (non-reference bases)
                        let mut mutation_count = 0u64;
                        let ref_base_upper = ref_base.to_ascii_uppercase();
                        match ref_base_upper {
                            'A' => mutation_count = base_c + base_g + base_t,
                            'C' => mutation_count = base_a + base_g + base_t,
                            'G' => mutation_count = base_a + base_c + base_t,
                            'T' | 'U' => mutation_count = base_a + base_c + base_g,
                            _ => {}
                        }

                        // Write to corresponding output file
                        if let Some(ref mut file) = output_files[idx] {
                            // Convert to relative position if requested
                            let output_position = if use_relative_position {
                                position - region.start + 1 // 1-based relative position
                            } else {
                                position
                            };
                            
                            writeln!(
                                file,
                                "{},{},{},{},{},{},{},{},{},{},{},{},{}",
                                chr,
                                &region.strand,
                                output_position,
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
                        }
                        
                        total_depths[idx] += depth;
                        depth_counts[idx] += 1;
                        extracted_counts[idx] += 1;
                    }
                }
            }
            Err(_) => {
                // Skip invalid lines but continue processing
                continue;
            }
        }
    }
    
    println!("Processing completed, total {} lines", line_count);
    
    // Close all output files and rename to final names
    // Files will be automatically flushed and closed when dropped
    let mut results = Vec::new();
    for (idx, region) in gene_regions.iter().enumerate() {
        // Flush and close the file by dropping it (take from Option and drop)
        let _ = output_files[idx].take();
        
        let gene_display_name = region.gene_name.as_ref()
            .unwrap_or(&region.gene_id)
            .replace('/', "_")
            .replace(' ', "_")
            .replace(':', "_")
            .replace('|', "_")
            .replace('\\', "_");
        
        if extracted_counts[idx] > 0 {
            let avg_depth = if depth_counts[idx] > 0 {
                total_depths[idx] as f64 / depth_counts[idx] as f64
            } else {
                0.0
            };

            let final_output = format!("{}_{}_depth{:.1}.csv", 
                output_prefix, 
                gene_display_name,
                avg_depth
            );
            
            std::fs::rename(&temp_outputs[idx], &final_output)?;
            
            results.push((final_output.clone(), extracted_counts[idx], avg_depth));
            
            println!(
                "  ✓ Gene {}: extracted {} rows, average depth {:.2}, output file: {}",
                gene_display_name, extracted_counts[idx], avg_depth, final_output
            );
        } else {
            let _ = std::fs::remove_file(&temp_outputs[idx]);
            println!(
                "  ✗ Gene {}: no data found, skipped",
                gene_display_name
            );
        }
    }
    
    Ok(results)
}

