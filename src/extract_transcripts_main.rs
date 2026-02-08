use rayon::prelude::*;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::time::Instant;

#[derive(Debug, Clone)]
struct TranscriptRecord {
    seqname: String,
    source: String,
    start: u64,
    end: u64,
    strand: String,
    transcript_id: String,
    gene_id: String,
}

#[derive(Debug)]
struct GenomeSequence {
    sequences: HashMap<String, String>,
}

impl GenomeSequence {
    fn new() -> Self {
        Self {
            sequences: HashMap::new(),
        }
    }

    fn load_from_fasta(&mut self, fasta_path: &str) -> Result<(), Box<dyn std::error::Error>> {
        let file = File::open(fasta_path)?;
        let reader = BufReader::new(file);
        let mut current_seq_name = String::new();
        let mut current_sequence = String::new();

        for line in reader.lines() {
            let line = line?;
            if line.starts_with('>') {
                // Save previous sequence
                if !current_seq_name.is_empty() {
                    self.sequences
                        .insert(current_seq_name.clone(), current_sequence.clone());
                }
                // Start new sequence
                current_seq_name = line[1..]
                    .split_whitespace()
                    .next()
                    .unwrap_or("")
                    .to_string();
                current_sequence.clear();
            } else {
                current_sequence.push_str(&line);
            }
        }

        // Save last sequence
        if !current_seq_name.is_empty() {
            self.sequences.insert(current_seq_name, current_sequence);
        }

        println!("Loaded {} chromosome sequences", self.sequences.len());
        Ok(())
    }

    fn get_sequence(&self, seqname: &str, start: u64, end: u64) -> Option<String> {
        self.sequences.get(seqname).map(|seq| {
            let start_idx = (start - 1) as usize; // Convert to 0-based
            let end_idx = end as usize;
            if start_idx < seq.len() && end_idx <= seq.len() {
                seq[start_idx..end_idx].to_string()
            } else {
                String::new()
            }
        })
    }
}

fn parse_gtf_line(line: &str) -> Option<TranscriptRecord> {
    let fields: Vec<&str> = line.split('\t').collect();
    if fields.len() < 9 {
        return None;
    }

    let seqname = fields[0];
    let source = fields[1];
    let feature = fields[2];
    let start: u64 = fields[3].parse().ok()?;
    let end: u64 = fields[4].parse().ok()?;
    let strand = fields[6];

    if feature != "exon" {
        return None;
    }

    // Parse attributes
    let attributes = fields[8];
    let mut transcript_id = String::new();
    let mut gene_id = String::new();

    for attr in attributes.split(';') {
        let attr = attr.trim();
        if attr.starts_with("transcript_id") {
            transcript_id = attr.split('"').nth(1)?.to_string();
        } else if attr.starts_with("gene_id") {
            gene_id = attr.split('"').nth(1)?.to_string();
        }
    }

    if transcript_id.is_empty() || gene_id.is_empty() {
        return None;
    }

    Some(TranscriptRecord {
        seqname: seqname.to_string(),
        source: source.to_string(),
        start,
        end,
        strand: strand.to_string(),
        transcript_id,
        gene_id,
    })
}

fn load_gene_mapping(
    mapping_path: &str,
) -> Result<HashMap<String, String>, Box<dyn std::error::Error>> {
    let mut mapping = HashMap::new();
    let file = File::open(mapping_path)?;
    let reader = BufReader::new(file);

    for line in reader.lines() {
        let line = line?;
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() >= 2 {
            mapping.insert(parts[0].to_string(), parts[1].to_string());
        }
    }

    println!("Loaded {} gene mappings", mapping.len());
    Ok(mapping)
}

fn load_gene_list(gene_list_path: &str) -> Result<Vec<String>, Box<dyn std::error::Error>> {
    let mut genes = Vec::new();
    let file = File::open(gene_list_path)?;
    let reader = BufReader::new(file);

    for line in reader.lines() {
        let line = line?;
        if !line.trim().is_empty() {
            genes.push(line.trim().to_string());
        }
    }

    println!("Loaded {} genes from list", genes.len());
    Ok(genes)
}

fn extract_transcript_sequence(
    transcript_records: &[TranscriptRecord],
    genome: &GenomeSequence,
) -> Option<String> {
    if transcript_records.is_empty() {
        return None;
    }

    let mut sorted_records = transcript_records.to_vec();
    sorted_records.sort_by_key(|r| r.start);

    let _seqname = &sorted_records[0].seqname;
    let strand = &sorted_records[0].strand;
    let mut transcript_sequence = String::new();

    for record in &sorted_records {
        if let Some(exon_seq) = genome.get_sequence(&record.seqname, record.start, record.end) {
            transcript_sequence.push_str(&exon_seq);
        } else {
            return None;
        }
    }

    // Reverse complement if negative strand
    if strand == "-" {
        transcript_sequence = reverse_complement(&transcript_sequence);
    }

    Some(transcript_sequence)
}

fn reverse_complement(sequence: &str) -> String {
    let mut result = String::with_capacity(sequence.len());
    for base in sequence.chars().rev() {
        let complement = match base {
            'A' => 'T',
            'T' => 'A',
            'C' => 'G',
            'G' => 'C',
            'N' => 'N',
            _ => 'N',
        };
        result.push(complement);
    }
    result
}

fn process_gene_batch(
    gene_batch: Vec<String>,
    gene_mapping: &HashMap<String, String>,
    transcript_records: &HashMap<String, Vec<TranscriptRecord>>,
    genome: &GenomeSequence,
    min_length: usize,
) -> Vec<(String, String)> {
    let mut results = Vec::new();

    for gene_id in gene_batch {
        if let Some(transcript_id) = gene_mapping.get(&gene_id) {
            if let Some(records) = transcript_records.get(transcript_id) {
                if let Some(sequence) = extract_transcript_sequence(records, genome) {
                    if sequence.len() >= min_length {
                        results.push((transcript_id.clone(), sequence));
                    }
                }
            }
        }
    }

    results
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().collect();
    if args.len() != 8 {
        eprintln!("Usage: {} <mapping_file> <gtf_file> <fasta_file> <gene_list> <output_prefix> <min_length> <threads>", args[0]);
        std::process::exit(1);
    }

    let mapping_file = &args[1];
    let gtf_file = &args[2];
    let fasta_file = &args[3];
    let gene_list_file = &args[4];
    let output_prefix = &args[5];
    let min_length: usize = args[6].parse()?;
    let num_threads: usize = args[7].parse()?;

    // Set number of threads for rayon
    rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build_global()
        .unwrap();

    let start_time = Instant::now();

    println!("Loading gene mapping...");
    let gene_mapping = load_gene_mapping(mapping_file)?;

    println!("Loading gene list...");
    let gene_list = load_gene_list(gene_list_file)?;

    println!("Loading genome sequences...");
    let mut genome = GenomeSequence::new();
    genome.load_from_fasta(fasta_file)?;

    println!("Loading GTF records...");
    let mut transcript_records: HashMap<String, Vec<TranscriptRecord>> = HashMap::new();
    let file = File::open(gtf_file)?;
    let reader = BufReader::new(file);

    for line in reader.lines() {
        let line = line?;
        if let Some(record) = parse_gtf_line(&line) {
            transcript_records
                .entry(record.transcript_id.clone())
                .or_insert_with(Vec::new)
                .push(record);
        }
    }

    println!("Loaded {} transcript records", transcript_records.len());

    // Filter genes that have mappings
    let target_genes: Vec<String> = gene_list
        .into_iter()
        .filter(|gene_id| gene_mapping.contains_key(gene_id))
        .collect();

    println!("Found {} target genes", target_genes.len());

    // Split genes into batches for parallel processing
    let batch_size = (target_genes.len() + num_threads - 1) / num_threads;
    let gene_batches: Vec<Vec<String>> = target_genes
        .chunks(batch_size)
        .map(|chunk| chunk.to_vec())
        .collect();

    println!(
        "Split into {} batches of ~{} genes each",
        gene_batches.len(),
        batch_size
    );

    // Parallel processing
    println!("Starting parallel processing...");
    let results: Vec<Vec<(String, String)>> = gene_batches
        .par_iter()
        .map(|batch| {
            process_gene_batch(
                batch.clone(),
                &gene_mapping,
                &transcript_records,
                &genome,
                min_length,
            )
        })
        .collect();

    // Flatten results
    let mut all_sequences = Vec::new();
    for batch_result in results {
        all_sequences.extend(batch_result);
    }

    // Write output files
    println!("Writing output files...");
    let fasta_output = format!("{}.fa", output_prefix);
    let mut fasta_file = File::create(&fasta_output)?;

    for (transcript_id, sequence) in &all_sequences {
        writeln!(fasta_file, ">{}", transcript_id)?;
        for chunk in sequence.as_bytes().chunks(80) {
            writeln!(fasta_file, "{}", String::from_utf8_lossy(chunk))?;
        }
    }

    let gtf_output = format!("{}.gtf", output_prefix);
    let mut gtf_file = File::create(&gtf_output)?;
    writeln!(gtf_file, "##gff-version 3")?;
    writeln!(gtf_file, "# Extracted longest transcript records")?;

    for (transcript_id, _) in &all_sequences {
        if let Some(records) = transcript_records.get(transcript_id) {
            for record in records {
                writeln!(
                    gtf_file,
                    "{}\t{}\texon\t{}\t{}\t.\t{}\t.\tgene_id \"{}\"; transcript_id \"{}\";",
                    record.seqname,
                    record.source,
                    record.start,
                    record.end,
                    record.strand,
                    record.gene_id,
                    record.transcript_id
                )?;
            }
        }
    }

    let elapsed = start_time.elapsed();
    println!("Completed in {:.2?}", elapsed);
    println!("Extracted {} transcripts", all_sequences.len());
    println!("FASTA file: {}", fasta_output);
    println!("GTF file: {}", gtf_output);

    Ok(())
}
