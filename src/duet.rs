use bio::io::fasta;
use clap::Args;
use rust_htslib::bam::{self, record::Cigar, Read};
use std::cmp::Ordering;
use std::collections::{HashMap, HashSet};
use std::error::Error;
use std::fmt;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::Path;
use std::time::Instant;
use rayon::prelude::*;
use plotters::prelude::*;

use crate::Logger;
use crate::progress::format_time_used;

const READ_SUPPORT_SCALE: f64 = 50.0;
const REACTIVITY_SCALE: f64 = 1.0;
const HIGH_REACTIVITY_THRESHOLD: f64 = 0.7;
const NEIGHBOR_RANGE: i32 = 2;
const NEIGHBOR_BOOST: f64 = 0.15;
const CONTIGUITY_MAX_GAP: u32 = 4;
const CONTIGUITY_MIN_LEN: usize = 3;
const READ_WEIGHT: f64 = 0.6;
// Minimum number of shared positions required to merge two clusters into the same global ensemble
const MIN_SHARED_POSITIONS: usize = 5;
// Minimum fraction of shared positions relative to the smaller cluster's total positions
const MIN_SHARED_FRACTION: f64 = 0.15;
// Minimum reactivity pattern correlation for merging clusters (based on DREEM principle:
// different structures have different mutation profiles, so clusters with dissimilar
// reactivity patterns should not be merged even if they share positions)
const MIN_REACTIVITY_CORRELATION: f64 = 0.60;

#[derive(Args, Debug)]
pub struct DuetArgs {
    /// Normalized reactivity CSV file (output from `modtector norm`)
    #[arg(short = 'i', long = "input")]
    pub input: String,
    /// BAM file containing aligned reads for co-variation analysis
    #[arg(short = 'b', long = "bam")]
    pub bam: String,
    /// Reference FASTA used for alignment
    #[arg(short = 'f', long = "fasta")]
    pub fasta: String,
    /// Output CSV file containing window-level ensemble assessments
    #[arg(short = 'o', long = "output")]
    pub output: String,
    /// Density threshold for DBSCAN clustering in standardized feature space
    #[arg(long = "epsilon", default_value_t = 0.75)]
    pub epsilon: f64,
    /// Minimum number of neighbours to form a core point in DBSCAN
    #[arg(long = "min-samples", default_value_t = 5)]
    pub min_samples: usize,
    /// Sliding window size (nt)
    #[arg(long = "window-size", default_value_t = 100)]
    pub window_size: u32,
    /// Sliding window step (nt)
    #[arg(long = "window-step", default_value_t = 50)]
    pub window_step: u32,
    /// Optional log file path
    #[arg(short = 'l', long = "log")]
    pub log: Option<String>,
    /// Optional summary CSV output (ensemble statistics). Defaults to `<output>_summary.csv`
    #[arg(long = "summary-output")]
    pub summary_output: Option<String>,
    /// Number of threads for parallel processing
    #[arg(short = 't', long = "threads", default_value_t = {
        std::thread::available_parallelism().map(|n| n.get()).unwrap_or(1)
    })]
    pub threads: usize,
    /// Minimum reactivity pattern correlation for merging clusters (0.0-1.0). 
    /// Higher values require more similar reactivity patterns to merge clusters.
    /// Based on DREEM principle: different structures have different mutation profiles.
    #[arg(long = "min-reactivity-correlation", default_value_t = 0.60)]
    pub min_reactivity_correlation: f64,
    /// Minimum fraction of shared positions relative to smaller cluster (0.0-1.0).
    /// Higher values require more overlap to merge clusters.
    #[arg(long = "min-shared-fraction", default_value_t = 0.30)]
    pub min_shared_fraction: f64,
    /// Minimum fraction of shared positions relative to larger cluster (0.0-1.0).
    /// Prevents large clusters from merging with small clusters.
    #[arg(long = "min-max-fraction", default_value_t = 0.20)]
    pub min_max_fraction: f64,
    /// Maximum signal positions for small ensemble merging (third pass).
    /// Ensembles with signal positions <= this value will be merged into larger ensembles.
    #[arg(long = "max-signal-positions-for-merge", default_value_t = 10)]
    pub max_signal_positions_for_merge: usize,
    /// Minimum reactivity correlation for merging small ensembles (0.0-1.0).
    /// Small ensembles will be merged if correlation >= this value with overlapping ensembles.
    #[arg(long = "min-small-ensemble-correlation", default_value_t = 0.50)]
    pub min_small_ensemble_correlation: f64,
}

pub fn validate_args(args: &DuetArgs) -> Result<(), Box<dyn Error>> {
    if args.input.trim().is_empty() {
        return Err("Error: Input file path cannot be empty".into());
    }
    if !Path::new(&args.input).exists() {
        return Err(format!("Error: Input file does not exist: {}", args.input).into());
    }
    if args.bam.trim().is_empty() {
        return Err("Error: BAM file path cannot be empty".into());
    }
    if !Path::new(&args.bam).exists() {
        return Err(format!("Error: BAM file does not exist: {}", args.bam).into());
    }
    if args.fasta.trim().is_empty() {
        return Err("Error: FASTA file path cannot be empty".into());
    }
    if !Path::new(&args.fasta).exists() {
        return Err(format!("Error: FASTA file does not exist: {}", args.fasta).into());
    }
    if args.output.trim().is_empty() || !args.output.ends_with(".csv") {
        return Err("Error: Output file path must end with .csv".into());
    }
    if args.epsilon <= 0.0 {
        return Err("Error: epsilon must be greater than 0".into());
    }
    if args.min_samples == 0 {
        return Err("Error: min-samples must be greater than 0".into());
    }
    if args.window_size == 0 {
        return Err("Error: window-size must be greater than 0".into());
    }
    if args.window_step == 0 {
        return Err("Error: window-step must be greater than 0".into());
    }
    if let Some(summary_path) = &args.summary_output {
        if summary_path.trim().is_empty() || !summary_path.ends_with(".csv") {
            return Err("Error: Summary output path must end with .csv".into());
        }
    }
    Ok(())
}

#[derive(Default, Clone)]
struct PositionStats {
    stop_reads: usize,
    mutation_reads: usize,
    co_stop_mut: usize,
    co_mut_stop: usize,
}

#[derive(Clone)]
struct ReactivityEntry {
    chr: String,
    strand: String,
    position: u32,
    _base: char,
    stop_reactivity: f64,
    mutation_reactivity: f64,
    stats: PositionStats,
}

#[derive(Clone)]
struct WindowReadFeature {
    stop_positions: Vec<u32>,
    mutation_positions: Vec<u32>,
}

#[derive(Clone)]
struct WindowDescriptor {
    chr: String,
    strand: String,
    start: u32,
    end: u32,
}

struct WindowPlan {
    descriptors: Vec<WindowDescriptor>,
    position_to_windows: HashMap<String, HashMap<u32, Vec<usize>>>,
}

impl WindowPlan {
    fn len(&self) -> usize {
        self.descriptors.len()
    }

    fn descriptor(&self, idx: usize) -> &WindowDescriptor {
        &self.descriptors[idx]
    }

    fn window_indices(&self, transcript: &str, pos: u32) -> Option<&Vec<usize>> {
        self.position_to_windows
            .get(transcript)
            .and_then(|mapping| mapping.get(&pos))
    }
}

#[derive(Clone, Copy, Default)]
struct ReadFeature {
    stop_count: usize,
    mutation_count: usize,
}

impl ReadFeature {
    fn co_occurrence(&self) -> bool {
        self.stop_count > 0 && self.mutation_count > 0
    }

    fn to_vec(self) -> Vec<f64> {
        vec![
            self.stop_count as f64,
            self.mutation_count as f64,
            if self.co_occurrence() { 1.0 } else { 0.0 },
        ]
    }
}

#[derive(Default)]
struct WindowReadAccumulator {
    stop_positions: HashSet<u32>,
    mutation_positions: HashSet<u32>,
}

#[derive(Default, Clone)]
struct WindowReactivityStats {
    mean_stop: f64,
    mean_mutation: f64,
}

enum WindowStatus {
    Ok,
    InsufficientReads,
    LowDensity,
}

impl fmt::Display for WindowStatus {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            WindowStatus::Ok => write!(f, "OK"),
            WindowStatus::InsufficientReads => write!(f, "InsufficientReads"),
            WindowStatus::LowDensity => write!(f, "LowDensity"),
        }
    }
}

struct ClusterSummary {
    id: usize,
    count: usize,
    occupancy: f64,
    mean_stop_count: f64,
    mean_mutation_count: f64,
    dual_signal_fraction: f64,
    is_noise: bool,
    confidence: f64,
    global_ids: Vec<usize>,
    cluster_key: Option<ClusterKey>,
}

#[derive(Clone, Copy, Hash, Eq, PartialEq)]
struct ClusterKey {
    window_index: usize,
    local_id: usize,
}

#[derive(Default, Clone)]
struct PositionContribution {
    stop_reads: usize,
    mutation_reads: usize,
}

struct ClusterDetail {
    key: ClusterKey,
    chr: String,
    strand: String,
    _window_start: u32,
    _window_end: u32,
    read_count: usize,
    dual_signal_reads: usize,
    positions: HashSet<u32>,
    position_stats: HashMap<u32, PositionContribution>,
    position_confidence: HashMap<u32, f64>,
    cluster_confidence: f64,
}

struct GlobalEnsemble {
    id: usize,
    chr: String,
    strand: String,
    window_count: usize,
    cluster_count: usize,
    read_count: usize,
    total_stop_reads: usize,
    total_mutation_reads: usize,
    dual_signal_reads: usize,
    unique_positions: usize,
    confidence: f64,
}

struct GlobalBaseRecord {
    ensemble_id: usize,
    chr: String,
    strand: String,
    position: u32,
    stop_reads: usize,
    mutation_reads: usize,
    stop_reactivity: f64,
    mutation_reactivity: f64,
    position_confidence: f64,
}

struct WindowResult {
    descriptor: WindowDescriptor,
    status: WindowStatus,
    read_count: usize,
    ensemble_count: usize,
    primary_id: Option<usize>,
    primary_occupancy: f64,
    noise_fraction: f64,
    dual_signal_fraction: f64,
    mean_stop_reactivity: f64,
    mean_mutation_reactivity: f64,
    clusters: Vec<ClusterSummary>,
}

pub fn analyze_reactivity(args: &DuetArgs, logger: &mut Logger) -> Result<(), Box<dyn Error>> {
    let start_time = Instant::now();
    
    logger.log("=== ModDetector Duet Function Log ===")?;
    logger.log(&format!("Software Version: v{}", crate::VERSION))?;
    logger.log(&format!(
        "Runtime: {}",
        chrono::Utc::now().format("%Y-%m-%d %H:%M:%S")
    ))?;
    logger.log(&format!("Reactivity Input: {}", args.input))?;
    logger.log(&format!("BAM Input: {}", args.bam))?;
    logger.log(&format!("Reference FASTA: {}", args.fasta))?;
    logger.log(&format!("Output File: {}", args.output))?;
    logger.log(&format!("epsilon: {:.3}", args.epsilon))?;
    logger.log(&format!("min_samples: {}", args.min_samples))?;
    logger.log(&format!("window_size: {}", args.window_size))?;
    logger.log(&format!("window_step: {}", args.window_step))?;
    logger.log(&format!("threads: {}", args.threads))?;

    // Set rayon thread pool size
    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()
        .unwrap_or(());

    let load_start = Instant::now();
    println!("[Duet] Loading normalized reactivity data");
    let mut positions = load_reactivity_table(&args.input)?;
    if positions.is_empty() {
        return Err("Error: No valid reactivity records were loaded".into());
    }
    logger.log(&format!("Loaded {} reactivity positions", positions.len()))?;
    let load_elapsed = load_start.elapsed();
    println!("{}", format_time_used(load_elapsed));

    let plan_start = Instant::now();
    println!("[Duet] Planning sliding windows");
    let plan = build_window_plan(&positions, args.window_size, args.window_step);
    if plan.len() == 0 {
        return Err("Error: No windows generated for provided reactivity data".into());
    }
    logger.log(&format!("Constructed {} sliding windows", plan.len()))?;

    let window_reactivity = compute_window_reactivity_stats(&plan, &positions);
    let plan_elapsed = plan_start.elapsed();
    println!("{}", format_time_used(plan_elapsed));

    let ref_start = Instant::now();
    println!("[Duet] Loading reference sequence(s)");
    let reference = load_reference_sequences(&args.fasta)?;
    if reference.is_empty() {
        return Err("Error: Reference FASTA is empty or unreadable".into());
    }
    let ref_elapsed = ref_start.elapsed();
    println!("{}", format_time_used(ref_elapsed));

    let bam_start = Instant::now();
    println!("[Duet] Scanning BAM for read-level co-variation");
    let mut window_features: Vec<Vec<ReadFeature>> = vec![Vec::new(); plan.len()];
    let mut read_features: Vec<Vec<WindowReadFeature>> = vec![Vec::new(); plan.len()];
    accumulate_covariation(
        &args.bam,
        &reference,
        &mut positions,
        &plan,
        &mut window_features,
        &mut read_features,
    )?;
    let bam_elapsed = bam_start.elapsed();
    println!("{}", format_time_used(bam_elapsed));

    let analyze_start = Instant::now();
    println!("[Duet] Aggregating window ensembles");
    let (mut window_results, cluster_details, (high_conf_positions, total_positions, avg_confidence)) = analyze_windows(
        &plan,
        &window_features,
        &read_features,
        &window_reactivity,
        args.epsilon,
        args.min_samples,
        &positions,
    );
    
    let high_conf_fraction = if total_positions > 0 {
        high_conf_positions as f64 / total_positions as f64
    } else {
        0.0
    };
    println!("[Duet]   {} positions ({:.2}%) assigned high confidence (reactivity >= {:.2})", 
             high_conf_positions, high_conf_fraction * 100.0, HIGH_REACTIVITY_THRESHOLD);
    println!("[Duet]   Average position confidence: {:.3} (total positions: {})", 
             avg_confidence, total_positions);
    let analyze_elapsed = analyze_start.elapsed();
    println!("{}", format_time_used(analyze_elapsed));

    let global_start = Instant::now();
    // First aggregation: window-level clusters -> initial ensembles
    let (first_ensembles, cluster_to_first, first_base_records) =
        aggregate_global_ensembles(
            &cluster_details, 
            &positions,
            args.min_reactivity_correlation,
            args.min_shared_fraction,
            args.min_max_fraction,
        );
    println!("[Duet] First aggregation: {} ensembles", first_ensembles.len());
    
    // Second aggregation: merge initial ensembles further based on signal positions only
    let (second_ensembles, first_to_second, second_base_records) =
        aggregate_ensembles_second_pass(
            &first_ensembles,
            &first_base_records,
            &positions,
            args.min_reactivity_correlation,
            args.min_shared_fraction,
            args.min_max_fraction,
        );
    println!("[Duet] Second aggregation: {} ensembles (from {} initial)", 
              second_ensembles.len(), first_ensembles.len());
    
    // Third aggregation: merge small fragments into larger ensembles
    let (global_ensembles, second_to_global, global_base_records) =
        aggregate_small_ensembles(
            &second_ensembles,
            &second_base_records,
            &positions,
            args.max_signal_positions_for_merge,
            args.min_small_ensemble_correlation,
        );
    println!("[Duet] Third aggregation: {} ensembles (from {} second-pass)", 
              global_ensembles.len(), second_ensembles.len());
    
    // Build cluster_to_global mapping from first_to_second and second_to_global
    let mut cluster_to_global: HashMap<ClusterKey, usize> = HashMap::new();
    for (cluster_key, first_id) in cluster_to_first {
        if let Some(&second_id) = first_to_second.get(&first_id) {
            if let Some(&global_id) = second_to_global.get(&second_id) {
                cluster_to_global.insert(cluster_key, global_id);
            }
        }
    }
    let global_elapsed = global_start.elapsed();
    println!("[Duet] Global ensemble aggregation: {}", format_time_used(global_elapsed));

    annotate_window_results(&mut window_results, &cluster_to_global);

    let output_start = Instant::now();
    write_outputs(
        args,
        logger,
        &window_results,
        &global_ensembles,
        &global_base_records,
    )?;
    let output_elapsed = output_start.elapsed();
    println!("[Duet] Output writing: {}", format_time_used(output_elapsed));

    let total_elapsed = start_time.elapsed();
    println!("[Duet] Total time: {}", format_time_used(total_elapsed));
    logger.log(&format!("Total time: {:.2}s", total_elapsed.as_secs_f64()))?;
    
    Ok(())
}

fn load_reactivity_table(path: &str) -> Result<HashMap<String, ReactivityEntry>, Box<dyn Error>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut lines = reader.lines();
    let header = lines.next().ok_or("Error: Reactivity CSV is empty")??;
    if !header.to_lowercase().contains("reactivity") {
        eprintln!(
            "[Duet] Warning: Reactivity header does not contain expected columns: {}",
            header
        );
    }

    let mut map: HashMap<String, ReactivityEntry> = HashMap::new();
    for (idx, line) in lines.enumerate() {
        let raw = line?;
        if raw.trim().is_empty() {
            continue;
        }
        let fields: Vec<&str> = raw.split(',').collect();
        if fields.len() < 6 {
            continue;
        }
        let chr = fields[0].trim().to_string();
        let strand = fields[1].trim().to_string();
        let position = fields[2].trim().parse::<u32>().unwrap_or(0);
        let base = fields[3].trim().chars().next().unwrap_or('N');
        let stop_reactivity = fields[4].trim().parse::<f64>().unwrap_or(0.0);
        let mutation_reactivity = fields[5].trim().parse::<f64>().unwrap_or(0.0);
        if !stop_reactivity.is_finite() || !mutation_reactivity.is_finite() {
            continue;
        }
        let key = make_key(&chr, &strand, position);
        map.entry(key).or_insert(ReactivityEntry {
            chr,
            strand,
            position,
            _base: base,
            stop_reactivity,
            mutation_reactivity,
            stats: PositionStats::default(),
        });

        if idx > 0 && idx % 200_000 == 0 {
            println!("[Duet]   parsed {} rows", idx);
        }
    }
    Ok(map)
}

fn load_reference_sequences(path: &str) -> Result<HashMap<String, Vec<u8>>, Box<dyn Error>> {
    let reader = fasta::Reader::from_file(path)?;
    let mut sequences = HashMap::new();
    for record in reader.records() {
        let record = record?;
        let id = record.id().to_string();
        let seq = record
            .seq()
            .iter()
            .map(|b| b.to_ascii_uppercase())
            .collect();
        sequences.insert(id, seq);
    }
    Ok(sequences)
}

fn build_window_plan(
    positions: &HashMap<String, ReactivityEntry>,
    window_size: u32,
    window_step: u32,
) -> WindowPlan {
    let mut grouped: HashMap<String, Vec<ReactivityEntry>> = HashMap::new();
    for entry in positions.values() {
        let key = format!("{}|{}", entry.chr, entry.strand);
        grouped.entry(key).or_default().push(entry.clone());
    }

    let mut descriptors: Vec<WindowDescriptor> = Vec::new();
    let mut position_to_windows: HashMap<String, HashMap<u32, Vec<usize>>> = HashMap::new();

    for (transcript_key, mut entries) in grouped {
        if entries.is_empty() {
            continue;
        }
        entries.sort_by_key(|e| e.position);
        let start_pos = entries.first().map(|e| e.position).unwrap_or(0);
        let end_pos = entries.last().map(|e| e.position).unwrap_or(start_pos);
        if start_pos == 0 && end_pos == 0 {
            continue;
        }
        let parts: Vec<&str> = transcript_key.split('|').collect();
        if parts.len() != 2 {
            continue;
        }
        let chr = parts[0].to_string();
        let strand = parts[1].to_string();

        let mut window_defs: Vec<(usize, u32, u32)> = Vec::new();
        let mut current_start = start_pos;
        loop {
            let end = current_start.saturating_add(window_size.saturating_sub(1));
            let window_end = if end_pos == 0 { end } else { end.min(end_pos) };
            let descriptor = WindowDescriptor {
                chr: chr.clone(),
                strand: strand.clone(),
                start: current_start,
                end: window_end,
            };
            let index = descriptors.len();
            descriptors.push(descriptor);
            window_defs.push((index, current_start, window_end));

            if current_start >= end_pos {
                break;
            }
            let next_start = current_start.saturating_add(window_step);
            if next_start <= current_start {
                break;
            }
            if next_start > end_pos {
                current_start = end_pos;
            } else {
                current_start = next_start;
            }
        }

        let mut mapping: HashMap<u32, Vec<usize>> = HashMap::new();
        for entry in entries {
            let mut covering: Vec<usize> = Vec::new();
            for &(idx, w_start, w_end) in &window_defs {
                if entry.position >= w_start && entry.position <= w_end {
                    covering.push(idx);
                }
            }
            if !covering.is_empty() {
                mapping.insert(entry.position, covering);
            }
        }
        position_to_windows.insert(transcript_key, mapping);
    }

    WindowPlan {
        descriptors,
        position_to_windows,
    }
}

fn compute_window_reactivity_stats(
    plan: &WindowPlan,
    positions: &HashMap<String, ReactivityEntry>,
) -> Vec<WindowReactivityStats> {
    let mut stats = vec![WindowReactivityStats::default(); plan.len()];
    let mut counts = vec![0usize; plan.len()];

    // Parallel iteration over positions - collect updates first
    let entries: Vec<_> = positions.values().collect();
    let updates: Vec<(usize, f64, f64)> = entries
        .par_iter()
        .flat_map(|entry| {
            let transcript_key = format!("{}|{}", entry.chr, entry.strand);
            if let Some(indices) = plan.window_indices(&transcript_key, entry.position) {
                indices.iter().map(|&idx| (idx, entry.stop_reactivity, entry.mutation_reactivity)).collect::<Vec<_>>()
            } else {
                Vec::new()
            }
        })
        .collect();

    // Aggregate updates (sequential to avoid race conditions)
    for (idx, stop_val, mut_val) in updates {
        stats[idx].mean_stop += stop_val;
        stats[idx].mean_mutation += mut_val;
        counts[idx] += 1;
    }

    for (idx, count) in counts.iter().enumerate() {
        if *count > 0 {
            stats[idx].mean_stop /= *count as f64;
            stats[idx].mean_mutation /= *count as f64;
        }
    }

    stats
}

fn accumulate_covariation(
    bam_path: &str,
    reference: &HashMap<String, Vec<u8>>,
    positions: &mut HashMap<String, ReactivityEntry>,
    plan: &WindowPlan,
    window_features: &mut [Vec<ReadFeature>],
    read_features: &mut [Vec<WindowReadFeature>],
) -> Result<(), Box<dyn Error>> {
    let mut reader = bam::Reader::from_path(bam_path)?;
    let header = reader.header().clone();
    let target_count = header.target_count();
    let mut tid_to_name: HashMap<i32, String> = HashMap::with_capacity(target_count as usize);
    for tid in 0..target_count {
        let name = String::from_utf8_lossy(header.tid2name(tid)).to_string();
        tid_to_name.insert(tid as i32, name);
    }

    let mut read_counter: usize = 0;
    let mut dual_signal_read_counter: usize = 0;
    for result in reader.records() {
        let record = result?;
        if record.is_unmapped() {
            continue;
        }
        let tid = record.tid();
        let chrom = match tid_to_name.get(&tid) {
            Some(name) => name.clone(),
            None => continue,
        };
        let ref_seq = match reference.get(&chrom) {
            Some(seq) => seq,
            None => continue,
        };

        let strand = if record.is_reverse() { "-" } else { "+" };
        let transcript_key = format!("{}|{}", chrom, strand);
        let Some(pos_mapping) = plan.position_to_windows.get(&transcript_key) else {
            continue;
        };

        // Determine stop position (if any) and ensure it overlaps a tracked window.
        let mut stop_positions: Vec<u32> = Vec::new();
        let stop_pos = if strand == "+" {
            if record.pos() < 0 {
                None
            } else {
                Some((record.pos() + 1) as u32)
            }
        } else {
            let end_pos = record.cigar().end_pos();
            if end_pos <= 0 {
                None
            } else {
                Some(end_pos as u32)
            }
        };
        if let Some(pos) = stop_pos {
            if pos_mapping.contains_key(&pos) {
                stop_positions.push(pos);
            }
        }

        // Enumerate mutation positions (per-window)
        let mut mutation_positions_set: HashSet<u32> = HashSet::new();
        let read_bases = record.seq().as_bytes();
        let mut ref_pos = record.pos();
        let mut read_pos: usize = 0;
        for cigar in record.cigar().iter() {
            match cigar {
                Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) => {
                    for i in 0..(*len as i64) {
                        let ref_index = ref_pos + i;
                        if ref_index < 0 {
                            continue;
                        }
                        let ref_index_usize = ref_index as usize;
                        if ref_index_usize >= ref_seq.len() {
                            break;
                        }
                        let read_index = read_pos + i as usize;
                        if read_index >= read_bases.len() {
                            break;
                        }
                        let ref_base = ref_seq[ref_index_usize];
                        let read_base = read_bases[read_index].to_ascii_uppercase();
                        if read_base != ref_base && read_base != b'N' && ref_base != b'N' {
                            let position = (ref_index + 1) as u32;
                            if pos_mapping.contains_key(&position) {
                                mutation_positions_set.insert(position);
                            }
                        }
                    }
                    ref_pos += *len as i64;
                    read_pos += *len as usize;
                }
                Cigar::Ins(len) | Cigar::SoftClip(len) => {
                    read_pos += *len as usize;
                }
                Cigar::Del(len) | Cigar::RefSkip(len) => {
                    ref_pos += *len as i64;
                }
                Cigar::HardClip(_) | Cigar::Pad(_) => {}
            }
        }
        let mutation_positions: Vec<u32> = mutation_positions_set.into_iter().collect();

        if stop_positions.is_empty() && mutation_positions.is_empty() {
            continue;
        }

        // Update per-position statistics
        for pos in &stop_positions {
            let key = make_key(&chrom, strand, *pos);
            if let Some(entry) = positions.get_mut(&key) {
                entry.stats.stop_reads += 1;
                entry.stats.co_stop_mut += mutation_positions.len();
            }
        }
        for pos in &mutation_positions {
            let key = make_key(&chrom, strand, *pos);
            if let Some(entry) = positions.get_mut(&key) {
                entry.stats.mutation_reads += 1;
                entry.stats.co_mut_stop += stop_positions.len();
            }
        }

        // Aggregate read-level features per window
        let mut accumulators: HashMap<usize, WindowReadAccumulator> = HashMap::new();
        for pos in &stop_positions {
            if let Some(indices) = pos_mapping.get(pos) {
                for &idx in indices {
                    accumulators
                        .entry(idx)
                        .or_insert_with(WindowReadAccumulator::default)
                        .stop_positions
                        .insert(*pos);
                }
            }
        }
        for pos in &mutation_positions {
            if let Some(indices) = pos_mapping.get(pos) {
                for &idx in indices {
                    accumulators
                        .entry(idx)
                        .or_insert_with(WindowReadAccumulator::default)
                        .mutation_positions
                        .insert(*pos);
                }
            }
        }

        // Check if this read has dual signals (both stop and mutation)
        let has_dual_signal = !stop_positions.is_empty() && !mutation_positions.is_empty();
        if has_dual_signal {
            dual_signal_read_counter += 1;
        }

        for (idx, acc) in accumulators {
            let stop_positions: Vec<u32> = acc.stop_positions.into_iter().collect();
            let mutation_positions: Vec<u32> = acc.mutation_positions.into_iter().collect();
            if stop_positions.is_empty() && mutation_positions.is_empty() {
                continue;
            }
            if let Some(bucket) = window_features.get_mut(idx) {
                bucket.push(ReadFeature {
                    stop_count: stop_positions.len(),
                    mutation_count: mutation_positions.len(),
                });
            }
            if let Some(bucket) = read_features.get_mut(idx) {
                bucket.push(WindowReadFeature {
                    stop_positions,
                    mutation_positions,
                });
            }
        }

        read_counter += 1;
        if read_counter % 1_000_000 == 0 {
            println!("[Duet]   processed {} reads", read_counter);
        }
    }

    println!("[Duet] Processed {} aligned reads", read_counter);
    let dual_signal_fraction = if read_counter > 0 {
        dual_signal_read_counter as f64 / read_counter as f64
    } else {
        0.0
    };
    println!("[Duet]   {} reads ({:.2}%) have dual signals (both stop and mutation)", 
             dual_signal_read_counter, dual_signal_fraction * 100.0);
    Ok(())
}

fn analyze_windows(
    plan: &WindowPlan,
    window_features: &[Vec<ReadFeature>],
    window_read_positions: &[Vec<WindowReadFeature>],
    window_reactivity: &[WindowReactivityStats],
    epsilon: f64,
    min_samples: usize,
    positions: &HashMap<String, ReactivityEntry>,
) -> (Vec<WindowResult>, Vec<ClusterDetail>, (usize, usize, f64)) {
    // Returns: (window_results, cluster_details, (high_confidence_positions, total_positions, avg_confidence))
    // Parallel processing of windows
    let results: Vec<(WindowResult, Vec<ClusterDetail>)> = (0..plan.len())
        .into_par_iter()
        .map(|idx| {
        let descriptor = plan.descriptor(idx).clone();
        let features = &window_features[idx];
        let position_features = &window_read_positions[idx];
        let reactivity = &window_reactivity[idx];
        let read_count = features.len();
        let dual_signal_count = features.iter().filter(|f| f.co_occurrence()).count();
        let dual_signal_fraction = if read_count > 0 {
            dual_signal_count as f64 / read_count as f64
        } else {
            0.0
        };

        if read_count == 0 {
            return (WindowResult {
                descriptor,
                status: WindowStatus::InsufficientReads,
                read_count,
                ensemble_count: 0,
                primary_id: None,
                primary_occupancy: 0.0,
                noise_fraction: 0.0,
                dual_signal_fraction: 0.0,
                mean_stop_reactivity: reactivity.mean_stop,
                mean_mutation_reactivity: reactivity.mean_mutation,
                clusters: Vec::new(),
            }, Vec::new());
        }

        if read_count < min_samples {
            let mut clusters: Vec<ClusterSummary> = Vec::new();
            let (stop_sum, mut_sum) = features.iter().fold((0usize, 0usize), |acc, feat| {
                (acc.0 + feat.stop_count, acc.1 + feat.mutation_count)
            });
            let avg_stop = stop_sum as f64 / read_count as f64;
            let avg_mut = mut_sum as f64 / read_count as f64;
            clusters.push(ClusterSummary {
                id: 1,
                count: read_count,
                occupancy: 1.0,
                mean_stop_count: avg_stop,
                mean_mutation_count: avg_mut,
                dual_signal_fraction,
                is_noise: false,
                confidence: 0.0,
                global_ids: Vec::new(),
                cluster_key: None,
            });

            return (WindowResult {
                descriptor,
                status: WindowStatus::InsufficientReads,
                read_count,
                ensemble_count: 1,
                primary_id: Some(1),
                primary_occupancy: 1.0,
                noise_fraction: 0.0,
                dual_signal_fraction,
                mean_stop_reactivity: reactivity.mean_stop,
                mean_mutation_reactivity: reactivity.mean_mutation,
                clusters,
            }, Vec::new());
        }

        let feature_matrix: Vec<Vec<f64>> =
            features.iter().copied().map(ReadFeature::to_vec).collect();
        let standardized = standardize_features(&feature_matrix);
        let raw_labels = dbscan(&standardized, epsilon, min_samples);
        let mut labels: Vec<usize> = if raw_labels.is_empty() {
            vec![0; read_count]
        } else {
            raw_labels
                .into_iter()
                .map(|label| if label >= 1 { label as usize } else { 0usize })
                .collect()
        };

        if labels.len() != read_count {
            labels = vec![0; read_count];
        }

        let unique_clusters: HashSet<usize> =
            labels.iter().copied().filter(|&id| id >= 1).collect();
        let status = if unique_clusters.is_empty() {
            WindowStatus::LowDensity
        } else {
            WindowStatus::Ok
        };

        let mut cluster_map: HashMap<usize, ClusterSummary> = HashMap::new();
        let mut detail_map: HashMap<usize, ClusterDetail> = HashMap::new();
        let mut noise_count = 0usize;
        for i in 0..read_count {
            let label = labels[i];
            let feature = features[i];
            let pos_feature = position_features
                .get(i)
                .cloned()
                .unwrap_or(WindowReadFeature {
                    stop_positions: Vec::new(),
                    mutation_positions: Vec::new(),
                });

            let entry = cluster_map.entry(label).or_insert_with(|| ClusterSummary {
                id: label,
                count: 0,
                occupancy: 0.0,
                mean_stop_count: 0.0,
                mean_mutation_count: 0.0,
                dual_signal_fraction: 0.0,
                is_noise: label == 0,
                confidence: 0.0,
                global_ids: Vec::new(),
                cluster_key: None,
            });
            entry.count += 1;
            entry.mean_stop_count += feature.stop_count as f64;
            entry.mean_mutation_count += feature.mutation_count as f64;
            if feature.co_occurrence() {
                entry.dual_signal_fraction += 1.0;
            }

            if label == 0 {
                noise_count += 1;
                continue;
            }

            let detail_entry = detail_map.entry(label).or_insert_with(|| ClusterDetail {
                key: ClusterKey {
                    window_index: idx,
                    local_id: label,
                },
                chr: descriptor.chr.clone(),
                strand: descriptor.strand.clone(),
                _window_start: descriptor.start,
                _window_end: descriptor.end,
                read_count: 0,
                dual_signal_reads: 0,
                positions: HashSet::new(),
                position_stats: HashMap::new(),
                position_confidence: HashMap::new(),
                cluster_confidence: 0.0,
            });
            detail_entry.read_count += 1;
            if feature.co_occurrence() {
                detail_entry.dual_signal_reads += 1;
            }
            for pos in pos_feature.stop_positions {
                detail_entry.positions.insert(pos);
                detail_entry
                    .position_stats
                    .entry(pos)
                    .or_insert_with(PositionContribution::default)
                    .stop_reads += 1;
            }
            for pos in pos_feature.mutation_positions {
                detail_entry.positions.insert(pos);
                detail_entry
                    .position_stats
                    .entry(pos)
                    .or_insert_with(PositionContribution::default)
                    .mutation_reads += 1;
            }
        }

        let mut clusters: Vec<ClusterSummary> = cluster_map
            .into_iter()
            .map(|(id, mut summary)| {
                if summary.count > 0 {
                    summary.occupancy = summary.count as f64 / read_count as f64;
                    summary.mean_stop_count /= summary.count as f64;
                    summary.mean_mutation_count /= summary.count as f64;
                    summary.dual_signal_fraction /= summary.count as f64;
                }
                summary.id = id;
                summary
            })
            .collect();

        clusters.sort_by(|a, b| {
            b.occupancy
                .partial_cmp(&a.occupancy)
                .unwrap_or(Ordering::Equal)
        });

        for cluster in clusters.iter_mut() {
            if !cluster.is_noise && cluster.count > 0 {
                if let Some(detail) = detail_map.get_mut(&cluster.id) {
                    let confidence = score_cluster_detail(detail, positions);
                    detail.cluster_confidence = confidence;
                    cluster.confidence = confidence;
                    cluster.cluster_key = Some(detail.key);
                }
            }
        }

        let ensemble_count = clusters
            .iter()
            .filter(|c| !c.is_noise && c.count > 0)
            .count();
        let primary_cluster = clusters
            .iter()
            .find(|c| !c.is_noise && c.count > 0)
            .map(|c| c.id);
        let primary_occupancy = clusters
            .iter()
            .find(|c| !c.is_noise && c.count > 0)
            .map(|c| c.occupancy)
            .unwrap_or(0.0);
        let noise_fraction = if read_count > 0 {
            noise_count as f64 / read_count as f64
        } else {
            0.0
        };

        let window_cluster_details: Vec<ClusterDetail> = detail_map.into_values().collect();

        (WindowResult {
            descriptor,
            status,
            read_count,
            ensemble_count,
            primary_id: primary_cluster,
            primary_occupancy,
            noise_fraction,
            dual_signal_fraction,
            mean_stop_reactivity: reactivity.mean_stop,
            mean_mutation_reactivity: reactivity.mean_mutation,
            clusters,
        }, window_cluster_details)
        })
        .collect();

    // Separate results and cluster_details, and collect confidence statistics
    let mut window_results = Vec::with_capacity(results.len());
    let mut cluster_details = Vec::new();
    let mut high_confidence_positions = 0usize;
    let mut total_positions = 0usize;
    let mut confidence_sum = 0.0f64;
    
    for (result, details) in results {
        window_results.push(result);
        for detail in &details {
            total_positions += detail.position_confidence.len();
            for confidence in detail.position_confidence.values() {
                confidence_sum += confidence;
                if *confidence >= HIGH_REACTIVITY_THRESHOLD {
                    high_confidence_positions += 1;
                }
            }
        }
        cluster_details.extend(details);
    }

    let avg_confidence = if total_positions > 0 {
        confidence_sum / total_positions as f64
    } else {
        0.0
    };

    (window_results, cluster_details, (high_confidence_positions, total_positions, avg_confidence))
}

fn write_outputs(
    args: &DuetArgs,
    logger: &mut Logger,
    window_results: &[WindowResult],
    global_ensembles: &[GlobalEnsemble],
    global_bases: &[GlobalBaseRecord],
) -> Result<(), Box<dyn Error>> {
    let mut output = File::create(&args.output)?;
    writeln!(
        output,
        "ChrID,Strand,WindowStart,WindowEnd,Status,ReadCount,EnsembleCount,PrimaryEnsembleID,PrimaryOccupancy,NoiseFraction,DualSignalFraction,MeanStopReactivity,MeanMutationReactivity,ClusterSummary,GlobalEnsembles"
    )?;

    for result in window_results {
        let descriptor = &result.descriptor;
        let primary_id = result
            .primary_id
            .map(|id| id.to_string())
            .unwrap_or_else(|| "NA".to_string());
        let mut summary_parts: Vec<String> = result
            .clusters
            .iter()
            .filter(|c| !c.is_noise && c.count > 0)
            .map(|c| {
                format!(
                    "c{}(occ={:.3},conf={:.3},stop={:.2},mut={:.2})",
                    c.id, c.occupancy, c.confidence, c.mean_stop_count, c.mean_mutation_count
                )
            })
            .collect();
        if let Some(noise) = result.clusters.iter().find(|c| c.is_noise && c.count > 0) {
            summary_parts.push(format!("noise={:.3}", noise.occupancy));
        }
        if summary_parts.is_empty() {
            summary_parts.push("NA".to_string());
        }

        let mut global_parts: Vec<String> = result
            .clusters
            .iter()
            .filter(|c| !c.global_ids.is_empty())
            .map(|c| {
                let ids: Vec<String> = c.global_ids.iter().map(|id| format!("E{}", id)).collect();
                format!("c{}->{}", c.id, ids.join("|"))
            })
            .collect();
        if global_parts.is_empty() {
            global_parts.push("NA".to_string());
        }

        writeln!(
            output,
            "{},{},{},{},{},{},{},{},{:.6},{:.6},{:.6},{:.6},{:.6},{},{}",
            descriptor.chr,
            descriptor.strand,
            descriptor.start,
            descriptor.end,
            result.status,
            result.read_count,
            result.ensemble_count,
            primary_id,
            result.primary_occupancy,
            result.noise_fraction,
            result.dual_signal_fraction,
            result.mean_stop_reactivity,
            result.mean_mutation_reactivity,
            summary_parts.join(";"),
            global_parts.join(";")
        )?;
    }

    let summary_path = args
        .summary_output
        .clone()
        .unwrap_or_else(|| default_summary_path(&args.output));
    let mut summary = File::create(&summary_path)?;
    writeln!(
        summary,
        "ChrID,Strand,WindowStart,WindowEnd,Status,ReadCount,EnsembleID,GlobalEnsembleID,IsNoise,Members,Occupancy,ClusterConfidence,MeanStopCount,MeanMutationCount,DualSignalFraction,MeanStopReactivity,MeanMutationReactivity"
    )?;

    for result in window_results {
        let descriptor = &result.descriptor;
        if result.clusters.is_empty() {
            writeln!(
                summary,
                "{},{},{},{},{},{},NA,NA,false,0,0.000,0.000,0.000,0.000,0.000,{:.6},{:.6}",
                descriptor.chr,
                descriptor.strand,
                descriptor.start,
                descriptor.end,
                result.status,
                result.read_count,
                result.mean_stop_reactivity,
                result.mean_mutation_reactivity
            )?;
            continue;
        }

        for cluster in &result.clusters {
            writeln!(
                summary,
                "{},{},{},{},{},{},{},{},{},{},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6}",
                descriptor.chr,
                descriptor.strand,
                descriptor.start,
                descriptor.end,
                result.status,
                result.read_count,
                cluster.id,
                cluster
                    .global_ids
                    .first()
                    .map(|id| format!("E{}", id))
                    .unwrap_or_else(|| "NA".to_string()),
                cluster.is_noise,
                cluster.count,
                cluster.occupancy,
                cluster.confidence,
                cluster.mean_stop_count,
                cluster.mean_mutation_count,
                cluster.dual_signal_fraction,
                result.mean_stop_reactivity,
                result.mean_mutation_reactivity
            )?;
        }
    }

    logger.log(&format!(
        "Duet window assessment written to {}",
        args.output
    ))?;
    logger.log(&format!("Duet window summary written to {}", summary_path))?;
    println!("[Duet] Output written to {}", args.output);
    println!("[Duet] Summary written to {}", summary_path);

    let global_path = if let Some(stripped) = args.output.strip_suffix(".csv") {
        format!("{}_global.csv", stripped)
        } else {
        format!("{}_global.csv", args.output)
    };
    let mut global_file = File::create(&global_path)?;
    writeln!(
        global_file,
        "ChrID,Strand,GlobalEnsembleID,WindowCount,ClusterCount,ReadCount,TotalStopReads,TotalMutationReads,DualSignalReads,UniquePositions,Confidence"
    )?;
    for ensemble in global_ensembles {
        writeln!(
            global_file,
            "{},{},E{},{},{},{},{},{},{},{},{}",
            ensemble.chr,
            ensemble.strand,
            ensemble.id,
            ensemble.window_count,
            ensemble.cluster_count,
            ensemble.read_count,
            ensemble.total_stop_reads,
            ensemble.total_mutation_reads,
            ensemble.dual_signal_reads,
            ensemble.unique_positions,
            ensemble.confidence
        )?;
    }

    let base_path = if let Some(stripped) = args.output.strip_suffix(".csv") {
        format!("{}_global_per_base.csv", stripped)
        } else {
        format!("{}_global_per_base.csv", args.output)
    };
    let mut base_file = File::create(&base_path)?;
        writeln!(
        base_file,
        "ChrID,Strand,Position,GlobalEnsembleID,StopReads,MutationReads,StopReactivity,MutationReactivity,PositionConfidence"
    )?;
    for record in global_bases {
        writeln!(
            base_file,
            "{},{},{},E{},{},{},{:.6},{:.6},{:.6}",
            record.chr,
            record.strand,
            record.position,
            record.ensemble_id,
            record.stop_reads,
            record.mutation_reads,
            record.stop_reactivity,
            record.mutation_reactivity,
            record.position_confidence
        )?;
    }

    logger.log(&format!("Duet global ensembles written to {}", global_path))?;
    logger.log(&format!(
        "Duet global per-base details written to {}",
        base_path
    ))?;
    println!("[Duet] Global ensembles written to {}", global_path);
    println!("[Duet] Global per-base details written to {}", base_path);

    // Generate heatmap for global ensembles reactivity distribution
    let heatmap_path = if let Some(stripped) = args.output.strip_suffix(".csv") {
        format!("{}_global_heatmap.svg", stripped)
    } else {
        format!("{}_global_heatmap.svg", args.output)
    };
    
    if !global_ensembles.is_empty() && !global_bases.is_empty() {
        match plot_global_ensemble_heatmap(global_ensembles, global_bases, &heatmap_path) {
            Ok(_) => {
                logger.log(&format!("Duet global ensemble heatmap written to {}", heatmap_path))?;
                println!("[Duet] Global ensemble heatmap written to {}", heatmap_path);
            }
            Err(e) => {
                eprintln!("[Duet Warning] Failed to generate heatmap: {}", e);
                logger.log(&format!("Warning: Failed to generate heatmap: {}", e))?;
            }
        }
    }

    // Generate t-SNE/PCA visualization for global ensembles reactivity distribution
    let tsne_path = if let Some(stripped) = args.output.strip_suffix(".csv") {
        format!("{}_global_tsne.svg", stripped)
    } else {
        format!("{}_global_tsne.svg", args.output)
    };
    
    if !global_ensembles.is_empty() && !global_bases.is_empty() {
        match plot_global_ensemble_tsne(global_ensembles, global_bases, &tsne_path) {
            Ok(_) => {
                logger.log(&format!("Duet global ensemble t-SNE plot written to {}", tsne_path))?;
                println!("[Duet] Global ensemble t-SNE plot written to {}", tsne_path);
            }
            Err(e) => {
                eprintln!("[Duet Warning] Failed to generate t-SNE plot: {}", e);
                logger.log(&format!("Warning: Failed to generate t-SNE plot: {}", e))?;
            }
        }
    }

    // Generate overlap track plots for non-full-length ensembles
    let track_dir = if let Some(stripped) = args.output.strip_suffix(".csv") {
        format!("{}_tracks", stripped)
    } else {
        format!("{}_tracks", args.output)
    };
    
    if !global_ensembles.is_empty() && !global_bases.is_empty() {
        match plot_ensemble_overlap_tracks(global_ensembles, global_bases, &track_dir) {
            Ok(_) => {
                logger.log(&format!("Duet ensemble overlap tracks written to {}", track_dir))?;
                println!("[Duet] Ensemble overlap tracks written to {}", track_dir);
            }
            Err(e) => {
                eprintln!("[Duet Warning] Failed to generate overlap tracks: {}", e);
                logger.log(&format!("Warning: Failed to generate overlap tracks: {}", e))?;
            }
        }
    }

    Ok(())
}

fn standardize_features(data: &[Vec<f64>]) -> Vec<Vec<f64>> {
    if data.is_empty() || data[0].is_empty() {
        return data.to_vec();
    }
    let dims = data[0].len();
    let mut means = vec![0.0; dims];
    let mut stds = vec![0.0; dims];
    let n = data.len() as f64;
    for vec in data {
        for (idx, value) in vec.iter().enumerate() {
            means[idx] += value;
        }
    }
    for mean in &mut means {
        *mean /= n;
    }
    for vec in data {
        for (idx, value) in vec.iter().enumerate() {
            let diff = value - means[idx];
            stds[idx] += diff * diff;
        }
    }
    for std in &mut stds {
        *std = (*std / n).sqrt();
        if *std < 1e-6 {
            *std = 1.0;
        }
    }
    let mut standardized = Vec::with_capacity(data.len());
    for vec in data {
        let mut transformed = Vec::with_capacity(dims);
        for idx in 0..dims {
            transformed.push((vec[idx] - means[idx]) / stds[idx]);
        }
        standardized.push(transformed);
    }
    standardized
}

fn dbscan(data: &[Vec<f64>], eps: f64, min_points: usize) -> Vec<isize> {
    if data.is_empty() {
        return Vec::new();
    }
    let mut labels = vec![-1isize; data.len()];
    let mut cluster_id: isize = 0;
    let mut visited = vec![false; data.len()];
    for point_idx in 0..data.len() {
        if visited[point_idx] {
            continue;
        }
        visited[point_idx] = true;
        let mut neighbours = region_query(data, point_idx, eps);
        if neighbours.len() < min_points {
            labels[point_idx] = -1;
            continue;
        }
        cluster_id += 1;
        labels[point_idx] = cluster_id;
        let mut i = 0;
        while i < neighbours.len() {
            let neighbour_idx = neighbours[i];
            if !visited[neighbour_idx] {
                visited[neighbour_idx] = true;
                let neighbour_neighbours = region_query(data, neighbour_idx, eps);
                if neighbour_neighbours.len() >= min_points {
                    neighbours.extend(neighbour_neighbours);
                }
            }
            if labels[neighbour_idx] == -1 {
                labels[neighbour_idx] = cluster_id;
            }
            i += 1;
        }
    }
    labels
}

fn region_query(data: &[Vec<f64>], point_idx: usize, eps: f64) -> Vec<usize> {
    let mut neighbours = Vec::new();
    for (idx, other) in data.iter().enumerate() {
        if euclidean_distance(&data[point_idx], other) <= eps {
            neighbours.push(idx);
        }
    }
    neighbours
}

fn euclidean_distance(a: &[f64], b: &[f64]) -> f64 {
    a.iter()
        .zip(b.iter())
        .map(|(x, y)| {
            let diff = x - y;
            diff * diff
        })
        .sum::<f64>()
        .sqrt()
}

fn default_summary_path(output_path: &str) -> String {
    if let Some(stripped) = output_path.strip_suffix(".csv") {
        format!("{}_summary.csv", stripped)
            } else {
        format!("{}_summary.csv", output_path)
    }
}

fn make_key(chr: &str, strand: &str, position: u32) -> String {
    format!("{chr}|{strand}|{position}")
}

fn annotate_window_results(
    window_results: &mut [WindowResult],
    cluster_to_global: &HashMap<ClusterKey, usize>,
) {
    for result in window_results.iter_mut() {
        for cluster in result.clusters.iter_mut() {
            if let Some(key) = cluster.cluster_key {
                if let Some(id) = cluster_to_global.get(&key) {
                    cluster.global_ids = vec![*id];
                }
            }
        }
    }
}

fn aggregate_global_ensembles(
    cluster_details: &[ClusterDetail],
    positions: &HashMap<String, ReactivityEntry>,
    min_reactivity_correlation: f64,
    min_shared_fraction: f64,
    min_max_fraction: f64,
) -> (
    Vec<GlobalEnsemble>,
    HashMap<ClusterKey, usize>,
    Vec<GlobalBaseRecord>,
) {
    if cluster_details.is_empty() {
        return (Vec::new(), HashMap::new(), Vec::new());
    }

    // Build position map in parallel first
    let position_map_entries: Vec<((String, String, u32), usize)> = cluster_details
        .par_iter()
        .enumerate()
        .flat_map(|(idx, detail)| {
            let positions_vec: Vec<u32> = detail.positions.iter().copied().collect();
            positions_vec.into_par_iter().map(move |pos| {
                ((detail.chr.clone(), detail.strand.clone(), pos), idx)
            })
        })
        .collect();
    
    let mut position_map: HashMap<(String, String, u32), Vec<usize>> = HashMap::new();
    for (key, idx) in position_map_entries {
        position_map.entry(key).or_default().push(idx);
    }

    // Build cluster overlaps from position map (sequential to avoid race conditions)
    let mut cluster_overlaps: HashMap<(usize, usize), usize> = HashMap::new();
    for indices in position_map.values() {
        if indices.len() > 1 {
            for i in 0..indices.len() {
                for j in (i + 1)..indices.len() {
                    let pair = (indices[i], indices[j]);
                    *cluster_overlaps.entry(pair).or_insert(0) += 1;
                }
            }
        }
    }

    // Only union clusters that share enough positions AND don't have conflicting exclusive positions
    // Use a more conservative approach: require higher overlap fraction to merge clusters
    let mut uf = UnionFind::new(cluster_details.len());
    let mut merge_count = 0usize;
    let mut conflict_count = 0usize;
    
    for ((idx1, idx2), shared_count) in cluster_overlaps.iter() {
        let detail1 = &cluster_details[*idx1];
        let detail2 = &cluster_details[*idx2];
        let min_positions = detail1.positions.len().min(detail2.positions.len());
        let max_positions = detail1.positions.len().max(detail2.positions.len());
        
        // More conservative: require higher overlap fraction (configurable via command line)
        // AND shared positions should be substantial relative to larger cluster to ensure substantial overlap
        let meets_count_threshold = *shared_count >= MIN_SHARED_POSITIONS;
        let meets_fraction_threshold = if min_positions > 0 {
            (*shared_count as f64 / min_positions as f64) >= min_shared_fraction
        } else {
            false
        };
        // Additional check: shared positions should be substantial relative to larger cluster
        let meets_max_fraction = if max_positions > 0 {
            (*shared_count as f64 / max_positions as f64) >= min_max_fraction
        } else {
            false
        };
        
        // Check for conflicting exclusive positions
        // If two clusters have different positions that overlap spatially (within a small window),
        // they represent different structural choices and should not be merged
        let has_conflicting_positions = {
            let pos1_exclusive: Vec<u32> = detail1.positions.difference(&detail2.positions).copied().collect();
            let pos2_exclusive: Vec<u32> = detail2.positions.difference(&detail1.positions).copied().collect();
            
            // Check if exclusive positions from different clusters are spatially close
            // (within 5nt), indicating they might be alternative choices at the same location
            let mut has_conflict = false;
            for &p1 in &pos1_exclusive {
                for &p2 in &pos2_exclusive {
                    // If positions are close (within 5nt) and both clusters have significant support,
                    // they likely represent alternative structural choices
                    if p1.abs_diff(p2) <= 5 {
                        // Check if both positions have substantial read support
                        let p1_support = detail1.position_stats.get(&p1)
                            .map(|c| c.stop_reads + c.mutation_reads)
                            .unwrap_or(0);
                        let p2_support = detail2.position_stats.get(&p2)
                            .map(|c| c.stop_reads + c.mutation_reads)
                            .unwrap_or(0);
                        
                        // If both have at least 3 reads, consider them conflicting
                        if p1_support >= 3 && p2_support >= 3 {
                            has_conflict = true;
                            conflict_count += 1;
                            break;
                        }
                    }
                }
                if has_conflict {
                    break;
                }
            }
            has_conflict
        };
        
        // Check reactivity pattern similarity on shared positions (DREEM principle:
        // different structures have different mutation profiles)
        // Only merge if clusters show similar reactivity patterns on shared positions
        let reactivity_similar = {
            let shared_positions: Vec<u32> = detail1.positions
                .intersection(&detail2.positions)
                .copied()
                .collect();
            
            if shared_positions.len() < 3 {
                // Not enough shared positions to assess pattern similarity
                true  // Allow merge if we can't assess similarity
            } else {
                // Collect relative reactivity values for each cluster at shared positions
                // We compare the pattern: which positions are more/less reactive in each cluster
                let mut react1 = Vec::new();
                let mut react2 = Vec::new();
                
                // Normalize by cluster read counts to get relative reactivity
                let total_reads1 = detail1.read_count.max(1) as f64;
                let total_reads2 = detail2.read_count.max(1) as f64;
                
                for &pos in &shared_positions {
                    // Get read support for this position in each cluster
                    let contrib1 = detail1.position_stats.get(&pos)
                        .map(|c| (c.stop_reads + c.mutation_reads) as f64)
                        .unwrap_or(0.0);
                    let contrib2 = detail2.position_stats.get(&pos)
                        .map(|c| (c.stop_reads + c.mutation_reads) as f64)
                        .unwrap_or(0.0);
                    
                    // Normalize by cluster size to get relative reactivity pattern
                    // This reflects how important this position is within each cluster
                    let rel_react1 = contrib1 / total_reads1;
                    let rel_react2 = contrib2 / total_reads2;
                    
                    react1.push(rel_react1);
                    react2.push(rel_react2);
                }
                
                // Calculate Pearson correlation coefficient to measure pattern similarity
                if react1.len() >= 3 {
                    let mean1 = react1.iter().sum::<f64>() / react1.len() as f64;
                    let mean2 = react2.iter().sum::<f64>() / react2.len() as f64;
                    
                    let mut numerator = 0.0;
                    let mut var1 = 0.0;
                    let mut var2 = 0.0;
                    
                    for i in 0..react1.len() {
                        let diff1 = react1[i] - mean1;
                        let diff2 = react2[i] - mean2;
                        numerator += diff1 * diff2;
                        var1 += diff1 * diff1;
                        var2 += diff2 * diff2;
                    }
                    
                    if var1 > 1e-10 && var2 > 1e-10 {
                        let correlation = numerator / (var1.sqrt() * var2.sqrt());
                        correlation >= min_reactivity_correlation
                    } else {
                        // If variance is too low, patterns are too similar to distinguish
                        true
                    }
                } else {
                    true
                }
            }
        };
        
        // Only union if they share enough positions AND don't have conflicting positions
        // AND meet the maximum fraction threshold AND have similar reactivity patterns
        if (meets_count_threshold || meets_fraction_threshold) && meets_max_fraction && !has_conflicting_positions && reactivity_similar {
            uf.union(*idx1, *idx2);
            merge_count += 1;
        }
    }
    
    eprintln!("[Duet Debug] Cluster pairs considered: {}, merged: {}, conflicts detected: {}", 
              cluster_overlaps.len(), merge_count, conflict_count);

    let mut components: HashMap<usize, Vec<usize>> = HashMap::new();
    for idx in 0..cluster_details.len() {
        let root = uf.find(idx);
        components.entry(root).or_default().push(idx);
    }

    let mut global_ensembles: Vec<GlobalEnsemble> = Vec::new();
    let mut cluster_to_global: HashMap<ClusterKey, usize> = HashMap::new();
    let mut global_bases: Vec<GlobalBaseRecord> = Vec::new();

    let mut next_id = 1usize;
    for indices in components.values() {
        let first = &cluster_details[indices[0]];
        let mut window_set: HashSet<usize> = HashSet::new();
        let mut read_total = 0usize;
        let mut dual_total = 0usize;
        let mut base_map: HashMap<u32, PositionContribution> = HashMap::new();
        let mut position_confidence_sum: HashMap<u32, (f64, f64)> = HashMap::new();
        let mut cluster_conf_sum = 0.0;
        let mut cluster_conf_weight = 0.0;

        for &idx in indices {
            let detail = &cluster_details[idx];
            window_set.insert(detail.key.window_index);
            read_total += detail.read_count;
            dual_total += detail.dual_signal_reads;
            for (pos, contrib) in &detail.position_stats {
                let entry = base_map
                    .entry(*pos)
                    .or_insert_with(PositionContribution::default);
                entry.stop_reads += contrib.stop_reads;
                entry.mutation_reads += contrib.mutation_reads;
            }
            for (pos, conf) in &detail.position_confidence {
                let weight = detail
                    .position_stats
                    .get(pos)
                    .map(|c| (c.stop_reads + c.mutation_reads) as f64)
                    .unwrap_or(0.0);
                let entry = position_confidence_sum.entry(*pos).or_insert((0.0, 0.0));
                entry.0 += conf * weight;
                entry.1 += weight;
            }
            cluster_to_global.insert(detail.key, next_id);
            cluster_conf_sum += detail.cluster_confidence * detail.read_count as f64;
            cluster_conf_weight += detail.read_count as f64;
        }

        let total_stop_reads = base_map.values().map(|c| c.stop_reads).sum();
        let total_mutation_reads = base_map.values().map(|c| c.mutation_reads).sum();
        let unique_positions = base_map.len();
        let ensemble_confidence = if cluster_conf_weight > 0.0 {
            cluster_conf_sum / cluster_conf_weight
    } else {
            0.0
        };

        global_ensembles.push(GlobalEnsemble {
            id: next_id,
            chr: first.chr.clone(),
            strand: first.strand.clone(),
            window_count: window_set.len(),
            cluster_count: indices.len(),
            read_count: read_total,
            total_stop_reads,
            total_mutation_reads,
            dual_signal_reads: dual_total,
            unique_positions,
            confidence: ensemble_confidence,
        });

        for (pos, contrib) in &base_map {
            let key = format!("{}|{}|{}", first.chr, first.strand, pos);
            let entry = positions.get(&key);
            let (stop_reactivity, mutation_reactivity) = entry
                .map(|r| (r.stop_reactivity, r.mutation_reactivity))
                .unwrap_or((0.0, 0.0));
            let confidence = position_confidence_sum
                .get(pos)
                .and_then(|(sum, weight)| {
                    if *weight > 0.0 {
                        Some(sum / weight)
                    } else {
                        None
                    }
                })
                .unwrap_or(0.0);
            global_bases.push(GlobalBaseRecord {
                ensemble_id: next_id,
                chr: first.chr.clone(),
                strand: first.strand.clone(),
                position: *pos,
                stop_reads: contrib.stop_reads,
                mutation_reads: contrib.mutation_reads,
                stop_reactivity,
                mutation_reactivity,
                position_confidence: confidence,
            });
        }

        next_id += 1;
    }

    global_bases.sort_by(|a, b| {
        a.ensemble_id
            .cmp(&b.ensemble_id)
            .then(a.position.cmp(&b.position))
    });

    (global_ensembles, cluster_to_global, global_bases)
}

// Second aggregation: merge ensembles based on signal positions only
// Only considers positions where both ensembles have signals (reactivity > 0)
fn aggregate_ensembles_second_pass(
    first_ensembles: &[GlobalEnsemble],
    first_base_records: &[GlobalBaseRecord],
    positions: &HashMap<String, ReactivityEntry>,
    min_reactivity_correlation: f64,
    min_shared_fraction: f64,
    min_max_fraction: f64,
) -> (
    Vec<GlobalEnsemble>,
    HashMap<usize, usize>, // first_id -> global_id
    Vec<GlobalBaseRecord>,
) {
    if first_ensembles.is_empty() {
        return (Vec::new(), HashMap::new(), Vec::new());
    }

    // Build reactivity matrix for ensembles (only signal positions)
    let mut all_positions: Vec<u32> = first_base_records
        .iter()
        .map(|b| b.position)
        .collect::<HashSet<_>>()
        .into_iter()
        .collect();
    all_positions.sort();

    // Build reactivity matrix: only include positions with signals
    let reactivity_matrix: Vec<Vec<f64>> = first_ensembles.iter().map(|ensemble| {
        let mut react_vec = vec![0.0; all_positions.len()];
        let pos_to_idx: HashMap<u32, usize> = all_positions.iter().enumerate()
            .map(|(idx, &pos)| (pos, idx))
            .collect();
        
        for base in first_base_records.iter().filter(|b| b.ensemble_id == ensemble.id) {
            if let Some(&idx) = pos_to_idx.get(&base.position) {
                let reactivity = base.stop_reactivity.max(base.mutation_reactivity);
                // Only include positions with signals (reactivity > 0)
                if reactivity > 1e-6 {
                    react_vec[idx] = reactivity;
                }
            }
        }
        react_vec
    }).collect();

    // Build dendrogram (hierarchical clustering tree) for merge order
    let dendrogram = build_dendrogram(&reactivity_matrix);

    // Merge ensembles based on dendrogram order (from closest to farthest)
    let mut uf = UnionFind::new(first_ensembles.len());
    let mut merge_count = 0usize;
    let mut skipped_low_correlation = 0usize;

    // Process merges in order of similarity (from dendrogram)
    for merge_pair in &dendrogram {
        let (idx1, idx2, _distance) = *merge_pair;
        
        // Skip if already merged
        if uf.find(idx1) == uf.find(idx2) {
            continue;
        }

        let ensemble1 = &first_ensembles[idx1];
        let ensemble2 = &first_ensembles[idx2];

        // Get positions where BOTH ensembles have signals
        let mut signal_positions_1: HashSet<u32> = HashSet::new();
        let mut signal_positions_2: HashSet<u32> = HashSet::new();
        
        for base in first_base_records.iter().filter(|b| b.ensemble_id == ensemble1.id) {
            let reactivity = base.stop_reactivity.max(base.mutation_reactivity);
            if reactivity > 1e-6 {
                signal_positions_1.insert(base.position);
            }
        }
        
        for base in first_base_records.iter().filter(|b| b.ensemble_id == ensemble2.id) {
            let reactivity = base.stop_reactivity.max(base.mutation_reactivity);
            if reactivity > 1e-6 {
                signal_positions_2.insert(base.position);
            }
        }

        // Get shared signal positions (both have signals)
        let shared_signal_positions: Vec<u32> = signal_positions_1
            .intersection(&signal_positions_2)
            .copied()
            .collect();

        let shared_count = shared_signal_positions.len();
        let min_signal_positions = signal_positions_1.len().min(signal_positions_2.len());
        let max_signal_positions = signal_positions_1.len().max(signal_positions_2.len());

        // Need at least 3 shared signal positions to merge
        if shared_count < 3 {
            continue;
        }

        // Check overlap fraction thresholds (similar to first aggregation)
        let meets_fraction_threshold = if min_signal_positions > 0 {
            (shared_count as f64 / min_signal_positions as f64) >= min_shared_fraction
        } else {
            false
        };
        let meets_max_fraction = if max_signal_positions > 0 {
            (shared_count as f64 / max_signal_positions as f64) >= min_max_fraction
        } else {
            false
        };

        // Must meet at least one fraction threshold
        if !meets_fraction_threshold && !meets_max_fraction {
            continue;
        }

        // Calculate reactivity correlation on shared signal positions
        let correlation = calculate_ensemble_signal_correlation(
            ensemble1,
            ensemble2,
            &shared_signal_positions,
            first_base_records,
        );

        // Merge if correlation is high enough AND meets fraction thresholds
        if correlation >= min_reactivity_correlation && (meets_fraction_threshold || meets_max_fraction) {
            uf.union(idx1, idx2);
            merge_count += 1;
        } else {
            skipped_low_correlation += 1;
        }
    }

    eprintln!("[Duet Debug 2nd] Ensemble pairs in dendrogram: {}, merged: {}, skipped (low correlation): {}", 
              dendrogram.len(), merge_count, skipped_low_correlation);

    // Build final ensembles from merged components
    let mut components: HashMap<usize, Vec<usize>> = HashMap::new();
    for idx in 0..first_ensembles.len() {
        let root = uf.find(idx);
        components.entry(root).or_default().push(idx);
    }

    let mut global_ensembles: Vec<GlobalEnsemble> = Vec::new();
    let mut first_to_global: HashMap<usize, usize> = HashMap::new();
    let mut global_bases: Vec<GlobalBaseRecord> = Vec::new();

    let mut next_id = 1usize;
    for indices in components.values() {
        let first_ensemble = &first_ensembles[indices[0]];
        let mut window_set: HashSet<usize> = HashSet::new();
        let mut read_total = 0usize;
        let mut dual_total = 0usize;
        let mut base_map: HashMap<u32, PositionContribution> = HashMap::new();
        let mut position_confidence_sum: HashMap<u32, (f64, f64)> = HashMap::new();
        let mut ensemble_conf_sum = 0.0;
        let mut ensemble_conf_weight = 0.0;

        // Collect all positions from all merged ensembles
        let mut all_merged_positions: HashSet<u32> = HashSet::new();
        for &idx in indices {
            let ensemble = &first_ensembles[idx];
            first_to_global.insert(ensemble.id, next_id);
            window_set.insert(ensemble.window_count);
            read_total += ensemble.read_count;
            dual_total += ensemble.dual_signal_reads;
            ensemble_conf_sum += ensemble.confidence * ensemble.read_count as f64;
            ensemble_conf_weight += ensemble.read_count as f64;

            for base in first_base_records.iter().filter(|b| b.ensemble_id == ensemble.id) {
                all_merged_positions.insert(base.position);
                let contrib = base_map.entry(base.position).or_insert(PositionContribution {
                    stop_reads: 0,
                    mutation_reads: 0,
                });
                contrib.stop_reads += base.stop_reads;
                contrib.mutation_reads += base.mutation_reads;

                let conf_entry = position_confidence_sum.entry(base.position).or_insert((0.0, 0.0));
                let weight = (base.stop_reads + base.mutation_reads) as f64;
                conf_entry.0 += base.position_confidence * weight;
                conf_entry.1 += weight;
            }
        }

        let unique_positions = all_merged_positions.len();
        let ensemble_confidence = if ensemble_conf_weight > 0.0 {
            ensemble_conf_sum / ensemble_conf_weight
        } else {
            0.0
        };

        global_ensembles.push(GlobalEnsemble {
            id: next_id,
            chr: first_ensemble.chr.clone(),
            strand: first_ensemble.strand.clone(),
            window_count: window_set.len(),
            cluster_count: indices.len(),
            read_count: read_total,
            total_stop_reads: base_map.values().map(|c| c.stop_reads).sum(),
            total_mutation_reads: base_map.values().map(|c| c.mutation_reads).sum(),
            dual_signal_reads: dual_total,
            unique_positions,
            confidence: ensemble_confidence,
        });

        // Include all positions from merged ensembles
        for &pos in &all_merged_positions {
            let contrib = base_map.get(&pos).cloned().unwrap_or(PositionContribution {
                stop_reads: 0,
                mutation_reads: 0,
            });
            
            let key = format!("{}|{}|{}", first_ensemble.chr, first_ensemble.strand, pos);
            
            // Calculate weighted average reactivity from all merged ensembles
            let mut reactivity_sum_stop = 0.0;
            let mut reactivity_sum_mutation = 0.0;
            let mut weight_stop = 0.0;
            let mut weight_mutation = 0.0;
            
            for &idx in indices {
                let ensemble = &first_ensembles[idx];
                if let Some(base) = first_base_records.iter()
                    .find(|b| b.ensemble_id == ensemble.id && b.position == pos) {
                    if base.stop_reads > 0 {
                        reactivity_sum_stop += base.stop_reactivity * base.stop_reads as f64;
                        weight_stop += base.stop_reads as f64;
                    }
                    if base.mutation_reads > 0 {
                        reactivity_sum_mutation += base.mutation_reactivity * base.mutation_reads as f64;
                        weight_mutation += base.mutation_reads as f64;
                    }
                }
            }
            
            let stop_reactivity = if weight_stop > 0.0 {
                reactivity_sum_stop / weight_stop
            } else {
                positions.get(&key)
                    .map(|r| r.stop_reactivity)
                    .unwrap_or(0.0)
            };
            
            let mutation_reactivity = if weight_mutation > 0.0 {
                reactivity_sum_mutation / weight_mutation
            } else {
                positions.get(&key)
                    .map(|r| r.mutation_reactivity)
                    .unwrap_or(0.0)
            };
            
            let confidence = position_confidence_sum
                .get(&pos)
                .and_then(|(sum, weight)| {
                    if *weight > 0.0 {
                        Some(sum / weight)
                    } else {
                        None
                    }
                })
                .unwrap_or(0.0);
            
            global_bases.push(GlobalBaseRecord {
                ensemble_id: next_id,
                chr: first_ensemble.chr.clone(),
                strand: first_ensemble.strand.clone(),
                position: pos,
                stop_reads: contrib.stop_reads,
                mutation_reads: contrib.mutation_reads,
                stop_reactivity,
                mutation_reactivity,
                position_confidence: confidence,
            });
        }

        next_id += 1;
    }

    global_bases.sort_by(|a, b| {
        a.ensemble_id
            .cmp(&b.ensemble_id)
            .then(a.position.cmp(&b.position))
    });

    (global_ensembles, first_to_global, global_bases)
}

// Third aggregation: merge small fragment ensembles into larger ensembles
// Small ensembles are identified by having few signal positions or poor continuity
// They are merged into the ensemble with highest correlation among top 2 overlapping ensembles
fn aggregate_small_ensembles(
    second_ensembles: &[GlobalEnsemble],
    second_base_records: &[GlobalBaseRecord],
    positions: &HashMap<String, ReactivityEntry>,
    max_signal_positions: usize,
    min_correlation: f64,
) -> (
    Vec<GlobalEnsemble>,
    HashMap<usize, usize>, // second_id -> global_id
    Vec<GlobalBaseRecord>,
) {
    if second_ensembles.is_empty() {
        return (Vec::new(), HashMap::new(), Vec::new());
    }

    // Identify small ensembles: those with few signal positions
    let mut small_ensemble_indices: Vec<usize> = Vec::new();
    let mut large_ensemble_indices: Vec<usize> = Vec::new();
    
    for (idx, ensemble) in second_ensembles.iter().enumerate() {
        // Count signal positions (positions with reactivity > 0)
        let signal_count = second_base_records.iter()
            .filter(|b| b.ensemble_id == ensemble.id)
            .filter(|b| b.stop_reactivity.max(b.mutation_reactivity) > 1e-6)
            .count();
        
        if signal_count <= max_signal_positions {
            small_ensemble_indices.push(idx);
        } else {
            large_ensemble_indices.push(idx);
        }
    }

    eprintln!("[Duet Debug 3rd] Small ensembles: {} (signal positions <= {}), Large ensembles: {}", 
              small_ensemble_indices.len(), max_signal_positions, large_ensemble_indices.len());

    if small_ensemble_indices.is_empty() {
        // No small ensembles to merge, return as-is
        let mut second_to_global: HashMap<usize, usize> = HashMap::new();
        let mut global_ensembles: Vec<GlobalEnsemble> = Vec::new();
        let mut global_bases: Vec<GlobalBaseRecord> = Vec::new();
        
        let mut next_id = 1usize;
        for ensemble in second_ensembles {
            second_to_global.insert(ensemble.id, next_id);
            global_ensembles.push(GlobalEnsemble {
                id: next_id,
                chr: ensemble.chr.clone(),
                strand: ensemble.strand.clone(),
                window_count: ensemble.window_count,
                cluster_count: ensemble.cluster_count,
                read_count: ensemble.read_count,
                total_stop_reads: ensemble.total_stop_reads,
                total_mutation_reads: ensemble.total_mutation_reads,
                dual_signal_reads: ensemble.dual_signal_reads,
                unique_positions: ensemble.unique_positions,
                confidence: ensemble.confidence,
            });
            next_id += 1;
        }
        
        for base in second_base_records {
            if let Some(&global_id) = second_to_global.get(&base.ensemble_id) {
                global_bases.push(GlobalBaseRecord {
                    ensemble_id: global_id,
                    chr: base.chr.clone(),
                    strand: base.strand.clone(),
                    position: base.position,
                    stop_reads: base.stop_reads,
                    mutation_reads: base.mutation_reads,
                    stop_reactivity: base.stop_reactivity,
                    mutation_reactivity: base.mutation_reactivity,
                    position_confidence: base.position_confidence,
                });
            }
        }
        
        global_bases.sort_by(|a, b| {
            a.ensemble_id.cmp(&b.ensemble_id).then(a.position.cmp(&b.position))
        });
        
        return (global_ensembles, second_to_global, global_bases);
    }

    // Build position sets for each ensemble
    let mut ensemble_positions: Vec<HashSet<u32>> = Vec::new();
    let mut ensemble_signal_positions: Vec<HashSet<u32>> = Vec::new();
    
    for ensemble in second_ensembles {
        let mut pos_set = HashSet::new();
        let mut signal_set = HashSet::new();
        for base in second_base_records.iter().filter(|b| b.ensemble_id == ensemble.id) {
            pos_set.insert(base.position);
            let reactivity = base.stop_reactivity.max(base.mutation_reactivity);
            if reactivity > 1e-6 {
                signal_set.insert(base.position);
            }
        }
        ensemble_positions.push(pos_set);
        ensemble_signal_positions.push(signal_set);
    }

    // For each small ensemble, find top 2 overlapping large ensembles
    let mut small_to_target: HashMap<usize, usize> = HashMap::new(); // small_idx -> target_idx
    let mut merge_count = 0usize;
    let mut skipped_low_correlation = 0usize;

    for &small_idx in &small_ensemble_indices {
        let small_positions = &ensemble_positions[small_idx];
        let small_signal_positions = &ensemble_signal_positions[small_idx];
        
        // Find overlaps with large ensembles
        let mut overlaps: Vec<(usize, usize, f64)> = Vec::new(); // (large_idx, overlap_count, correlation)
        
        for &large_idx in &large_ensemble_indices {
            let large_positions = &ensemble_positions[large_idx];
            let large_signal_positions = &ensemble_signal_positions[large_idx];
            
            // Calculate overlap (shared positions)
            let overlap_count = small_positions.intersection(large_positions).count();
            
            if overlap_count > 0 {
                // Get shared signal positions
                let shared_signal_positions: Vec<u32> = small_signal_positions
                    .intersection(large_signal_positions)
                    .copied()
                    .collect();
                
                // Calculate correlation on shared signal positions
                let correlation = if shared_signal_positions.len() >= 3 {
                    calculate_ensemble_signal_correlation(
                        &second_ensembles[small_idx],
                        &second_ensembles[large_idx],
                        &shared_signal_positions,
                        second_base_records,
                    )
                } else if shared_signal_positions.len() > 0 {
                    // If fewer than 3 shared signal positions, use a lower threshold
                    // but still calculate correlation
                    calculate_ensemble_signal_correlation(
                        &second_ensembles[small_idx],
                        &second_ensembles[large_idx],
                        &shared_signal_positions,
                        second_base_records,
                    )
                } else {
                    0.0
                };
                
                overlaps.push((large_idx, overlap_count, correlation));
            }
        }
        
        // Sort by overlap count (descending), then by correlation (descending)
        overlaps.sort_by(|a, b| {
            b.1.cmp(&a.1).then(b.2.partial_cmp(&a.2).unwrap_or(std::cmp::Ordering::Equal))
        });
        
        // Take top 2 candidates
        if overlaps.len() >= 1 {
            let (target_idx, _overlap, correlation) = overlaps[0];
            
            // If we have 2 candidates, compare correlations
            if overlaps.len() >= 2 {
                let (alt_idx, _alt_overlap, alt_correlation) = overlaps[1];
                
                // Choose the one with higher correlation
                if alt_correlation > correlation {
                    if alt_correlation >= min_correlation {
                        small_to_target.insert(small_idx, alt_idx);
                        merge_count += 1;
                    } else {
                        skipped_low_correlation += 1;
                    }
                } else {
                    if correlation >= min_correlation {
                        small_to_target.insert(small_idx, target_idx);
                        merge_count += 1;
                    } else {
                        skipped_low_correlation += 1;
                    }
                }
            } else {
                // Only one candidate
                if correlation >= min_correlation {
                    small_to_target.insert(small_idx, target_idx);
                    merge_count += 1;
                } else {
                    skipped_low_correlation += 1;
                }
            }
        }
    }

    eprintln!("[Duet Debug 3rd] Small ensembles processed: {}, merged: {}, skipped (low correlation): {}", 
              small_ensemble_indices.len(), merge_count, skipped_low_correlation);

    // Build final ensembles: merge small ensembles into their targets
    let mut uf = UnionFind::new(second_ensembles.len());
    
    // First, merge small ensembles into their targets
    for (small_idx, &target_idx) in &small_to_target {
        uf.union(*small_idx, target_idx);
    }
    
    // Build components
    let mut components: HashMap<usize, Vec<usize>> = HashMap::new();
    for idx in 0..second_ensembles.len() {
        let root = uf.find(idx);
        components.entry(root).or_default().push(idx);
    }

    let mut global_ensembles: Vec<GlobalEnsemble> = Vec::new();
    let mut second_to_global: HashMap<usize, usize> = HashMap::new();
    let mut global_bases: Vec<GlobalBaseRecord> = Vec::new();

    let mut next_id = 1usize;
    for indices in components.values() {
        let second_ensemble = &second_ensembles[indices[0]];
        let mut window_set: HashSet<usize> = HashSet::new();
        let mut read_total = 0usize;
        let mut dual_total = 0usize;
        let mut base_map: HashMap<u32, PositionContribution> = HashMap::new();
        let mut position_confidence_sum: HashMap<u32, (f64, f64)> = HashMap::new();
        let mut ensemble_conf_sum = 0.0;
        let mut ensemble_conf_weight = 0.0;

        // Collect all positions from all merged ensembles
        let mut all_merged_positions: HashSet<u32> = HashSet::new();
        for &idx in indices {
            let ensemble = &second_ensembles[idx];
            second_to_global.insert(ensemble.id, next_id);
            window_set.insert(ensemble.window_count);
            read_total += ensemble.read_count;
            dual_total += ensemble.dual_signal_reads;
            ensemble_conf_sum += ensemble.confidence * ensemble.read_count as f64;
            ensemble_conf_weight += ensemble.read_count as f64;

            for base in second_base_records.iter().filter(|b| b.ensemble_id == ensemble.id) {
                all_merged_positions.insert(base.position);
                let contrib = base_map.entry(base.position).or_insert(PositionContribution {
                    stop_reads: 0,
                    mutation_reads: 0,
                });
                contrib.stop_reads += base.stop_reads;
                contrib.mutation_reads += base.mutation_reads;

                let conf_entry = position_confidence_sum.entry(base.position).or_insert((0.0, 0.0));
                let weight = (base.stop_reads + base.mutation_reads) as f64;
                conf_entry.0 += base.position_confidence * weight;
                conf_entry.1 += weight;
            }
        }

        let unique_positions = all_merged_positions.len();
        let ensemble_confidence = if ensemble_conf_weight > 0.0 {
            ensemble_conf_sum / ensemble_conf_weight
        } else {
            0.0
        };

        global_ensembles.push(GlobalEnsemble {
            id: next_id,
            chr: second_ensemble.chr.clone(),
            strand: second_ensemble.strand.clone(),
            window_count: window_set.len(),
            cluster_count: indices.len(),
            read_count: read_total,
            total_stop_reads: base_map.values().map(|c| c.stop_reads).sum(),
            total_mutation_reads: base_map.values().map(|c| c.mutation_reads).sum(),
            dual_signal_reads: dual_total,
            unique_positions,
            confidence: ensemble_confidence,
        });

        // Include all positions from merged ensembles
        for &pos in &all_merged_positions {
            let contrib = base_map.get(&pos).cloned().unwrap_or(PositionContribution {
                stop_reads: 0,
                mutation_reads: 0,
            });
            
            let key = format!("{}|{}|{}", second_ensemble.chr, second_ensemble.strand, pos);
            
            // Calculate weighted average reactivity from all merged ensembles
            let mut reactivity_sum_stop = 0.0;
            let mut reactivity_sum_mutation = 0.0;
            let mut weight_stop = 0.0;
            let mut weight_mutation = 0.0;
            
            for &idx in indices {
                let ensemble = &second_ensembles[idx];
                if let Some(base) = second_base_records.iter()
                    .find(|b| b.ensemble_id == ensemble.id && b.position == pos) {
                    if base.stop_reads > 0 {
                        reactivity_sum_stop += base.stop_reactivity * base.stop_reads as f64;
                        weight_stop += base.stop_reads as f64;
                    }
                    if base.mutation_reads > 0 {
                        reactivity_sum_mutation += base.mutation_reactivity * base.mutation_reads as f64;
                        weight_mutation += base.mutation_reads as f64;
                    }
                }
            }
            
            let stop_reactivity = if weight_stop > 0.0 {
                reactivity_sum_stop / weight_stop
            } else {
                positions.get(&key)
                    .map(|r| r.stop_reactivity)
                    .unwrap_or(0.0)
            };
            
            let mutation_reactivity = if weight_mutation > 0.0 {
                reactivity_sum_mutation / weight_mutation
            } else {
                positions.get(&key)
                    .map(|r| r.mutation_reactivity)
                    .unwrap_or(0.0)
            };
            
            let confidence = position_confidence_sum
                .get(&pos)
                .and_then(|(sum, weight)| {
                    if *weight > 0.0 {
                        Some(sum / weight)
                    } else {
                        None
                    }
                })
                .unwrap_or(0.0);
            
            global_bases.push(GlobalBaseRecord {
                ensemble_id: next_id,
                chr: second_ensemble.chr.clone(),
                strand: second_ensemble.strand.clone(),
                position: pos,
                stop_reads: contrib.stop_reads,
                mutation_reads: contrib.mutation_reads,
                stop_reactivity,
                mutation_reactivity,
                position_confidence: confidence,
            });
        }

        next_id += 1;
    }

    global_bases.sort_by(|a, b| {
        a.ensemble_id
            .cmp(&b.ensemble_id)
            .then(a.position.cmp(&b.position))
    });

    (global_ensembles, second_to_global, global_bases)
}

// Build dendrogram (hierarchical clustering tree) for guiding merge order
// Returns list of (idx1, idx2, distance) pairs in order of increasing distance
fn build_dendrogram(matrix: &[Vec<f64>]) -> Vec<(usize, usize, f64)> {
    if matrix.is_empty() {
        return Vec::new();
    }
    
    let n = matrix.len();
    if n == 1 {
        return Vec::new();
    }
    
    // Compute pairwise distances (Euclidean distance on signal positions only)
    let mut distances: Vec<Vec<f64>> = vec![vec![0.0; n]; n];
    for i in 0..n {
        for j in (i + 1)..n {
            // Only consider positions where both have signals
            let mut signal_positions = Vec::new();
            for (idx, (&val1, &val2)) in matrix[i].iter().zip(matrix[j].iter()).enumerate() {
                if val1 > 1e-6 && val2 > 1e-6 {
                    signal_positions.push((idx, val1, val2));
                }
            }
            
            if signal_positions.is_empty() {
                distances[i][j] = f64::INFINITY;
                distances[j][i] = f64::INFINITY;
            } else {
                let dist: f64 = signal_positions.iter()
                    .map(|(_, a, b)| (a - b) * (a - b))
                    .sum::<f64>()
                    .sqrt();
                distances[i][j] = dist;
                distances[j][i] = dist;
            }
        }
    }
    
    // Build dendrogram using hierarchical clustering
    let mut dendrogram: Vec<(usize, usize, f64)> = Vec::new();
    let mut clusters: Vec<Vec<usize>> = (0..n).map(|i| vec![i]).collect();
    let mut active: Vec<bool> = vec![true; n];
    
    // Merge clusters until we have a single cluster
    while clusters.iter().filter(|c| active[c[0]]).count() > 1 {
        // Find the two closest active clusters
        let mut min_dist = f64::INFINITY;
        let mut merge_i = 0;
        let mut merge_j = 1;
        let mut merge_idx1 = 0;
        let mut merge_idx2 = 1;
        
        let active_clusters: Vec<usize> = (0..clusters.len())
            .filter(|&i| active[i])
            .collect();
        
        for &i in &active_clusters {
            for &j in &active_clusters {
                if i >= j {
                    continue;
                }
                
                // Compute average linkage distance between clusters
                let mut total_dist = 0.0;
                let mut count = 0;
                for &idx1 in &clusters[i] {
                    for &idx2 in &clusters[j] {
                        if distances[idx1][idx2] < f64::INFINITY {
                            total_dist += distances[idx1][idx2];
                            count += 1;
                        }
                    }
                }
                let avg_dist = if count > 0 { total_dist / count as f64 } else { f64::INFINITY };
                
                if avg_dist < min_dist {
                    min_dist = avg_dist;
                    merge_i = i;
                    merge_j = j;
                    merge_idx1 = clusters[i][0];
                    merge_idx2 = clusters[j][0];
                }
            }
        }
        
        if min_dist == f64::INFINITY {
            break; // No more merges possible
        }
        
        // Record this merge in dendrogram
        dendrogram.push((merge_idx1, merge_idx2, min_dist));
        
        // Merge clusters
        let mut merged = clusters[merge_i].clone();
        merged.extend_from_slice(&clusters[merge_j]);
        clusters[merge_i] = merged;
        active[merge_j] = false;
    }
    
    dendrogram
}

// Calculate reactivity correlation between two ensembles on shared signal positions only
fn calculate_ensemble_signal_correlation(
    ensemble1: &GlobalEnsemble,
    ensemble2: &GlobalEnsemble,
    shared_signal_positions: &[u32],
    base_records: &[GlobalBaseRecord],
) -> f64 {
    if shared_signal_positions.len() < 3 {
        return 0.0;
    }

    let mut react1 = Vec::new();
    let mut react2 = Vec::new();

    for &pos in shared_signal_positions {
        let reactivity1 = base_records.iter()
            .find(|b| b.ensemble_id == ensemble1.id && b.position == pos)
            .map(|b| b.stop_reactivity.max(b.mutation_reactivity))
            .unwrap_or(0.0);

        let reactivity2 = base_records.iter()
            .find(|b| b.ensemble_id == ensemble2.id && b.position == pos)
            .map(|b| b.stop_reactivity.max(b.mutation_reactivity))
            .unwrap_or(0.0);

        // Both should have signals (already filtered, but double-check)
        if reactivity1 > 1e-6 && reactivity2 > 1e-6 {
            react1.push(reactivity1);
            react2.push(reactivity2);
        }
    }

    if react1.len() < 3 {
        return 0.0;
    }

    // Calculate Pearson correlation
    let mean1 = react1.iter().sum::<f64>() / react1.len() as f64;
    let mean2 = react2.iter().sum::<f64>() / react2.len() as f64;

    let mut numerator = 0.0;
    let mut var1 = 0.0;
    let mut var2 = 0.0;

    for i in 0..react1.len() {
        let diff1 = react1[i] - mean1;
        let diff2 = react2[i] - mean2;
        numerator += diff1 * diff2;
        var1 += diff1 * diff1;
        var2 += diff2 * diff2;
    }

    if var1 > 1e-10 && var2 > 1e-10 {
        numerator / (var1.sqrt() * var2.sqrt())
    } else {
        0.0
    }
}

fn plot_global_ensemble_heatmap(
    global_ensembles: &[GlobalEnsemble],
    global_bases: &[GlobalBaseRecord],
    output_path: &str,
) -> Result<(), Box<dyn Error>> {
    if global_ensembles.is_empty() || global_bases.is_empty() {
        return Err("No global ensemble data to plot".into());
    }

    // Group bases by ensemble_id and collect positions
    let mut ensemble_data: HashMap<usize, Vec<&GlobalBaseRecord>> = HashMap::new();
    for base in global_bases {
        ensemble_data.entry(base.ensemble_id).or_default().push(base);
    }

    // Sort ensembles by ID
    let mut ensemble_ids: Vec<usize> = ensemble_data.keys().copied().collect();
    ensemble_ids.sort();

    if ensemble_ids.is_empty() {
        return Err("No ensemble data found".into());
    }

    // Collect all unique positions across all ensembles
    let mut all_positions: Vec<u32> = global_bases
        .iter()
        .map(|b| b.position)
        .collect::<HashSet<_>>()
        .into_iter()
        .collect();
    all_positions.sort();

    if all_positions.is_empty() {
        return Err("No position data found".into());
    }

    // Build reactivity matrix: rows = ensembles, columns = positions
    // Use combined reactivity (max of stop and mutation reactivity)
    let num_ensembles = ensemble_ids.len();
    let num_positions = all_positions.len();
    
    let mut reactivity_matrix: Vec<Vec<f64>> = vec![vec![0.0; num_positions]; num_ensembles];
    let mut position_to_idx: HashMap<u32, usize> = HashMap::new();
    for (idx, &pos) in all_positions.iter().enumerate() {
        position_to_idx.insert(pos, idx);
    }

    for (row_idx, &ensemble_id) in ensemble_ids.iter().enumerate() {
        if let Some(bases) = ensemble_data.get(&ensemble_id) {
            for base in bases {
                if let Some(&col_idx) = position_to_idx.get(&base.position) {
                    // Use combined reactivity (max of stop and mutation)
                    let reactivity = base.stop_reactivity.max(base.mutation_reactivity);
                    reactivity_matrix[row_idx][col_idx] = reactivity;
                }
            }
        }
    }

    // Perform hierarchical clustering on rows (ensembles) to reorder them
    let row_order = hierarchical_cluster_rows(&reactivity_matrix);
    
    // Reorder ensemble_ids and reactivity_matrix according to clustering
    let mut reordered_ensemble_ids: Vec<usize> = row_order.iter().map(|&idx| ensemble_ids[idx]).collect();
    let mut reordered_matrix: Vec<Vec<f64>> = row_order.iter().map(|&idx| reactivity_matrix[idx].clone()).collect();
    
    ensemble_ids = reordered_ensemble_ids;
    reactivity_matrix = reordered_matrix;

    // Find max reactivity for color scaling
    let max_reactivity = reactivity_matrix
        .iter()
        .flat_map(|row| row.iter())
        .fold(0.0f64, |a, &b| a.max(b));

    if max_reactivity <= 0.0 {
        return Err("No reactivity data found".into());
    }

    // Calculate dimensions dynamically based on number of ensembles
    let title_height = 60;
    let label_width = 80;
    let label_height = 30;
    let colorbar_width = 30;
    let colorbar_margin = 20;
    let plot_margin = 20;
    
    // Minimum cell height for better visibility (at least 15 pixels per row)
    let min_cell_height = 15.0;
    let base_plot_width = 1200;
    let base_plot_height = 800;
    
    // Calculate required height based on number of ensembles
    let required_plot_height = (num_ensembles as f64 * min_cell_height).max(400.0);
    let total_height = (title_height + label_height + plot_margin * 2 + required_plot_height as i32) as u32;
    
    // Use fixed width, dynamic height
    let plot_width = base_plot_width - label_width - colorbar_width - colorbar_margin - plot_margin * 2;
    let plot_height = required_plot_height;
    
    let cell_width = plot_width as f64 / num_positions as f64;
    let cell_height = plot_height / num_ensembles as f64;
    
    // Create SVG with dynamic height
    let root = SVGBackend::new(output_path, (base_plot_width as u32, total_height)).into_drawing_area();
    root.fill(&WHITE)?;

    // Draw title
    root.draw(&Text::new(
        "Global Ensemble Reactivity Distribution",
        (600, 30),
        ("sans-serif", 30).into_font().color(&BLACK),
    ))?;

    // Calculate signal position counts for each ensemble
    let mut ensemble_signal_counts: HashMap<usize, usize> = HashMap::new();
    for &ensemble_id in &ensemble_ids {
        let signal_count = ensemble_data.get(&ensemble_id)
            .map(|bases| {
                bases.iter()
                    .filter(|b| b.stop_reactivity.max(b.mutation_reactivity) > 1e-6)
                    .count()
            })
            .unwrap_or(0);
        ensemble_signal_counts.insert(ensemble_id, signal_count);
    }

    // Draw ensemble labels (y-axis) with signal position counts
    for (row_idx, &ensemble_id) in ensemble_ids.iter().enumerate() {
        let y = title_height + plot_margin + ((row_idx as f64 + 0.5) * cell_height) as i32;
        let signal_count = ensemble_signal_counts.get(&ensemble_id).copied().unwrap_or(0);
        root.draw(&Text::new(
            format!("E{} ({})", ensemble_id, signal_count),
            (label_width / 2, y),
            ("sans-serif", 12).into_font().color(&BLACK),
        ))?;
    }

    // Draw position labels (x-axis) - sample every 10th position to avoid crowding
    let label_interval = (num_positions / 20).max(1);
    for (col_idx, &pos) in all_positions.iter().enumerate() {
        if col_idx % label_interval == 0 || col_idx == num_positions - 1 {
            let x = label_width + plot_margin + ((col_idx as f64 + 0.5) * cell_width) as i32;
            root.draw(&Text::new(
                pos.to_string(),
                (x, (title_height as f64 + plot_height + plot_margin as f64) as i32 + 15),
                ("sans-serif", 10).into_font().color(&BLACK),
            ))?;
        }
    }

    // Draw heatmap cells - draw all cells, even with zero reactivity
    for (row_idx, row) in reactivity_matrix.iter().enumerate() {
        for (col_idx, &reactivity) in row.iter().enumerate() {
            let x_start = (label_width + plot_margin) as f64 + (col_idx as f64 * cell_width);
            let y_start = (title_height + plot_margin) as f64 + (row_idx as f64 * cell_height);
            let x_end = x_start + cell_width;
            let y_end = y_start + cell_height;
            
            // Normalize reactivity to 0-1 for color mapping
            let normalized = if max_reactivity > 0.0 {
                (reactivity / max_reactivity).min(1.0).max(0.0)
            } else {
                0.0
            };
            
            // Use a color gradient: white (0) -> yellow -> orange -> red (1)
            // This provides better visibility than blue-white-red
            let color = if normalized < 0.33 {
                // White to yellow: (255, 255, 255) -> (255, 255, 0)
                let t = normalized / 0.33;
                RGBColor(
                    255,
                    255,
                    ((1.0 - t) * 255.0) as u8,
                )
            } else if normalized < 0.67 {
                // Yellow to orange: (255, 255, 0) -> (255, 165, 0)
                let t = (normalized - 0.33) / 0.34;
                RGBColor(
                    255,
                    (255.0 - t * 90.0) as u8,
                    0,
                )
            } else {
                // Orange to red: (255, 165, 0) -> (255, 0, 0)
                let t = (normalized - 0.67) / 0.33;
                RGBColor(
                    255,
                    ((165.0 * (1.0 - t)) as u8),
                    0,
                )
            };

            root.draw(&Rectangle::new(
                [(x_start as i32, y_start as i32), (x_end as i32, y_end as i32)],
                color.filled(),
            ))?;
        }
    }

    // Draw colorbar
    let colorbar_x = base_plot_width - colorbar_width - colorbar_margin;
    let colorbar_y = title_height + plot_margin;
    let colorbar_height_f = plot_height;
    let colorbar_step = colorbar_height_f / 100.0;

    for i in 0..100 {
        let y_start = colorbar_y as f64 + i as f64 * colorbar_step;
        let y_end = y_start + colorbar_step;
        let normalized = 1.0 - (i as f64 / 100.0);
        let color = if normalized < 0.33 {
            // White to yellow
            let t = normalized / 0.33;
            RGBColor(
                255,
                255,
                ((1.0 - t) * 255.0) as u8,
            )
        } else if normalized < 0.67 {
            // Yellow to orange
            let t = (normalized - 0.33) / 0.34;
            RGBColor(
                255,
                (255.0 - t * 90.0) as u8,
                0,
            )
        } else {
            // Orange to red
            let t = (normalized - 0.67) / 0.33;
            RGBColor(
                255,
                ((165.0 * (1.0 - t)) as u8),
                0,
            )
        };
        root.draw(&Rectangle::new(
            [(colorbar_x as i32, y_start as i32), ((colorbar_x + colorbar_width) as i32, y_end as i32)],
            color.filled(),
        ))?;
    }

    // Draw colorbar labels
    root.draw(&Text::new(
        format!("{:.2}", max_reactivity),
        (colorbar_x as i32 + colorbar_width / 2, colorbar_y as i32 - 5),
        ("sans-serif", 10).into_font().color(&BLACK),
    ))?;
    root.draw(&Text::new(
        "0.00",
        (colorbar_x as i32 + colorbar_width / 2, (colorbar_y as f64 + plot_height) as i32 + 15),
        ("sans-serif", 10).into_font().color(&BLACK),
    ))?;
    root.draw(&Text::new(
        "Reactivity",
        (colorbar_x as i32 - 20, (colorbar_y as f64 + plot_height / 2.0) as i32),
        ("sans-serif", 12).into_font().color(&BLACK),
    ))?;

    root.present()?;

    Ok(())
}

fn plot_global_ensemble_tsne(
    global_ensembles: &[GlobalEnsemble],
    global_bases: &[GlobalBaseRecord],
    output_path: &str,
) -> Result<(), Box<dyn Error>> {
    if global_ensembles.is_empty() || global_bases.is_empty() {
        return Err("No global ensemble data to plot".into());
    }

    // Group bases by ensemble_id
    let mut ensemble_data: HashMap<usize, Vec<&GlobalBaseRecord>> = HashMap::new();
    for base in global_bases {
        ensemble_data.entry(base.ensemble_id).or_default().push(base);
    }

    // Sort ensembles by ID
    let mut ensemble_ids: Vec<usize> = ensemble_data.keys().copied().collect();
    ensemble_ids.sort();

    if ensemble_ids.is_empty() {
        return Err("No ensemble data found".into());
    }

    // Collect all unique positions across all ensembles
    let mut all_positions: Vec<u32> = global_bases
        .iter()
        .map(|b| b.position)
        .collect::<HashSet<_>>()
        .into_iter()
        .collect();
    all_positions.sort();

    if all_positions.is_empty() {
        return Err("No position data found".into());
    }

    // Build reactivity matrix: rows = ensembles, columns = positions
    let num_ensembles = ensemble_ids.len();
    let num_positions = all_positions.len();
    
    let mut reactivity_matrix: Vec<Vec<f64>> = vec![vec![0.0; num_positions]; num_ensembles];
    let mut position_to_idx: HashMap<u32, usize> = HashMap::new();
    for (idx, &pos) in all_positions.iter().enumerate() {
        position_to_idx.insert(pos, idx);
    }

    for (row_idx, &ensemble_id) in ensemble_ids.iter().enumerate() {
        if let Some(bases) = ensemble_data.get(&ensemble_id) {
            for base in bases {
                if let Some(&col_idx) = position_to_idx.get(&base.position) {
                    // Use combined reactivity (max of stop and mutation)
                    let reactivity = base.stop_reactivity.max(base.mutation_reactivity);
                    reactivity_matrix[row_idx][col_idx] = reactivity;
                }
            }
        }
    }

    // Count signal positions for each ensemble (for point size)
    let signal_counts: Vec<usize> = reactivity_matrix.iter().map(|row| {
        row.iter().filter(|&&val| val > 0.0).count()
    }).collect();

    // Perform UMAP (simplified version) for dimensionality reduction
    let (umap1, umap2) = perform_umap(&reactivity_matrix)?;

    // Find ranges for scaling
    let pc1_min = umap1.iter().fold(f64::INFINITY, |a, &b| a.min(b));
    let pc1_max = umap1.iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b));
    let pc2_min = umap2.iter().fold(f64::INFINITY, |a, &b| a.min(b));
    let pc2_max = umap2.iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b));

    let pc1_range = pc1_max - pc1_min;
    let pc2_range = pc2_max - pc2_min;
    let margin = pc1_range.max(pc2_range) * 0.1;

    // Create SVG plot
    let root = SVGBackend::new(output_path, (1000, 800)).into_drawing_area();
    root.fill(&WHITE)?;

    // Calculate plot area
    let title_height = 60;
    let margin_size = 80;
    let plot_width = 1000 - margin_size * 2;
    let plot_height = 800 - title_height - margin_size * 2;

    // Draw title
    root.draw(&Text::new(
        "Global Ensemble Reactivity Distribution (UMAP)",
        (500, 30),
        ("sans-serif", 28).into_font().color(&BLACK),
    ))?;

    // Draw axes using Rectangle (simpler approach)
    let x_axis_y = title_height + plot_height + margin_size;
    root.draw(&Rectangle::new(
        [(margin_size, x_axis_y - 1), (margin_size + plot_width, x_axis_y + 1)],
        BLACK.filled(),
    ))?;
    root.draw(&Text::new(
        "UMAP1",
        (margin_size + plot_width / 2, x_axis_y + 40),
        ("sans-serif", 16).into_font().color(&BLACK),
    ))?;

    let y_axis_x = margin_size;
    root.draw(&Rectangle::new(
        [(y_axis_x - 1, title_height), (y_axis_x + 1, title_height + plot_height)],
        BLACK.filled(),
    ))?;
    root.draw(&Text::new(
        "UMAP2",
        (y_axis_x - 50, title_height + plot_height / 2),
        ("sans-serif", 16).into_font().color(&BLACK),
    ))?;

    // Draw grid lines using thin rectangles
    for i in 0..6 {
        let x = margin_size + ((i as f64 / 5.0) * plot_width as f64) as i32;
        root.draw(&Rectangle::new(
            [(x, title_height), (x + 1, title_height + plot_height)],
            RGBColor(200, 200, 200).filled(),
        ))?;
    }
    for i in 0..6 {
        let y = title_height + ((i as f64 / 5.0) * plot_height as f64) as i32;
        root.draw(&Rectangle::new(
            [(margin_size, y), (margin_size + plot_width, y + 1)],
            RGBColor(200, 200, 200).filled(),
        ))?;
    }

    // Get start positions for each ensemble (for coloring)
    let mut ensemble_start_positions: Vec<u32> = Vec::new();
    for &ensemble_id in &ensemble_ids {
        if let Some(ensemble) = global_ensembles.iter().find(|e| e.id == ensemble_id) {
            // Find the minimum position for this ensemble
            let start_pos = ensemble_data.get(&ensemble_id)
                .and_then(|bases| bases.iter().map(|b| b.position).min())
                .unwrap_or(0);
            ensemble_start_positions.push(start_pos);
        } else {
            ensemble_start_positions.push(0);
        }
    }
    
    // Find min/max start positions for color mapping
    let min_start = ensemble_start_positions.iter().min().copied().unwrap_or(0);
    let max_start = ensemble_start_positions.iter().max().copied().unwrap_or(1);
    let start_range = (max_start - min_start).max(1) as f64;
    
    // Generate colors based on start position
    let colors: Vec<RGBColor> = ensemble_start_positions.iter().map(|&start_pos| {
        // Normalize start position to 0-1
        let normalized = if start_range > 0.0 {
            ((start_pos - min_start) as f64 / start_range).min(1.0).max(0.0)
        } else {
            0.5
        };
        
        // Use a color gradient: blue (early) -> green -> yellow -> red (late)
        let (r, g, b) = if normalized < 0.33 {
            // Blue to green
            let t = normalized / 0.33;
            (0, (t * 255.0) as u8, ((1.0 - t) * 255.0) as u8)
        } else if normalized < 0.67 {
            // Green to yellow
            let t = (normalized - 0.33) / 0.34;
            ((t * 255.0) as u8, 255, 0)
        } else {
            // Yellow to red
            let t = (normalized - 0.67) / 0.33;
            (255, ((1.0 - t) * 255.0) as u8, 0)
        };
        
        RGBColor(r, g, b)
    }).collect();

    // Find min/max signal counts for scaling point sizes
    let min_signals = signal_counts.iter().min().copied().unwrap_or(1);
    let max_signals = signal_counts.iter().max().copied().unwrap_or(1);
    let signal_range = (max_signals - min_signals) as f64;
    let min_point_size = 5.0;
    let max_point_size = 25.0;
    let point_size_range = max_point_size - min_point_size;

    // Draw points
    for (idx, &ensemble_id) in ensemble_ids.iter().enumerate() {
        let x_normalized = (umap1[idx] - pc1_min + margin) / (pc1_range + 2.0 * margin);
        let y_normalized = 1.0 - (umap2[idx] - pc2_min + margin) / (pc2_range + 2.0 * margin);
        
        let x = margin_size + (x_normalized * plot_width as f64) as i32;
        let y = title_height + (y_normalized * plot_height as f64) as i32;
        
        let color = &colors[idx];
        
        // Calculate point size based on signal count
        let signal_count = signal_counts[idx];
        let size_normalized = if signal_range > 0.0 {
            ((signal_count - min_signals) as f64 / signal_range).min(1.0).max(0.0)
        } else {
            0.5
        };
        let point_size = min_point_size + size_normalized * point_size_range;
        
        // Draw circle with size proportional to signal count
        root.draw(&Circle::new(
            (x, y),
            point_size as i32,
            color.filled(),
        ))?;
        
        // Draw label (only if point is large enough)
        if point_size > 10.0 {
            root.draw(&Text::new(
                format!("E{}", ensemble_id),
                (x + (point_size as i32) + 2, y + 4),
                ("sans-serif", 11).into_font().color(&BLACK),
            ))?;
        }
    }

    // Draw legend
    let legend_x = 1000 - 200;
    let legend_y = title_height + 50;
    let legend_width = 180;
    let legend_height = 120;
    
    // Draw legend background
    root.draw(&Rectangle::new(
        [(legend_x, legend_y), (legend_x + legend_width, legend_y + legend_height)],
        RGBColor(250, 250, 250).filled(),
    ))?;
    // Draw legend border
    root.draw(&Rectangle::new(
        [(legend_x, legend_y), (legend_x + legend_width, legend_y + legend_height)],
        BLACK.stroke_width(1),
    ))?;
    
    // Legend title
    root.draw(&Text::new(
        "Legend",
        (legend_x + 10, legend_y + 20),
        ("sans-serif", 14).into_font().color(&BLACK),
    ))?;
    
    // Color legend: Start Position
    root.draw(&Text::new(
        "Color: Start Position",
        (legend_x + 10, legend_y + 45),
        ("sans-serif", 11).into_font().color(&BLACK),
    ))?;
    
    // Draw color gradient bar
    let colorbar_x = legend_x + 10;
    let colorbar_y = legend_y + 55;
    let colorbar_w = 160;
    let colorbar_h = 10;
    for i in 0..colorbar_w {
        let normalized = i as f64 / colorbar_w as f64;
        let (r, g, b) = if normalized < 0.33 {
            let t = normalized / 0.33;
            (0, (t * 255.0) as u8, ((1.0 - t) * 255.0) as u8)
        } else if normalized < 0.67 {
            let t = (normalized - 0.33) / 0.34;
            ((t * 255.0) as u8, 255, 0)
        } else {
            let t = (normalized - 0.67) / 0.33;
            (255, ((1.0 - t) * 255.0) as u8, 0)
        };
        root.draw(&Rectangle::new(
            [(colorbar_x + i, colorbar_y), (colorbar_x + i + 1, colorbar_y + colorbar_h)],
            RGBColor(r, g, b).filled(),
        ))?;
    }
    
    // Color bar labels
    root.draw(&Text::new(
        format!("{}", min_start),
        (colorbar_x, colorbar_y + 20),
        ("sans-serif", 9).into_font().color(&BLACK),
    ))?;
    root.draw(&Text::new(
        format!("{}", max_start),
        (colorbar_x + colorbar_w - 30, colorbar_y + 20),
        ("sans-serif", 9).into_font().color(&BLACK),
    ))?;
    
    // Size legend: Signal Positions
    root.draw(&Text::new(
        "Size: Signal Positions",
        (legend_x + 10, legend_y + 85),
        ("sans-serif", 11).into_font().color(&BLACK),
    ))?;
    
    // Draw example circles for size
    let size_legend_x = legend_x + 10;
    let size_legend_y = legend_y + 95;
    root.draw(&Circle::new(
        (size_legend_x + 5, size_legend_y),
        min_point_size as i32,
        RGBColor(100, 100, 100).filled(),
    ))?;
    root.draw(&Text::new(
        format!("{}", min_signals),
        (size_legend_x + 15, size_legend_y + 4),
        ("sans-serif", 9).into_font().color(&BLACK),
    ))?;
    
    root.draw(&Circle::new(
        (size_legend_x + 60, size_legend_y),
        max_point_size as i32,
        RGBColor(100, 100, 100).filled(),
    ))?;
    root.draw(&Text::new(
        format!("{}", max_signals),
        (size_legend_x + 70, size_legend_y + 4),
        ("sans-serif", 9).into_font().color(&BLACK),
    ))?;

    root.present()?;
    Ok(())
}

// Plot overlap tracks: for each non-full-length ensemble, plot its overlap region with full-length ensembles
fn plot_ensemble_overlap_tracks(
    global_ensembles: &[GlobalEnsemble],
    global_bases: &[GlobalBaseRecord],
    output_dir: &str,
) -> Result<(), Box<dyn Error>> {
    if global_ensembles.is_empty() || global_bases.is_empty() {
        return Err("No global ensemble data to plot".into());
    }

    // Group bases by ensemble_id
    let mut ensemble_data: HashMap<usize, Vec<&GlobalBaseRecord>> = HashMap::new();
    for base in global_bases {
        ensemble_data.entry(base.ensemble_id).or_default().push(base);
    }

    // Find full-length ensembles (those with maximum position coverage)
    let mut ensemble_ranges: Vec<(usize, u32, u32)> = Vec::new();
    for ensemble in global_ensembles {
        if let Some(bases) = ensemble_data.get(&ensemble.id) {
            let positions: Vec<u32> = bases.iter().map(|b| b.position).collect();
            if let (Some(&min_pos), Some(&max_pos)) = (positions.iter().min(), positions.iter().max()) {
                let range = max_pos - min_pos;
                ensemble_ranges.push((ensemble.id, min_pos, max_pos));
            }
        }
    }

    if ensemble_ranges.is_empty() {
        return Err("No valid ensemble ranges found".into());
    }

    // Find maximum range (full-length ensembles)
    let max_range = ensemble_ranges.iter()
        .map(|(_, min_pos, max_pos)| max_pos - min_pos)
        .max()
        .unwrap_or(0);

    // Identify full-length ensembles (within 5% of max range)
    let full_length_threshold = max_range as f64 * 0.95;
    let mut full_length_ensembles: Vec<usize> = Vec::new();
    let mut non_full_length_ensembles: Vec<usize> = Vec::new();

    for (id, min_pos, max_pos) in &ensemble_ranges {
        let range = max_pos - min_pos;
        if range as f64 >= full_length_threshold {
            full_length_ensembles.push(*id);
        } else {
            non_full_length_ensembles.push(*id);
        }
    }

    if full_length_ensembles.is_empty() {
        eprintln!("[Duet Warning] No full-length ensembles found, skipping track plots");
        return Ok(());
    }

    if non_full_length_ensembles.is_empty() {
        eprintln!("[Duet Info] All ensembles are full-length, no overlap tracks to plot");
        return Ok(());
    }

    // Create output directory
    std::fs::create_dir_all(output_dir)?;

    // Get all positions covered by full-length ensembles
    let mut full_length_positions: HashSet<u32> = HashSet::new();
    for &full_id in &full_length_ensembles {
        if let Some(bases) = ensemble_data.get(&full_id) {
            for base in bases {
                full_length_positions.insert(base.position);
            }
        }
    }

    // For each non-full-length ensemble, plot overlap region
    for &non_full_id in &non_full_length_ensembles {
        if let Some(bases) = ensemble_data.get(&non_full_id) {
            // Find overlap positions
            let mut overlap_positions: Vec<u32> = bases.iter()
                .map(|b| b.position)
                .filter(|pos| full_length_positions.contains(pos))
                .collect();
            overlap_positions.sort();

            if overlap_positions.is_empty() {
                continue; // Skip if no overlap
            }

            // Get reactivity data for non-full-length ensemble overlap positions
            let mut non_full_reactivity_data: Vec<(u32, f64)> = Vec::new();
            for &pos in &overlap_positions {
                if let Some(base) = bases.iter().find(|b| b.position == pos) {
                    let reactivity = base.stop_reactivity.max(base.mutation_reactivity);
                    non_full_reactivity_data.push((pos, reactivity));
                }
            }

            // Get reactivity data for full-length ensemble(s) at overlap positions
            // Use average reactivity if multiple full-length ensembles
            let mut full_reactivity_data: Vec<(u32, f64)> = Vec::new();
            for &pos in &overlap_positions {
                let mut reactivities: Vec<f64> = Vec::new();
                for &full_id in &full_length_ensembles {
                    if let Some(full_bases) = ensemble_data.get(&full_id) {
                        if let Some(base) = full_bases.iter().find(|b| b.position == pos) {
                            let reactivity = base.stop_reactivity.max(base.mutation_reactivity);
                            reactivities.push(reactivity);
                        }
                    }
                }
                // Use average reactivity from all full-length ensembles
                let avg_reactivity = if reactivities.is_empty() {
                    f64::NAN
                } else {
                    reactivities.iter().sum::<f64>() / reactivities.len() as f64
                };
                full_reactivity_data.push((pos, avg_reactivity));
            }

            // Plot combined track: full-length on top, non-full-length on bottom
            let track_path = format!("{}/E{}_overlap_track.svg", output_dir, non_full_id);
            plot_combined_track(
                &full_reactivity_data,
                &non_full_reactivity_data,
                &track_path,
                non_full_id,
                &full_length_ensembles,
            )?;
        }
    }

    Ok(())
}

// Plot combined track: full-length ensemble on top, specific ensemble on bottom
fn plot_combined_track(
    full_reactivity_data: &[(u32, f64)],
    non_full_reactivity_data: &[(u32, f64)],
    output_path: &str,
    ensemble_id: usize,
    full_length_ids: &[usize],
) -> Result<(), Box<dyn Error>> {
    if full_reactivity_data.is_empty() || non_full_reactivity_data.is_empty() {
        return Err("No reactivity data to plot".into());
    }

    let positions: Vec<u32> = non_full_reactivity_data.iter().map(|(pos, _)| *pos).collect();
    let full_reactivities: Vec<f64> = full_reactivity_data.iter().map(|(_, react)| *react).collect();
    let non_full_reactivities: Vec<f64> = non_full_reactivity_data.iter().map(|(_, react)| *react).collect();

    let min_pos = positions.iter().min().copied().unwrap_or(0);
    let max_pos = positions.iter().max().copied().unwrap_or(0);
    let max_reactivity = full_reactivities.iter()
        .chain(non_full_reactivities.iter())
        .filter(|r| r.is_finite() && **r >= 0.0)
        .fold(0.0f64, |a, &b| a.max(b));

    // Calculate dimensions
    let title_height = 60;
    let margin = 50;
    let plot_width = 1200;
    let track_height = 150;
    let track_spacing = 30;
    let plot_height = track_height * 2 + track_spacing;
    let total_width = plot_width + margin * 2;
    let total_height = plot_height + title_height + margin * 2;

    // Create SVG
    let root = SVGBackend::new(output_path, (total_width as u32, total_height as u32)).into_drawing_area();
    root.fill(&WHITE)?;

    // Draw title
    let full_ids_str: Vec<String> = full_length_ids.iter().map(|id| format!("E{}", id)).collect();
    root.draw(&Text::new(
        format!("Ensemble E{} Overlap Region (Full-length: {})", ensemble_id, full_ids_str.join(", ")),
        (total_width as i32 / 2, 30),
        ("sans-serif", 18).into_font().color(&BLACK),
    ))?;

    // Calculate bar width
    let pos_range = (max_pos - min_pos).max(1) as f64;
    let bar_width = (plot_width as f64 / pos_range).max(1.0).min(10.0);

    // Draw full-length ensemble track (top)
    let top_y_base = title_height + margin + track_height;
    let top_y_top = title_height + margin;
    
    // Draw label for full-length track
    root.draw(&Text::new(
        "Full-length ensemble(s)",
        (margin - 5, top_y_top + track_height / 2),
        ("sans-serif", 12).into_font().color(&BLACK),
    ))?;

    for (idx, (pos, reactivity)) in full_reactivity_data.iter().enumerate() {
        let x_start = margin as f64 + ((*pos - min_pos) as f64 / pos_range) * plot_width as f64;
        let x_end = x_start + bar_width;
        let bar_height = if max_reactivity > 0.0 && reactivity.is_finite() && *reactivity >= 0.0 {
            (reactivity / max_reactivity) * track_height as f64
        } else {
            0.0
        };
        let y_top = top_y_base as f64 - bar_height;

        let color = get_reactivity_color(*reactivity);
        root.draw(&Rectangle::new(
            [(x_start as i32, y_top as i32), (x_end as i32, top_y_base)],
            color.filled(),
        ))?;
    }

    // Draw non-full-length ensemble track (bottom)
    let bottom_y_base = title_height + margin + track_height * 2 + track_spacing;
    let bottom_y_top = title_height + margin + track_height + track_spacing;
    
    // Draw label for non-full-length track
    root.draw(&Text::new(
        format!("Ensemble E{}", ensemble_id),
        (margin - 5, bottom_y_top + track_height / 2),
        ("sans-serif", 12).into_font().color(&BLACK),
    ))?;

    for (idx, (pos, reactivity)) in non_full_reactivity_data.iter().enumerate() {
        let x_start = margin as f64 + ((*pos - min_pos) as f64 / pos_range) * plot_width as f64;
        let x_end = x_start + bar_width;
        let bar_height = if max_reactivity > 0.0 && reactivity.is_finite() && *reactivity >= 0.0 {
            (reactivity / max_reactivity) * track_height as f64
        } else {
            0.0
        };
        let y_top = bottom_y_base as f64 - bar_height;

        let color = get_reactivity_color(*reactivity);
        root.draw(&Rectangle::new(
            [(x_start as i32, y_top as i32), (x_end as i32, bottom_y_base)],
            color.filled(),
        ))?;
    }

    // Draw axes
    let x_axis_y = bottom_y_base;
    root.draw(&Rectangle::new(
        [(margin, x_axis_y), (margin + plot_width, x_axis_y + 2)],
        BLACK.filled(),
    ))?;

    let y_axis_x = margin;
    root.draw(&Rectangle::new(
        [(y_axis_x - 2, title_height + margin), (y_axis_x, bottom_y_base)],
        BLACK.filled(),
    ))?;

    // Draw position labels (sample every 10% of positions)
    let label_interval = (positions.len() / 10).max(1);
    for (idx, &pos) in positions.iter().enumerate() {
        if idx % label_interval == 0 || idx == positions.len() - 1 {
            let x = margin as f64 + ((pos - min_pos) as f64 / pos_range) * plot_width as f64;
            root.draw(&Text::new(
                pos.to_string(),
                (x as i32, x_axis_y + 15),
                ("sans-serif", 10).into_font().color(&BLACK),
            ))?;
        }
    }

    // Draw reactivity labels
    root.draw(&Text::new(
        "0.0",
        (y_axis_x - 30, x_axis_y),
        ("sans-serif", 10).into_font().color(&BLACK),
    ))?;
    root.draw(&Text::new(
        format!("{:.1}", max_reactivity),
        (y_axis_x - 40, title_height + margin + 10),
        ("sans-serif", 10).into_font().color(&BLACK),
    ))?;

    // Draw legend
    let legend_x = margin + plot_width - 200;
    let legend_y = title_height + 20;
    root.draw(&Text::new(
        "Reactivity:",
        (legend_x, legend_y),
        ("sans-serif", 12).into_font().color(&BLACK),
    ))?;
    
    let legend_items = vec![
        ("0-0.3", RGBColor(0, 0, 0)),
        ("0.3-0.5", RGBColor(0, 0, 255)),
        ("0.5-0.7", RGBColor(255, 165, 0)),
        ("0.7-1.0", RGBColor(139, 0, 0)),
        ("No data", RGBColor(128, 128, 128)),
    ];
    
    for (idx, (label, color)) in legend_items.iter().enumerate() {
        let y = legend_y + 20 + (idx * 20) as i32;
        root.draw(&Rectangle::new(
            [(legend_x, y - 8), (legend_x + 15, y + 2)],
            color.filled(),
        ))?;
        root.draw(&Text::new(
            *label,
            (legend_x + 20, y - 4),
            ("sans-serif", 10).into_font().color(&BLACK),
        ))?;
    }

    root.present()?;
    Ok(())
}

// Helper function to get color based on reactivity value
fn get_reactivity_color(reactivity: f64) -> RGBColor {
    if reactivity.is_nan() || reactivity < 0.0 {
        // No data: gray
        RGBColor(128, 128, 128)
    } else if reactivity < 0.3 {
        // 0-0.3: black
        RGBColor(0, 0, 0)
    } else if reactivity < 0.5 {
        // 0.3-0.5: blue
        RGBColor(0, 0, 255)
    } else if reactivity < 0.7 {
        // 0.5-0.7: orange
        RGBColor(255, 165, 0)
    } else {
        // 0.7-1.0: dark red
        RGBColor(139, 0, 0)
    }
}

// Plot a single reactivity track (barplot) - kept for backward compatibility
fn plot_single_track(
    reactivity_data: &[(u32, f64)],
    output_path: &str,
    ensemble_id: usize,
) -> Result<(), Box<dyn Error>> {
    if reactivity_data.is_empty() {
        return Err("No reactivity data to plot".into());
    }

    let positions: Vec<u32> = reactivity_data.iter().map(|(pos, _)| *pos).collect();
    let reactivities: Vec<f64> = reactivity_data.iter().map(|(_, react)| *react).collect();

    let min_pos = positions.iter().min().copied().unwrap_or(0);
    let max_pos = positions.iter().max().copied().unwrap_or(0);
    let max_reactivity = reactivities.iter().fold(0.0f64, |a, &b| a.max(b));

    // Calculate dimensions
    let title_height = 60;
    let margin = 50;
    let plot_width = 1200;
    let plot_height = 200;
    let total_width = plot_width + margin * 2;
    let total_height = plot_height + title_height + margin * 2;

    // Create SVG
    let root = SVGBackend::new(output_path, (total_width as u32, total_height as u32)).into_drawing_area();
    root.fill(&WHITE)?;

    // Draw title
    root.draw(&Text::new(
        format!("Ensemble E{} Overlap Region Reactivity", ensemble_id),
        (total_width as i32 / 2, 30),
        ("sans-serif", 20).into_font().color(&BLACK),
    ))?;

    // Calculate bar width
    let pos_range = (max_pos - min_pos).max(1) as f64;
    let bar_width = (plot_width as f64 / pos_range).max(1.0).min(10.0);

    // Draw reactivity bars
    for (pos, reactivity) in reactivity_data {
        let x_start = margin as f64 + ((*pos - min_pos) as f64 / pos_range) * plot_width as f64;
        let x_end = x_start + bar_width;
        let y_base = title_height + margin + plot_height;
        let bar_height = if max_reactivity > 0.0 {
            (reactivity / max_reactivity) * plot_height as f64
        } else {
            0.0
        };
        let y_top = y_base as f64 - bar_height;

        // Choose color based on reactivity value
        let color = if reactivity.is_nan() || *reactivity < 0.0 {
            // No data: gray
            RGBColor(128, 128, 128)
        } else if *reactivity < 0.3 {
            // 0-0.3: black
            RGBColor(0, 0, 0)
        } else if *reactivity < 0.5 {
            // 0.3-0.5: blue
            RGBColor(0, 0, 255)
        } else if *reactivity < 0.7 {
            // 0.5-0.7: orange
            RGBColor(255, 165, 0)
        } else {
            // 0.7-1.0: dark red
            RGBColor(139, 0, 0)
        };

        root.draw(&Rectangle::new(
            [(x_start as i32, y_top as i32), (x_end as i32, y_base)],
            color.filled(),
        ))?;
    }

    // Draw axes
    let x_axis_y = title_height + margin + plot_height;
    root.draw(&Rectangle::new(
        [(margin, x_axis_y), (margin + plot_width, x_axis_y + 2)],
        BLACK.filled(),
    ))?;

    let y_axis_x = margin;
    root.draw(&Rectangle::new(
        [(y_axis_x - 2, title_height + margin), (y_axis_x, title_height + margin + plot_height)],
        BLACK.filled(),
    ))?;

    // Draw position labels (sample every 10% of positions)
    let label_interval = (positions.len() / 10).max(1);
    for (idx, &pos) in positions.iter().enumerate() {
        if idx % label_interval == 0 || idx == positions.len() - 1 {
            let x = margin as f64 + ((pos - min_pos) as f64 / pos_range) * plot_width as f64;
            root.draw(&Text::new(
                pos.to_string(),
                (x as i32, x_axis_y + 15),
                ("sans-serif", 10).into_font().color(&BLACK),
            ))?;
        }
    }

    // Draw reactivity labels
    root.draw(&Text::new(
        "0.0",
        (y_axis_x - 30, x_axis_y),
        ("sans-serif", 10).into_font().color(&BLACK),
    ))?;
    root.draw(&Text::new(
        format!("{:.1}", max_reactivity),
        (y_axis_x - 40, title_height + margin + 10),
        ("sans-serif", 10).into_font().color(&BLACK),
    ))?;

    // Draw legend
    let legend_x = margin + plot_width - 200;
    let legend_y = title_height + 20;
    root.draw(&Text::new(
        "Reactivity:",
        (legend_x, legend_y),
        ("sans-serif", 12).into_font().color(&BLACK),
    ))?;
    
    let legend_items = vec![
        ("0-0.3", RGBColor(0, 0, 0)),
        ("0.3-0.5", RGBColor(0, 0, 255)),
        ("0.5-0.7", RGBColor(255, 165, 0)),
        ("0.7-1.0", RGBColor(139, 0, 0)),
        ("No data", RGBColor(128, 128, 128)),
    ];
    
    for (idx, (label, color)) in legend_items.iter().enumerate() {
        let y = legend_y + 20 + (idx * 20) as i32;
        root.draw(&Rectangle::new(
            [(legend_x, y - 8), (legend_x + 15, y + 2)],
            color.filled(),
        ))?;
        root.draw(&Text::new(
            *label,
            (legend_x + 20, y - 4),
            ("sans-serif", 10).into_font().color(&BLACK),
        ))?;
    }

    root.present()?;
    Ok(())
}

fn perform_umap(data: &[Vec<f64>]) -> Result<(Vec<f64>, Vec<f64>), Box<dyn Error>> {
    if data.is_empty() || data[0].is_empty() {
        return Err("Empty data matrix".into());
    }

    let n = data.len();
    let p = data[0].len();

    if n < 2 {
        return Err("Need at least 2 samples for UMAP".into());
    }

    // Simplified UMAP implementation
    // Step 1: Compute pairwise distances
    let mut distances: Vec<Vec<f64>> = vec![vec![0.0; n]; n];
    for i in 0..n {
        for j in (i + 1)..n {
            let dist: f64 = data[i].iter()
                .zip(data[j].iter())
                .map(|(a, b)| (a - b) * (a - b))
                .sum::<f64>()
                .sqrt();
            distances[i][j] = dist;
            distances[j][i] = dist;
        }
    }

    // Step 2: Find k-nearest neighbors (k = min(15, n-1))
    let k = (15.min(n - 1)).max(1);
    let mut knn_distances: Vec<Vec<f64>> = vec![Vec::new(); n];
    for i in 0..n {
        let mut dists: Vec<(usize, f64)> = distances[i].iter()
            .enumerate()
            .filter(|(idx, _)| *idx != i)
            .map(|(idx, &dist)| (idx, dist))
            .collect();
        dists.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());
        knn_distances[i] = dists.iter().take(k).map(|(_, d)| *d).collect();
    }

    // Step 3: Compute sigma for each point (using k-th nearest neighbor distance)
    let mut sigmas: Vec<f64> = knn_distances.iter().map(|dists| {
        if dists.is_empty() {
            1.0
        } else {
            dists.last().copied().unwrap_or(1.0).max(1e-10)
        }
    }).collect();

    // Step 4: Initialize low-dimensional embedding (using PCA as starting point)
    let (mut y1, mut y2) = perform_pca(data)?;
    
    // Step 5: Optimize embedding using simplified UMAP objective
    // Simplified: use stochastic gradient descent with simplified cost function
    let learning_rate = 0.1;
    let iterations = 200;
    let min_dist = 0.1;

    for iter in 0..iterations {
        let mut y1_new = y1.clone();
        let mut y2_new = y2.clone();

        for i in 0..n {
            let mut grad1 = 0.0;
            let mut grad2 = 0.0;

            for j in 0..n {
                if i == j {
                    continue;
                }

                let high_dist = distances[i][j];
                let low_dist = ((y1[i] - y1[j]).powi(2) + (y2[i] - y2[j]).powi(2)).sqrt().max(1e-10);

                // High-dimensional probability (simplified)
                let high_prob = (-high_dist / sigmas[i]).exp().min(1.0);

                // Low-dimensional probability (simplified)
                let low_prob = 1.0 / (1.0 + (low_dist - min_dist) / 0.1);

                // Gradient
                let coeff = 2.0 * (high_prob - low_prob) / low_dist;
                grad1 += coeff * (y1[i] - y1[j]);
                grad2 += coeff * (y2[i] - y2[j]);
            }

            y1_new[i] += learning_rate * grad1;
            y2_new[i] += learning_rate * grad2;
        }

        // Adaptive learning rate
        let current_lr = learning_rate * (1.0 - (iter as f64 / iterations as f64));
        for i in 0..n {
            y1[i] += current_lr * (y1_new[i] - y1[i]);
            y2[i] += current_lr * (y2_new[i] - y2[i]);
        }
    }

    Ok((y1, y2))
}

fn hierarchical_cluster_rows(matrix: &[Vec<f64>]) -> Vec<usize> {
    if matrix.is_empty() {
        return Vec::new();
    }
    
    let n = matrix.len();
    if n == 1 {
        return vec![0];
    }
    
    // Compute pairwise distances (Euclidean distance)
    let mut distances: Vec<Vec<f64>> = vec![vec![0.0; n]; n];
    for i in 0..n {
        for j in (i + 1)..n {
            let dist: f64 = matrix[i].iter()
                .zip(matrix[j].iter())
                .map(|(a, b)| (a - b) * (a - b))
                .sum::<f64>()
                .sqrt();
            distances[i][j] = dist;
            distances[j][i] = dist;
        }
    }
    
    // Hierarchical clustering using average linkage
    // Start with each row as its own cluster
    let mut clusters: Vec<Vec<usize>> = (0..n).map(|i| vec![i]).collect();
    let mut cluster_distances: Vec<Vec<f64>> = distances.clone();
    
    // Merge clusters until we have a single cluster
    while clusters.len() > 1 {
        // Find the two closest clusters
        let mut min_dist = f64::INFINITY;
        let mut merge_i = 0;
        let mut merge_j = 1;
        
        for i in 0..clusters.len() {
            for j in (i + 1)..clusters.len() {
                // Compute average linkage distance between clusters
                let mut total_dist = 0.0;
                let mut count = 0;
                for &idx1 in &clusters[i] {
                    for &idx2 in &clusters[j] {
                        total_dist += distances[idx1][idx2];
                        count += 1;
                    }
                }
                let avg_dist = if count > 0 { total_dist / count as f64 } else { f64::INFINITY };
                
                if avg_dist < min_dist {
                    min_dist = avg_dist;
                    merge_i = i;
                    merge_j = j;
                }
            }
        }
        
        // Merge clusters
        let mut merged = clusters[merge_i].clone();
        merged.extend_from_slice(&clusters[merge_j]);
        clusters.remove(merge_j);
        clusters[merge_i] = merged;
    }
    
    // Return the order from the final cluster
    clusters[0].clone()
}

fn perform_pca(data: &[Vec<f64>]) -> Result<(Vec<f64>, Vec<f64>), Box<dyn Error>> {
    if data.is_empty() || data[0].is_empty() {
        return Err("Empty data matrix".into());
    }

    let n = data.len();
    let p = data[0].len();

    // Center the data (subtract mean)
    let mut means = vec![0.0; p];
    for row in data.iter() {
        for (j, &val) in row.iter().enumerate() {
            means[j] += val;
        }
    }
    for mean in &mut means {
        *mean /= n as f64;
    }

    let mut centered: Vec<Vec<f64>> = data.iter().map(|row| {
        row.iter().enumerate().map(|(j, &val)| val - means[j]).collect()
    }).collect();

    // Compute covariance matrix (simplified: just compute first two principal components)
    // For efficiency, we'll use power iteration to find the first two principal components
    
    // Initialize random vectors
    let mut v1 = vec![1.0 / (p as f64).sqrt(); p];
    let mut v2 = vec![1.0 / (p as f64).sqrt(); p];
    
    // Power iteration for first PC
    for _ in 0..50 {
        let mut new_v1 = vec![0.0; p];
        for i in 0..n {
            let dot_product: f64 = centered[i].iter().zip(v1.iter()).map(|(a, b)| a * b).sum();
            for j in 0..p {
                new_v1[j] += dot_product * centered[i][j];
            }
        }
        let norm: f64 = new_v1.iter().map(|x| x * x).sum::<f64>().sqrt();
        if norm > 1e-10 {
            for x in &mut new_v1 {
                *x /= norm;
            }
        }
        v1 = new_v1;
    }
    
    // Remove first PC component from data for second PC
    let mut centered2 = centered.clone();
    for row in &mut centered2 {
        let dot_product: f64 = row.iter().zip(v1.iter()).map(|(a, b)| a * b).sum();
        for j in 0..p {
            row[j] -= dot_product * v1[j];
        }
    }
    
    // Power iteration for second PC
    for _ in 0..50 {
        let mut new_v2 = vec![0.0; p];
        for i in 0..n {
            let dot_product: f64 = centered2[i].iter().zip(v2.iter()).map(|(a, b)| a * b).sum();
            for j in 0..p {
                new_v2[j] += dot_product * centered2[i][j];
            }
        }
        let norm: f64 = new_v2.iter().map(|x| x * x).sum::<f64>().sqrt();
        if norm > 1e-10 {
            for x in &mut new_v2 {
                *x /= norm;
            }
        }
        v2 = new_v2;
    }
    
    // Project data onto principal components
    let pc1: Vec<f64> = centered.iter().map(|row| {
        row.iter().zip(v1.iter()).map(|(a, b)| a * b).sum()
    }).collect();
    
    let pc2: Vec<f64> = centered.iter().map(|row| {
        row.iter().zip(v2.iter()).map(|(a, b)| a * b).sum()
    }).collect();
    
    Ok((pc1, pc2))
}

fn filter_positions_by_contiguity(detail: &mut ClusterDetail) {
    if detail.position_stats.is_empty() {
        return;
    }
    let mut positions: Vec<u32> = detail.positions.iter().copied().collect();
    positions.sort_unstable();
    let mut keep: HashSet<u32> = HashSet::new();
    let mut i = 0;
    while i < positions.len() {
        let mut j = i;
        while j + 1 < positions.len() && positions[j + 1] - positions[j] <= CONTIGUITY_MAX_GAP {
            j += 1;
        }
        let group_len = j - i + 1;
        if group_len >= CONTIGUITY_MIN_LEN {
            for idx in i..=j {
                keep.insert(positions[idx]);
            }
        }
        i = j + 1;
    }

    if keep.is_empty() {
        detail.positions.clear();
        detail.position_stats.clear();
        detail.position_confidence.clear();
        return;
    }

    detail.positions.retain(|pos| keep.contains(pos));
    detail.position_stats.retain(|pos, _| keep.contains(pos));
    detail
        .position_confidence
        .retain(|pos, _| keep.contains(pos));
}

fn compute_position_confidence(
    chr: &str,
    strand: &str,
    position: u32,
    reads: usize,
    stop_reactivity: f64,
    mutation_reactivity: f64,
    positions_map: &HashMap<String, ReactivityEntry>,
) -> f64 {
    let read_score = 1.0 - (-((reads as f64) / READ_SUPPORT_SCALE)).exp();
    let mut reactivity_score = (stop_reactivity.max(mutation_reactivity) / REACTIVITY_SCALE)
        .min(1.0)
        .max(0.0);

    if reactivity_score < HIGH_REACTIVITY_THRESHOLD {
        let mut neighbor_high = false;
        for offset in -NEIGHBOR_RANGE..=NEIGHBOR_RANGE {
            if offset == 0 {
                continue;
            }
            let neighbor_pos = position as i32 + offset;
            if neighbor_pos <= 0 {
                continue;
            }
            let key = make_key(chr, strand, neighbor_pos as u32);
            if let Some(entry) = positions_map.get(&key) {
                let neighbor_reactivity = entry.stop_reactivity.max(entry.mutation_reactivity);
                if neighbor_reactivity >= HIGH_REACTIVITY_THRESHOLD {
                    neighbor_high = true;
                    break;
                }
            }
        }
        if neighbor_high {
            reactivity_score = (reactivity_score + NEIGHBOR_BOOST).min(1.0);
        }
    }

    let confidence = READ_WEIGHT * read_score + (1.0 - READ_WEIGHT) * reactivity_score;
    confidence.min(1.0)
}

fn score_cluster_detail(
    detail: &mut ClusterDetail,
    positions_map: &HashMap<String, ReactivityEntry>,
) -> f64 {
    filter_positions_by_contiguity(detail);
    if detail.position_stats.is_empty() {
        detail.position_confidence.clear();
        return 0.0;
    }

    let mut weighted_sum = 0.0;
    let mut total_weight = 0.0;
    detail.position_confidence.clear();

    for (&pos, contrib) in detail.position_stats.iter() {
        let key = make_key(&detail.chr, &detail.strand, pos);
        let entry = positions_map.get(&key);
        let stop_reactivity = entry.map(|r| r.stop_reactivity).unwrap_or(0.0);
        let mutation_reactivity = entry.map(|r| r.mutation_reactivity).unwrap_or(0.0);
        let reads = contrib.stop_reads + contrib.mutation_reads;
        let confidence = compute_position_confidence(
            &detail.chr,
            &detail.strand,
            pos,
            reads,
            stop_reactivity,
            mutation_reactivity,
            positions_map,
        );
        detail.position_confidence.insert(pos, confidence);
        let weight = (reads.max(1)) as f64;
        weighted_sum += confidence * weight;
        total_weight += weight;
    }

    if total_weight > 0.0 {
        weighted_sum / total_weight
    } else {
        0.0
    }
}

struct UnionFind {
    parent: Vec<usize>,
    rank: Vec<usize>,
}

impl UnionFind {
    fn new(size: usize) -> Self {
        Self {
            parent: (0..size).collect(),
            rank: vec![0; size],
        }
    }

    fn find(&mut self, x: usize) -> usize {
        if self.parent[x] != x {
            let root = self.find(self.parent[x]);
            self.parent[x] = root;
        }
        self.parent[x]
    }

    fn union(&mut self, a: usize, b: usize) {
        let root_a = self.find(a);
        let root_b = self.find(b);
        if root_a == root_b {
            return;
        }
        if self.rank[root_a] < self.rank[root_b] {
            self.parent[root_a] = root_b;
        } else if self.rank[root_a] > self.rank[root_b] {
            self.parent[root_b] = root_a;
        } else {
            self.parent[root_b] = root_a;
            self.rank[root_a] += 1;
        }
    }
}
