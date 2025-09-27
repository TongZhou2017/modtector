// Main.rs content copied from stone project

// Version information constants
const VERSION: &str = env!("CARGO_PKG_VERSION");

use rust_htslib::bam::Read;

use rust_htslib::bam::record::Seq;
use bio::io::fasta;
use std::error::Error;
use std::time::Instant;
use std::collections::HashMap;
use std::io::{Write, BufWriter};

use std::path::Path;
use clap::{Parser, Subcommand, Args};
use rayon::prelude::*;

mod compare;
mod plot;
mod norm;
mod reactivity;
mod evaluate;
mod progress;

/// Logger manager supporting dynamic progress display and detailed logging
pub struct Logger {
    writer: BufWriter<std::fs::File>,
    last_progress: String,
}

impl Logger {
    pub fn new(file: std::fs::File) -> Self {
        Self {
            writer: BufWriter::new(file),
            last_progress: String::new(),
        }
    }
    
    /// Record detailed log information
    pub fn log(&mut self, message: &str) -> std::io::Result<()> {
        let timestamp = chrono::Utc::now().format("%Y-%m-%d %H:%M:%S");
        writeln!(self.writer, "[{}] {}", timestamp, message)?;
        self.writer.flush()?;
        Ok(())
    }
    
    /// Display dynamic progress information (overwrite previous line)
    pub fn progress(&mut self, message: &str) -> std::io::Result<()> {
        // Clear previous line
        if !self.last_progress.is_empty() {
            print!("\r{}", " ".repeat(self.last_progress.len()));
        }
        
        // Display new progress
        print!("\r{}", message);
        std::io::stdout().flush()?;
        
        self.last_progress = message.to_string();
        Ok(())
    }
    
    /// Finish progress display
    pub fn finish_progress(&mut self) -> std::io::Result<()> {
        if !self.last_progress.is_empty() {
            println!(); // New line
            self.last_progress.clear();
        }
        Ok(())
    }
    
    /// Record log and display progress simultaneously
    pub fn log_and_progress(&mut self, message: &str) -> std::io::Result<()> {
        self.log(message)?;
        self.progress(message)?;
        Ok(())
    }
}

#[derive(Parser)]
#[command(author, version, about, long_about = None)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

    #[derive(Subcommand)]
enum Commands {
    /// Data processing: Generate pileup data from BAM files
    Count(CountArgs),
    /// Sample signal normalization/filtering
    Norm(norm::NormArgs),
    /// Calculate reactivity scores
    Reactivity(reactivity::ReactivityArgs),
    /// Compare modified and unmodified samples, output differential CSV
    Compare(CompareArgs),
    /// Evaluate reactivity accuracy
    Evaluate(EvaluateArgs),
    /// Plot signal distribution graphs
    Plot(PlotArgs),
}

#[derive(Args)]
struct CountArgs {
    /// BAM file path
    #[arg(short = 'b', long = "bam")]
    pub bam: String,
    /// Reference FASTA file path
    #[arg(short = 'f', long = "fasta")]
    pub fasta: String,
    /// Output CSV path
    #[arg(short = 'o', long = "output")]
    pub output: String,
    /// Number of parallel threads
    #[arg(short = 't', long = "threads")]
    pub threads: Option<usize>,
    /// Log file path (optional)
    #[arg(short = 'l', long = "log")]
    pub log: Option<String>,
    /// Window size for genome segmentation (bases). If not set, no windowing.
    #[arg(short = 'w', long = "window")]
    pub window: Option<usize>,
}

#[derive(Args)]
struct CompareArgs {
    /// Comparison mode (mod-vs-unmod, reactivity-groups, biological-replicates, reactivity-results)
    #[arg(short = 'm', long = "mode", default_value = "mod-vs-unmod")]
    pub mode: String,
    
    // Original mod and unmod comparison parameters
    /// Modified sample CSV (mod-vs-unmod mode)
    #[arg(short = 'M', long = "mod")]
    pub mod_csv: Option<String>,
    /// Unmodified sample CSV (mod-vs-unmod mode)
    #[arg(short = 'U', long = "unmod")]
    pub unmod_csv: Option<String>,
    
    // Multi-repeat comparison parameters
    /// Group 1 file list (biological-replicates mode, comma-separated)
    #[arg(long = "group1")]
    pub group1_files: Option<String>,
    /// Group 2 file list (biological-replicates mode, comma-separated)
    #[arg(long = "group2")]
    pub group2_files: Option<String>,
    
    // reactivity group comparison parameters
    /// Group 1 reactivity file list (reactivity-groups mode, comma-separated)
    #[arg(long = "reactivity-group1")]
    pub reactivity_group1: Option<String>,
    /// Group 2 reactivity file list (reactivity-groups mode, comma-separated)
    #[arg(long = "reactivity-group2")]
    pub reactivity_group2: Option<String>,
    
    // reactivity results comparison parameters
    /// First reactivity result file (reactivity-results mode)
    #[arg(long = "reactivity1")]
    pub reactivity1: Option<String>,
    /// Second reactivity result file (reactivity-results mode)
    #[arg(long = "reactivity2")]
    pub reactivity2: Option<String>,
    
    /// Output CSV path
    #[arg(short = 'o', long = "output")]
    pub output: String,
    
    // Filtering parameters
    /// Minimum depth
    #[arg(short = 'd', long = "min-depth", default_value_t = 10)]
    pub min_depth: usize,
    /// Minimum fold change
    #[arg(short = 'f', long = "min-fold", default_value_t = 2.0)]
    pub min_fold: f64,
    
    // Statistical testing parameters
    /// Statistical test type (t-test, mann-whitney, wilcoxon, chi-square, continuity, diffscan, deltashape)
    #[arg(short = 't', long = "test", default_value = "t-test")]
    pub test_type: String,
    /// Significance level
    #[arg(short = 'a', long = "alpha", default_value_t = 0.05)]
    pub alpha: f64,
    
    /// Log file path (optional)
    #[arg(short = 'l', long = "log")]
    pub log: Option<String>,
}

#[derive(Args)]
struct EvaluateArgs {
    // Input files
    /// reactivity signal file
    #[arg(short = 'r', long = "reactivity")]
    pub reactivity_file: String,
    /// Secondary structure file
    #[arg(short = 's', long = "structure")]
    pub structure_file: String,
    
    // Output configuration
    /// Output directory
    #[arg(short = 'o', long = "output")]
    pub output_dir: String,
    /// Log file path (default to result directory)
    #[arg(short = 'l', long = "log")]
    pub log: Option<String>,
    
    // Basic configuration
    /// Signal type (stop, mutation, or both)
    #[arg(short = 't', long = "signal-type", default_value = "stop")]
    pub signal_type: String,
    
    // Gene configuration
    /// Gene ID (for base matching)
    #[arg(short = 'g', long = "gene-id")]
    pub gene_id: String,
    /// Strand information (+ or -)
    #[arg(short = 'S', long = "strand", default_value = "+")]
    pub strand: String,
    
    // Evaluation configuration
    /// Use base matching (default true)
    #[arg(long = "base-matching", default_value_t = true)]
    pub use_base_matching: bool,
    /// Use auto-shift correction (default true)
    #[arg(long = "auto-shift", default_value_t = true)]
    pub auto_shift: bool,
    
    // Performance configuration
    /// Use optimized version (default false)
    #[arg(short = 'O', long = "optimized", default_value_t = false)]
    pub optimized: bool,
}

#[derive(Args)]
struct PlotArgs {
    // Input files
    /// Modified sample CSV (optional for SVG-only mode)
    #[arg(short = 'M', long = "mod")]
    pub mod_csv: Option<String>,
    /// Unmodified sample CSV (optional for SVG-only mode)
    #[arg(short = 'U', long = "unmod")]
    pub unmod_csv: Option<String>,
    
    // Output configuration
    /// Output directory
    #[arg(short = 'o', long = "output")]
    pub output: String,
    
    // Threshold configuration
    /// Coverage threshold (default 0.2)
    #[arg(short = 'c', long = "coverage")]
    pub coverage_threshold: Option<f64>,
    /// Depth threshold (default 50 reads)
    #[arg(short = 'd', long = "depth")]
    pub depth_threshold: Option<f64>,
    
    // Optional input files
    /// reactivity csv file (containing stop and mutation signals, optional)
    #[arg(short = 'r', long = "reactivity")]
    pub reactivity_file: Option<String>,
    /// Genome annotation GFF/GTF file (optional)
    #[arg(short = 'g', long = "gff")]
    pub gff: Option<String>,
    
    // SVG plotting configuration
    /// SVG template file for RNA structure visualization (optional)
    #[arg(long = "svg-template")]
    pub svg_template: Option<String>,
    /// Bases to include in SVG plot (e.g., AC for DMS-seq)
    #[arg(long = "svg-bases", default_value = "ACGT")]
    pub svg_bases: String,
    /// Signal type to plot (stop, mutation, or all for multiple signals)
    #[arg(long = "svg-signal", default_value = "all")]
    pub svg_signal: String,
    /// Strand to include in SVG plot (+, -, or both)
    #[arg(long = "svg-strand", default_value = "+")]
    pub svg_strand: String,
    /// Reference sequence file for alignment (optional)
    #[arg(long = "svg-ref")]
    pub svg_ref: Option<String>,
    /// Maximum shift for alignment search
    #[arg(long = "svg-max-shift", default_value_t = 5)]
    pub svg_max_shift: u32,
    
    // Performance configuration
    /// Number of parallel threads
    #[arg(short = 't', long = "threads")]
    pub threads: Option<usize>,
    
    // Log configuration
    /// Log file path (optional)
    #[arg(short = 'l', long = "log")]
    pub log: Option<String>,
}


#[derive(Debug)]
#[derive(Clone)]
struct Item {
    tid: String,
    strand: char,
    position: u32,
    refbase: char,
    rfc: usize,
    pipc: usize,
    depth: usize,
    ins: usize,
    del: usize,
    a: u32, 
    c: u32, 
    g: u32, 
    t: u32,
}

// Read FASTA file, return a hash map: TID -> sequence
fn load_reference_sequences(fasta_path: &str) -> Result<HashMap<String, Vec<char>>, Box<dyn Error>> {
    let mut reference_sequences = HashMap::new();
    let reader = fasta::Reader::from_file(fasta_path)?;

    for result in reader.records() {
        let record = result?;
        let seq_id = record.id().to_string();
        let seq: Vec<char> = record.seq().iter().map(|&x| x as char).collect();
        reference_sequences.insert(seq_id, seq);
    }

    Ok(reference_sequences)
}

fn findchar(qpos: usize, seq: Seq<'_>) -> u8 {
    let znum = qpos / 2;
    let ynum = qpos % 2;
    let bnum: u8 = match ynum {
        0 => seq.encoded[znum] / 16,
        _ => seq.encoded[znum] % 16,
    };
    bnum
}

fn chartonum(c: char) -> u8 {
    match c {
        'A' => return 1,
        'C' => return 2,
        'G' => return 4,
        'T' => return 8,
        _ => return 15,
    }
}

fn itemcreate(tid: String, position: u32, refbase: char) -> (Item, Item) {
    // Forward strand item instance
    let positive_strand_item = Item {
        tid: tid.clone(), 
        strand: '+',
        position,
        refbase,
        rfc: 0,
        pipc: 0,
        depth: 0,
        ins: 0,
        del: 0,
        a: 0,
        c: 0,
        g: 0,
        t: 0,
    };

    // Reverse strand item instance
    let negative_strand_item = Item {
        tid: tid.clone(),
        strand: '-',
        position,
        refbase,
        rfc: 0,
        pipc: 0,
        depth: 0,
        ins: 0,
        del: 0,
        a: 0,
        c: 0,
        g: 0,
        t: 0,
    };

    (positive_strand_item, negative_strand_item)
}

// Normalization function: normalize base_density/rt_stop
fn normalize_signal(
    base_density: &mut HashMap<String, Vec<f64>>,
    rt_stop: &mut HashMap<String, Vec<f64>>,
    head_skip: usize,
    tail_skip: usize,
    norm_percentile_start: f64, // 0.1 means top 10%
    norm_percentile_end: f64,   // 0.2 means top 20%
    use_median: bool,           // true: use median, false: use mean
) {
    for (tid, bd_vec) in base_density.iter_mut() {
        let len = bd_vec.len();
        if len <= head_skip + tail_skip + 10 { continue; }
        let start = head_skip;
        let end = len - tail_skip;
        let mut region: Vec<f64> = bd_vec[start..end].iter().cloned().collect();
        region.sort_by(|a, b| b.partial_cmp(a).unwrap()); // Descending order
        let n = region.len();
        let s = (n as f64 * norm_percentile_start).round() as usize;
        let e = (n as f64 * norm_percentile_end).round() as usize;
        if s >= e || e > n { continue; }
        let slice = &region[s..e];
        let scaling = if use_median {
            let mut v = slice.to_vec();
            v.sort_by(|a, b| a.partial_cmp(b).unwrap());
            let m = v.len();
            if m == 0 { 1.0 } else if m % 2 == 0 { (v[m/2-1] + v[m/2]) / 2.0 } else { v[m/2] }
        } else {
            let sum: f64 = slice.iter().sum();
            if slice.len() == 0 { 1.0 } else { sum / slice.len() as f64 }
        };
        let scaling = if scaling > 1e-6 { scaling } else { 1.0 };
        // Normalize base_density
        for v in bd_vec.iter_mut() {
            *v /= scaling;
        }
        // Normalize rt_stop
        if let Some(rt_vec) = rt_stop.get_mut(tid) {
            for v in rt_vec.iter_mut() {
                *v /= scaling;
            }
        }
    }
}

fn main() -> Result<(), Box<dyn Error>> {
    let cli = Cli::parse();
    match cli.command {
        Commands::Count(args) => {
            // Validate count command parameters
            validate_count_args(&args)?;
            
            // Set up log file
            let log_file = if let Some(log_path) = &args.log {
                std::fs::File::create(log_path)?
            } else {
                std::fs::File::create("count.log")?
            };
            let mut logger = Logger::new(log_file);
            
            let num_threads = args.threads.unwrap_or_else(|| {
                std::thread::available_parallelism()
                    .map(|n| n.get())
                    .unwrap_or(4)
            });
            
            // Record environment information and parameters
            logger.log("=== ModDetector Count Function Log ===")?;
            logger.log(&format!("Software Version: v{}", VERSION))?;
            logger.log(&format!("Runtime: {}", chrono::Utc::now().format("%Y-%m-%d %H:%M:%S")))?;
            logger.log(&format!("BAM File: {}", args.bam))?;
            logger.log(&format!("FASTA File: {}", args.fasta))?;
            logger.log(&format!("Output File: {}", args.output))?;
            logger.log(&format!("Threads: {}", num_threads))?;
            logger.log("Starting BAM file processing...")?;
            
            let result = process_sample_pileup_ultra_fast(&args.bam, &args.fasta, &args.output, num_threads, args.window);
            
            match &result {
                Ok(_) => logger.log("BAM file processing completed")?,
                Err(e) => logger.log(&format!("BAM file processing failed: {}", e))?,
            }
            
            result
        }
                    Commands::Norm(args) => {
                let log_file = if let Some(log_path) = &args.log {
                    std::fs::File::create(log_path)?
                } else {
                    std::fs::File::create("norm.log")?
                };
                let mut logger = Logger::new(log_file);
                norm::norm_csv(&args, &mut logger)
            }
                    Commands::Reactivity(args) => {
                let log_file = if let Some(log_path) = &args.log {
                    std::fs::File::create(log_path)?
                } else {
                    std::fs::File::create("reactivity.log")?
                };
                let mut logger = Logger::new(log_file);
                reactivity::calc_reactivity(&args, &mut logger)
            }
        Commands::Evaluate(args) => {
            // Validate evaluate command parameters
            validate_evaluate_args(&args)?;
            
            // Create output directory
            std::fs::create_dir_all(&args.output_dir)?;
            
            // Set up log file
            let log_file = if let Some(log_path) = &args.log {
                std::fs::File::create(log_path)?
            } else {
                std::fs::File::create(format!("{}/evaluation.log", args.output_dir))?
            };
            let mut logger = Logger::new(log_file);
            
            // Select evaluation strategy based on signal type
            if args.signal_type == "both" {
                // Evaluate both stop and mutation signals simultaneously
                evaluate::evaluate_both_signals_main(
                    &args.reactivity_file,
                    &args.structure_file,
                    &args.output_dir,
                    &args.gene_id,
                    &args.strand,
                    &mut logger,
                )
            } else if args.optimized {
                // Optimized version evaluation
                evaluate::evaluate_reactivity_accuracy_main_optimized(
                    &args.reactivity_file,
                    &args.structure_file,
                    &args.output_dir,
                    &args.signal_type,
                    &args.strand,
                    &mut logger,
                )
            } else if args.use_base_matching && args.auto_shift {
                // 使用碱基匹配和自动shift校正的完整评估
                evaluate::evaluate_reactivity_accuracy_with_base_matching_main(
                    &args.reactivity_file,
                    &args.structure_file,
                    &args.output_dir,
                    &args.signal_type,
                    &args.gene_id,
                    &args.strand,
                    &mut logger,
                )
            } else if args.auto_shift {
                // 仅使用自动shift校正
                evaluate::evaluate_reactivity_accuracy_with_auto_shift_main(
                    &args.reactivity_file,
                    &args.structure_file,
                    &args.output_dir,
                    &args.signal_type,
                    &args.gene_id,
                    &args.strand,
                    &mut logger,
                )
            } else {
                // 基础评估（不使用碱基匹配和自动shift校正）
                evaluate::evaluate_reactivity_accuracy_main(
                    &args.reactivity_file,
                    &args.structure_file,
                    &args.output_dir,
                    &args.signal_type,
                    &mut logger,
                )
            }
        }
        Commands::Compare(args) => {
            // 验证compare命令参数
            validate_compare_args(&args)?;
            
            let log_file = if let Some(log_path) = &args.log {
                std::fs::File::create(log_path)?
            } else {
                std::fs::File::create("compare.log")?
            };
            let mut logger = Logger::new(log_file);
            
            // 解析统计检验类型
            let test_type = match args.test_type.as_str() {
                "t-test" => compare::StatisticalTest::TTest,
                "mann-whitney" => compare::StatisticalTest::MannWhitneyU,
                "wilcoxon" => compare::StatisticalTest::WilcoxonSignedRank,
                "chi-square" => compare::StatisticalTest::ChiSquare,
                "continuity" => compare::StatisticalTest::ContinuityCorrection,
                "diffscan" => compare::StatisticalTest::DiffScan,
                "deltashape" => compare::StatisticalTest::DeltaSHAPE,
                _ => return Err("Invalid test type. Supported: t-test, mann-whitney, wilcoxon, chi-square, continuity, diffscan, deltashape".into()),
            };
            
            match args.mode.as_str() {
                "mod-vs-unmod" => {
                    let mod_csv = args.mod_csv.as_ref().ok_or("Mod CSV file is required for mod-vs-unmod mode")?;
                    let unmod_csv = args.unmod_csv.as_ref().ok_or("Unmod CSV file is required for mod-vs-unmod mode")?;
                    
                    compare::compare_samples_with_stats(
                        mod_csv,
                        unmod_csv,
                        &args.output,
                        args.min_depth,
                        args.min_fold,
                        Some(test_type),
                        args.alpha,
                        &mut logger,
                    )
                },
                "biological-replicates" => {
                    let group1_files = args.group1_files.as_ref().ok_or("Group1 files are required for biological-replicates mode")?;
                    let group2_files = args.group2_files.as_ref().ok_or("Group2 files are required for biological-replicates mode")?;
                    
                    let group1_list: Vec<String> = group1_files.split(',').map(|s| s.trim().to_string()).collect();
                    let group2_list: Vec<String> = group2_files.split(',').map(|s| s.trim().to_string()).collect();
                    
                    compare::compare_biological_replicates(
                        &group1_list,
                        &group2_list,
                        &args.output,
                        test_type,
                        args.alpha,
                        args.min_depth,
                        &mut logger,
                    )
                },
                "reactivity-groups" => {
                    let group1_files = args.reactivity_group1.as_ref().ok_or("Reactivity group1 files are required for reactivity-groups mode")?;
                    let group2_files = args.reactivity_group2.as_ref().ok_or("Reactivity group2 files are required for reactivity-groups mode")?;
                    
                    let group1_list: Vec<String> = group1_files.split(',').map(|s| s.trim().to_string()).collect();
                    let group2_list: Vec<String> = group2_files.split(',').map(|s| s.trim().to_string()).collect();
                    
                    compare::compare_reactivity_groups(
                        &group1_list,
                        &group2_list,
                        &args.output,
                        test_type,
                        args.alpha,
                        &mut logger,
                    )
                },
                "reactivity-results" => {
                    let reactivity1 = args.reactivity1.as_ref().ok_or("First reactivity file is required for reactivity-results mode")?;
                    let reactivity2 = args.reactivity2.as_ref().ok_or("Second reactivity file is required for reactivity-results mode")?;
                    
                    compare::compare_reactivity_results(
                        reactivity1,
                        reactivity2,
                        &args.output,
                        test_type,
                        args.alpha,
                        &mut logger,
                    )
                },
                _ => return Err("Invalid mode. Supported: mod-vs-unmod, biological-replicates, reactivity-groups, reactivity-results".into()),
            }
        }
        Commands::Plot(args) => {
            // 验证plot命令参数
            validate_plot_args(&args)?;
            
            let log_file = if let Some(log_path) = &args.log {
                std::fs::File::create(log_path)?
            } else {
                std::fs::File::create("plot.log")?
            };
            let mut logger = Logger::new(log_file);
            
            // 检查是否只想画SVG
            let svg_only_mode = args.svg_template.is_some() && args.reactivity_file.is_some();
            
            if !svg_only_mode {
                // 执行常规的plot功能
                let mod_csv = args.mod_csv.as_ref().unwrap(); // 已经在validate_plot_args中验证过
                let unmod_csv = args.unmod_csv.as_ref().unwrap(); // 已经在validate_plot_args中验证过
                
                plot::plot_signal_distributions(
                    mod_csv,
                    unmod_csv,
                    &args.output,
                    args.threads,
                    args.coverage_threshold,
                    args.depth_threshold,
                    args.reactivity_file.as_deref(),
                    args.gff.as_deref(),
                    &mut logger,
                )?;
            } else {
                logger.log("SVG-only mode: Skipping regular plot generation")?;
            }
            
            // 如果提供了SVG模板，执行SVG绘图
            if let Some(svg_template) = &args.svg_template {
                if let Some(reactivity_file) = &args.reactivity_file {
                    logger.log("Starting SVG reactivity plotting...")?;
                    
                    // 确保输出目录存在
                    std::fs::create_dir_all(&args.output)?;
                    
                    // 根据信号类型决定绘图方式
                    if args.svg_signal == "all" {
                        // 绘制所有信号类型
                        logger.log("Auto-detecting signal types and plotting all...")?;
                        plot::plot_multiple_signals_to_svg(
                            reactivity_file,
                            svg_template,
                            &args.output,
                            &args.svg_bases,
                            &args.svg_strand,
                            args.svg_ref.as_deref(),
                            args.svg_max_shift,
                            None, // 使用默认颜色范围
                        )?;
                        logger.log("All signal types plotted successfully")?;
                    } else {
                        // 绘制指定信号类型
                        logger.log(&format!("Plotting signal type: {}", args.svg_signal))?;
                        let svg_output = format!("{}/rna_structure_colored_{}.svg", args.output, args.svg_signal);
                        
                        if let Some(ref_file) = &args.svg_ref {
                            // 使用参考序列对齐
                            plot::plot_reactivity_to_svg_with_alignment(
                                reactivity_file,
                                svg_template,
                                &svg_output,
                                &args.svg_bases,
                                &args.svg_signal,
                                &args.svg_strand,
                                ref_file,
                                args.svg_max_shift,
                                None, // 使用默认颜色范围
                            )?;
                            logger.log(&format!("SVG plot with alignment created: {}", svg_output))?;
                        } else {
                            // 不使用对齐
                            plot::plot_reactivity_to_svg(
                                reactivity_file,
                                svg_template,
                                &svg_output,
                                &args.svg_bases,
                                &args.svg_signal,
                                &args.svg_strand,
                                None, // 使用默认颜色范围
                            )?;
                            logger.log(&format!("SVG plot created: {}", svg_output))?;
                        }
                    }
                } else {
                    logger.log("Warning: SVG template provided but no reactivity file specified. Skipping SVG plotting.")?;
                }
            }
            
            Ok(())
        }
    }
}


/// 验证count命令参数
fn validate_count_args(args: &CountArgs) -> Result<(), Box<dyn Error>> {
    // 验证BAM文件
    if args.bam.trim().is_empty() {
        return Err("错误: BAM文件路径不能为空".into());
    }
    if !Path::new(&args.bam).exists() {
        return Err(format!("错误: BAM文件不存在: {}", args.bam).into());
    }
    if !args.bam.ends_with(".bam") {
        return Err(format!("错误: BAM文件路径必须以.bam结尾: {}", args.bam).into());
    }
    
    // 验证FASTA文件
    if args.fasta.trim().is_empty() {
        return Err("错误: FASTA文件路径不能为空".into());
    }
    if !Path::new(&args.fasta).exists() {
        return Err(format!("错误: FASTA文件不存在: {}", args.fasta).into());
    }
    if !args.fasta.ends_with(".fa") && !args.fasta.ends_with(".fasta") {
        return Err(format!("错误: FASTA文件路径必须以.fa或.fasta结尾: {}", args.fasta).into());
    }
    
    // 验证输出路径
    if args.output.trim().is_empty() {
        return Err("错误: 输出文件路径不能为空".into());
    }
    if !args.output.ends_with(".csv") {
        return Err(format!("错误: 输出文件路径必须以.csv结尾: {}", args.output).into());
    }
    
    // 验证线程数
    if let Some(threads) = args.threads {
        if threads == 0 {
            return Err("错误: 线程数不能为0".into());
        }
        if threads > 64 {
            return Err(format!("错误: 线程数不能超过64 (当前: {})", threads).into());
        }
    }
    // 验证窗口大小（可选）
    if let Some(w) = args.window {
        if w == 0 {
            return Err("错误: 窗口大小不能为0".into());
        }
        if w > 10_000_000 {
            return Err(format!("错误: 窗口大小不能超过10000000 (当前: {})", w).into());
        }
    }
    
    Ok(())
}

/// 验证compare命令参数
fn validate_compare_args(args: &CompareArgs) -> Result<(), Box<dyn Error>> {
    // 验证比较模式
    match args.mode.as_str() {
        "mod-vs-unmod" => {
            // 验证修饰样本文件
            let mod_csv = args.mod_csv.as_ref().ok_or("错误: mod-vs-unmod模式需要指定修饰样本文件")?;
            if mod_csv.trim().is_empty() {
                return Err("错误: 修饰样本文件路径不能为空".into());
            }
            if !Path::new(mod_csv).exists() {
                return Err(format!("错误: 修饰样本文件不存在: {}", mod_csv).into());
            }
            if !mod_csv.ends_with(".csv") {
                return Err(format!("错误: 修饰样本文件路径必须以.csv结尾: {}", mod_csv).into());
            }
            
            // 验证未修饰样本文件
            let unmod_csv = args.unmod_csv.as_ref().ok_or("错误: mod-vs-unmod模式需要指定未修饰样本文件")?;
            if unmod_csv.trim().is_empty() {
                return Err("错误: 未修饰样本文件路径不能为空".into());
            }
            if !Path::new(unmod_csv).exists() {
                return Err(format!("错误: 未修饰样本文件不存在: {}", unmod_csv).into());
            }
            if !unmod_csv.ends_with(".csv") {
                return Err(format!("错误: 未修饰样本文件路径必须以.csv结尾: {}", unmod_csv).into());
            }
        },
        "biological-replicates" => {
            // 验证组1文件
            let group1_files = args.group1_files.as_ref().ok_or("错误: biological-replicates模式需要指定组1文件")?;
            if group1_files.trim().is_empty() {
                return Err("错误: 组1文件列表不能为空".into());
            }
            for file in group1_files.split(',') {
                let file = file.trim();
                if file.is_empty() {
                    continue;
                }
                if !Path::new(file).exists() {
                    return Err(format!("错误: 组1文件不存在: {}", file).into());
                }
                if !file.ends_with(".csv") {
                    return Err(format!("错误: 组1文件路径必须以.csv结尾: {}", file).into());
                }
            }
            
            // 验证组2文件
            let group2_files = args.group2_files.as_ref().ok_or("错误: biological-replicates模式需要指定组2文件")?;
            if group2_files.trim().is_empty() {
                return Err("错误: 组2文件列表不能为空".into());
            }
            for file in group2_files.split(',') {
                let file = file.trim();
                if file.is_empty() {
                    continue;
                }
                if !Path::new(file).exists() {
                    return Err(format!("错误: 组2文件不存在: {}", file).into());
                }
                if !file.ends_with(".csv") {
                    return Err(format!("错误: 组2文件路径必须以.csv结尾: {}", file).into());
                }
            }
        },
        "reactivity-groups" => {
            // 验证reactivity组1文件
            let group1_files = args.reactivity_group1.as_ref().ok_or("错误: reactivity-groups模式需要指定组1文件")?;
            if group1_files.trim().is_empty() {
                return Err("错误: reactivity组1文件列表不能为空".into());
            }
            for file in group1_files.split(',') {
                let file = file.trim();
                if file.is_empty() {
                    continue;
                }
                if !Path::new(file).exists() {
                    return Err(format!("错误: reactivity组1文件不存在: {}", file).into());
                }
                if !file.ends_with(".csv") {
                    return Err(format!("错误: reactivity组1文件路径必须以.csv结尾: {}", file).into());
                }
            }
            
            // 验证reactivity组2文件
            let group2_files = args.reactivity_group2.as_ref().ok_or("错误: reactivity-groups模式需要指定组2文件")?;
            if group2_files.trim().is_empty() {
                return Err("错误: reactivity组2文件列表不能为空".into());
            }
            for file in group2_files.split(',') {
                let file = file.trim();
                if file.is_empty() {
                    continue;
                }
                if !Path::new(file).exists() {
                    return Err(format!("错误: reactivity组2文件不存在: {}", file).into());
                }
                if !file.ends_with(".csv") {
                    return Err(format!("错误: reactivity组2文件路径必须以.csv结尾: {}", file).into());
                }
            }
        },
        "reactivity-results" => {
            // 验证第一个reactivity文件
            let reactivity1 = args.reactivity1.as_ref().ok_or("错误: reactivity-results模式需要指定第一个reactivity文件")?;
            if reactivity1.trim().is_empty() {
                return Err("错误: 第一个reactivity文件路径不能为空".into());
            }
            if !Path::new(reactivity1).exists() {
                return Err(format!("错误: 第一个reactivity文件不存在: {}", reactivity1).into());
            }
            if !reactivity1.ends_with(".csv") {
                return Err(format!("错误: 第一个reactivity文件路径必须以.csv结尾: {}", reactivity1).into());
            }
            
            // 验证第二个reactivity文件
            let reactivity2 = args.reactivity2.as_ref().ok_or("错误: reactivity-results模式需要指定第二个reactivity文件")?;
            if reactivity2.trim().is_empty() {
                return Err("错误: 第二个reactivity文件路径不能为空".into());
            }
            if !Path::new(reactivity2).exists() {
                return Err(format!("错误: 第二个reactivity文件不存在: {}", reactivity2).into());
            }
            if !reactivity2.ends_with(".csv") {
                return Err(format!("错误: 第二个reactivity文件路径必须以.csv结尾: {}", reactivity2).into());
            }
        },
        _ => return Err("错误: 无效的比较模式。支持的模式: mod-vs-unmod, biological-replicates, reactivity-groups, reactivity-results".into()),
    }
    
    // 验证输出路径
    if args.output.trim().is_empty() {
        return Err("错误: 输出文件路径不能为空".into());
    }
    if !args.output.ends_with(".csv") {
        return Err(format!("错误: 输出文件路径必须以.csv结尾: {}", args.output).into());
    }
    
    // 验证最小深度
    if args.min_depth == 0 {
        return Err("错误: 最小深度不能为0".into());
    }
    if args.min_depth > 10000 {
        return Err(format!("错误: 最小深度不能超过10000 (当前: {})", args.min_depth).into());
    }
    
    // 验证最小fold change
    if args.min_fold <= 0.0 {
        return Err("错误: 最小fold change必须大于0".into());
    }
    if args.min_fold > 1000.0 {
        return Err(format!("错误: 最小fold change不能超过1000 (当前: {})", args.min_fold).into());
    }
    
    // 验证统计检验类型
    match args.test_type.as_str() {
        "t-test" | "mann-whitney" | "wilcoxon" | "chi-square" | "continuity" | "diffscan" | "deltashape" => {},
        _ => return Err("错误: 无效的统计检验类型。支持的检验: t-test, mann-whitney, wilcoxon, chi-square, continuity, diffscan, deltashape".into()),
    }
    
    // 验证显著性水平
    if args.alpha <= 0.0 || args.alpha >= 1.0 {
        return Err("错误: 显著性水平必须在0和1之间".into());
    }
    
    Ok(())
}

/// 验证evaluate命令参数
fn validate_evaluate_args(args: &EvaluateArgs) -> Result<(), Box<dyn Error>> {
    // 验证reactivity文件
    if args.reactivity_file.trim().is_empty() {
        return Err("错误: reactivity文件路径不能为空".into());
    }
    if !Path::new(&args.reactivity_file).exists() {
        return Err(format!("错误: reactivity文件不存在: {}", args.reactivity_file).into());
    }
    if !args.reactivity_file.ends_with(".csv") {
        return Err(format!("错误: reactivity文件路径必须以.csv结尾: {}", args.reactivity_file).into());
    }
    
    // 验证二级结构文件
    if args.structure_file.trim().is_empty() {
        return Err("错误: 二级结构文件路径不能为空".into());
    }
    if !Path::new(&args.structure_file).exists() {
        return Err(format!("错误: 二级结构文件不存在: {}", args.structure_file).into());
    }
    if !args.structure_file.ends_with(".dp") && !args.structure_file.ends_with(".ct") {
        return Err(format!("错误: 二级结构文件路径必须以.dp或.ct结尾: {}", args.structure_file).into());
    }
    
    // 验证输出目录
    if args.output_dir.trim().is_empty() {
        return Err("错误: 输出目录路径不能为空".into());
    }
    
    // 验证信号类型
    if args.signal_type != "stop" && args.signal_type != "mutation" && args.signal_type != "both" {
        return Err(format!("错误: 信号类型必须是'stop'、'mutation'或'both'，当前: {}", args.signal_type).into());
    }
    
    // 验证基因ID
    if args.gene_id.trim().is_empty() {
        return Err("错误: 基因ID不能为空".into());
    }
    
    // 验证链信息
    if args.strand != "+" && args.strand != "-" {
        return Err(format!("错误: 链信息必须是'+'或'-'，当前: {}", args.strand).into());
    }
    
    // 验证参数组合
    if args.optimized && (args.use_base_matching || args.auto_shift) {
        return Err("错误: 优化模式不能与base-matching或auto-shift参数同时使用".into());
    }
    
    if !args.use_base_matching && args.auto_shift {
        return Err("错误: auto-shift参数需要与base-matching参数同时使用".into());
    }
    
    Ok(())
}

/// 验证plot命令参数
fn validate_plot_args(args: &PlotArgs) -> Result<(), Box<dyn Error>> {
    // 检查是否只想画SVG
    let svg_only_mode = args.svg_template.is_some() && args.reactivity_file.is_some();
    
    // 如果不是SVG-only模式，验证修饰样本文件
    if !svg_only_mode {
        let mod_csv = args.mod_csv.as_ref().ok_or("错误: 修饰样本文件路径不能为空（SVG-only模式除外）")?;
        if mod_csv.trim().is_empty() {
            return Err("错误: 修饰样本文件路径不能为空".into());
        }
        if !Path::new(mod_csv).exists() {
            return Err(format!("错误: 修饰样本文件不存在: {}", mod_csv).into());
        }
        if !mod_csv.ends_with(".csv") {
            return Err(format!("错误: 修饰样本文件路径必须以.csv结尾: {}", mod_csv).into());
        }
        
        // 验证未修饰样本文件
        let unmod_csv = args.unmod_csv.as_ref().ok_or("错误: 未修饰样本文件路径不能为空（SVG-only模式除外）")?;
        if unmod_csv.trim().is_empty() {
            return Err("错误: 未修饰样本文件路径不能为空".into());
        }
        if !Path::new(unmod_csv).exists() {
            return Err(format!("错误: 未修饰样本文件不存在: {}", unmod_csv).into());
        }
        if !unmod_csv.ends_with(".csv") {
            return Err(format!("错误: 未修饰样本文件路径必须以.csv结尾: {}", unmod_csv).into());
        }
    }
    
    // 验证输出目录
    if args.output.trim().is_empty() {
        return Err("错误: 输出目录路径不能为空".into());
    }
    
    // 验证coverage阈值
    if let Some(coverage) = args.coverage_threshold {
        if coverage < 0.0 || coverage > 1.0 {
            return Err(format!("错误: coverage阈值必须在0.0到1.0之间，当前: {}", coverage).into());
        }
    }
    
    // 验证depth阈值
    if let Some(depth) = args.depth_threshold {
        if depth <= 0.0 {
            return Err("错误: depth阈值必须大于0".into());
        }
        if depth > 10000.0 {
            return Err(format!("错误: depth阈值不能超过10000，当前: {}", depth).into());
        }
    }
    
    // 验证线程数
    if let Some(threads) = args.threads {
        if threads == 0 {
            return Err("错误: 线程数不能为0".into());
        }
        if threads > 64 {
            return Err(format!("错误: 线程数不能超过64 (当前: {})", threads).into());
        }
    }
    
    // 验证可选的reactivity文件
    if let Some(reactivity_file) = &args.reactivity_file {
        if !Path::new(reactivity_file).exists() {
            return Err(format!("错误: reactivity文件不存在: {}", reactivity_file).into());
        }
        if !reactivity_file.ends_with(".csv") {
            return Err(format!("错误: reactivity文件路径必须以.csv结尾: {}", reactivity_file).into());
        }
    }
    
    // 验证GFF文件
    if let Some(gff) = &args.gff {
        if !Path::new(gff).exists() {
            return Err(format!("错误: GFF文件不存在: {}", gff).into());
        }
        if !gff.ends_with(".gff") && !gff.ends_with(".gff3") && !gff.ends_with(".gtf") {
            return Err(format!("错误: GFF文件路径必须以.gff、.gff3或.gtf结尾: {}", gff).into());
        }
    }
    
    // 验证SVG模板文件
    if let Some(svg_template) = &args.svg_template {
        if !Path::new(svg_template).exists() {
            return Err(format!("错误: SVG模板文件不存在: {}", svg_template).into());
        }
        if !svg_template.ends_with(".svg") {
            return Err(format!("错误: SVG模板文件路径必须以.svg结尾: {}", svg_template).into());
        }
    }
    
    // 验证SVG参考序列文件
    if let Some(svg_ref) = &args.svg_ref {
        if !Path::new(svg_ref).exists() {
            return Err(format!("错误: SVG参考序列文件不存在: {}", svg_ref).into());
        }
    }
    
    // 验证SVG信号类型参数
    if !matches!(args.svg_signal.as_str(), "all" | "stop" | "mutation" | "score") {
        return Err(format!("错误: 无效的SVG信号类型: {}. 支持的类型: all, stop, mutation, score", args.svg_signal).into());
    }
    
    // 验证SVG链类型参数
    if !matches!(args.svg_strand.as_str(), "+" | "-" | "both") {
        return Err(format!("错误: 无效的SVG链类型: {}. 支持的类型: +, -, both", args.svg_strand).into());
    }
    
    // 验证SVG碱基过滤字符串
    if args.svg_bases.is_empty() {
        return Err("错误: SVG碱基过滤字符串不能为空".into());
    }
    for base in args.svg_bases.chars() {
        if !"ACGT".contains(base.to_ascii_uppercase()) {
            return Err(format!("错误: SVG碱基过滤字符串包含无效碱基: {}", base).into());
        }
    }
    
    // 验证SVG最大shift
    if args.svg_max_shift > 100 {
        return Err(format!("错误: SVG最大shift不能超过100，当前: {}", args.svg_max_shift).into());
    }
    
    Ok(())
}


/// 超高效的多线程pileup处理，使用优化的I/O策略
#[derive(Clone)]
struct WindowTask { region: String, start: u32, end: u32 }

fn process_sample_pileup_ultra_fast(bam_path: &str, fasta_path: &str, output_path: &str, num_threads: usize, window: Option<usize>) -> Result<(), Box<dyn Error>> {
    if !Path::new(bam_path).exists() {
        return Err(format!("BAM文件不存在: {}", bam_path).into());
    }
    if !Path::new(fasta_path).exists() {
        return Err(format!("FASTA文件不存在: {}", fasta_path).into());
    }
    
    if let Some(w) = window {
        println!("[loading data] Threads={}, Window={}bp", num_threads, w);
    } else {
        println!("[loading data] Threads={} (no window)", num_threads);
    }
    println!("    BAM={}", bam_path);
    println!("    FASTA={}", fasta_path);
    println!();
    
    let start = Instant::now();
    
    // 加载参考序列
    let reference_sequences = load_reference_sequences(fasta_path)?;
    
    // 设置rayon线程池
    rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build_global()
        .unwrap();
    
    // 构建任务：窗口化或区域化
    let (num_chunks, progress_total, chunk_iter): (usize, usize, Vec<Vec<WindowTask>>) = if let Some(wsize) = window {
        let mut windows: Vec<WindowTask> = Vec::new();
        let overlap_size = 3u32; // 3bp重叠窗口，确保边界位点不丢失
        
        for (ref_name, ref_seq) in &reference_sequences {
            let ref_len = ref_seq.len() as u32;
            let step = wsize as u32;
            if step >= ref_len {
                windows.push(WindowTask { region: ref_name.clone(), start: 0, end: ref_len });
            } else {
                let mut s: u32 = 0;
                while s < ref_len {
                    let e = (s + step).min(ref_len);
                    windows.push(WindowTask { region: ref_name.clone(), start: s, end: e });
                    
                    // 计算下一个窗口的起始位置，考虑重叠
                    if e >= ref_len {
                        break; // 已到达序列末尾
                    }
                    s = e.saturating_sub(overlap_size); // 3bp重叠，避免下溢
                }
            }
        }

        // 近似覆盖度加权与负载均衡：
        // 1) 根据哈希打散不同染色体/位置，减少重区域聚集
        // 2) 创建远多于线程数的小批任务，动态并行
        use std::collections::hash_map::DefaultHasher;
        use std::hash::{Hash, Hasher};

        let mut keyed: Vec<(u64, WindowTask)> = windows
            .into_iter()
            .map(|w| {
                let mut hasher = DefaultHasher::new();
                // 将region、start作为key进行哈希，达到跨染色体交错与打散
                w.region.hash(&mut hasher);
                w.start.hash(&mut hasher);
                let key = hasher.finish();
                (key, w)
            })
            .collect();

        // 按哈希key排序，产生稳定的打散顺序
        keyed.sort_by_key(|(k, _)| *k);
        let windows_shuffled: Vec<WindowTask> = keyed.into_iter().map(|(_, w)| w).collect();

        // 计算小批任务大小：目标是任务数为线程数的20-50倍
        let target_multiplicity: usize = 32; // 经验值，可后续暴露为参数
        let mut batch_size = (windows_shuffled.len() / (num_threads.saturating_mul(target_multiplicity))).max(1);
        // 避免批过小导致调度/合并开销过大
        if batch_size > 4096 { batch_size = 4096; }

        let chunks: Vec<Vec<WindowTask>> = windows_shuffled
            .chunks(batch_size)
            .map(|c| c.to_vec())
            .collect();

        (chunks.len(), chunks.iter().map(|c| c.len()).sum(), chunks)
    } else {
        // 无窗口：每个region一个窗口
        let mut windows: Vec<WindowTask> = Vec::new();
        for (ref_name, ref_seq) in &reference_sequences {
            let ref_len = ref_seq.len() as u32;
            windows.push(WindowTask { region: ref_name.clone(), start: 0, end: ref_len });
        }
        let chunk_size = (windows.len() + num_threads - 1) / num_threads.max(1);
        let chunks: Vec<Vec<WindowTask>> = windows
            .chunks(chunk_size.max(1))
            .map(|chunk| chunk.to_vec())
            .collect();
        (chunks.len(), windows.len(), chunks)
    };
    
    // 使用原子计数器跟踪进度
    use std::sync::atomic::{AtomicUsize, Ordering};
    use std::sync::Arc;
    
    let completed = Arc::new(AtomicUsize::new(0));
    
    // 使用动态进度显示
    print!("\r[Running] Dividing {} tasks into {} chunks", 
           progress_total, num_chunks);
    std::io::stdout().flush()?;
    
    // 并行处理每个区域块
    let chunk_results: Vec<Result<(HashMap<String, HashMap<u32, Item>>, HashMap<String, HashMap<u32, Item>>, HashMap<String, Vec<u32>>), Box<dyn Error + Send + Sync>>> = chunk_iter
        .into_par_iter()
        .map(|window_chunk| {
            let result = process_region_chunk_ultra_fast(bam_path, &reference_sequences, window_chunk);
            
            // 更新进度
            let completed_count = completed.fetch_add(1, Ordering::Relaxed) + 1;
            let percentage = (completed_count as f64 * 100.0) / num_chunks as f64;
            print!("\r[Running] Processing {}/{} chunks ({:.1}%)", completed_count, num_chunks, percentage);
            std::io::stdout().flush().ok();
            
            result
        })
        .collect();
    
    // 覆盖显示完成状态
    println!("\r[Done] Processing {} chunks", num_chunks);
    let mut pos_items_map: HashMap<String, HashMap<u32, Item>> = HashMap::new();
    let mut neg_items_map: HashMap<String, HashMap<u32, Item>> = HashMap::new();
    let mut region_positions: HashMap<String, Vec<u32>> = HashMap::new();
    
    // 合并函数：同一位点的统计进行字段累加，而非覆盖
    // 重叠窗口会产生重复位点，此函数确保统计信息正确累加
    fn merge_item_into(target: &mut Item, src: &Item) {
        target.a += src.a;
        target.c += src.c;
        target.g += src.g;
        target.t += src.t;
        target.depth += src.depth;
        target.ins += src.ins;
        target.del += src.del;
        target.rfc += src.rfc;
        target.pipc += src.pipc;
        // strand、position、refbase 一致即可，若不一致保持target为主
    }

    for (chunk_idx, chunk_result) in chunk_results.into_iter().enumerate() {
        let (pos_items, neg_items, positions) = match chunk_result {
            Ok(result) => result,
            Err(e) => return Err(format!("并行处理错误 (块 {}): {}", chunk_idx, e).into()),
        };
        
        // 合并正链数据（累加）
        for (region, items) in pos_items {
            let target_map = pos_items_map.entry(region).or_insert_with(HashMap::new);
            for (pos, item) in items {
                if let Some(existing) = target_map.get_mut(&pos) {
                    merge_item_into(existing, &item);
                } else {
                    target_map.insert(pos, item);
                }
            }
        }
        
        // 合并负链数据（累加）
        for (region, items) in neg_items {
            let target_map = neg_items_map.entry(region).or_insert_with(HashMap::new);
            for (pos, item) in items {
                if let Some(existing) = target_map.get_mut(&pos) {
                    merge_item_into(existing, &item);
                } else {
                    target_map.insert(pos, item);
                }
            }
        }
        
        // 合并位置信息
        for (region, mut positions_vec) in positions {
            let target_vec = region_positions.entry(region).or_insert_with(Vec::new);
            target_vec.append(&mut positions_vec);
        }
    }
    
    // 对位置信息进行排序和去重（处理重叠窗口产生的重复位点）
    for (_, positions_vec) in region_positions.iter_mut() {
        positions_vec.sort();
        positions_vec.dedup();
    }
    
    // 使用缓冲写入优化输出性能
    let total_regions = reference_sequences.len();
    let mut write_progress = crate::progress::SimpleProgress::new(total_regions);
    
    print!("\r[Writing] Writing result file...");
    std::io::stdout().flush()?;
    let output_file = std::fs::File::create(output_path)?;
    let mut output = std::io::BufWriter::with_capacity(1024 * 1024, output_file); // 1MB缓冲区
    
    use std::io::Write;
    writeln!(output, "ChrID,pipe_truncation_Strand,pipe_truncation_ChrPos,rf_mutation_Base,rf_mutation_Count,pipe_truncation_count,depth,rf_mutation_ins,rf_mutation_del,base_A,base_C,base_G,base_T")?;
    
    // 批量写入所有区域的数据
    let mut write_buffer = String::with_capacity(1024 * 1024); // 1MB写入缓冲区
    
    for (i, (ref_name, ref_seq)) in reference_sequences.iter().enumerate() {
        write_progress.update(i)?;
        let empty_vec = Vec::new();
        let actual_positions = region_positions.get(ref_name).unwrap_or(&empty_vec);
        
        if actual_positions.is_empty() {
            // 如果该region没有任何位点，跳过输出（默认过滤空信号）
            continue;
        } else {
            // 输出该region所有实际有位点的数据
            for &pos in actual_positions {
                let ref_base = if pos < ref_seq.len() as u32 {
                    ref_seq[pos as usize]
                } else {
                    'N'
                };
                
                // 获取正链数据
                let pos_item = if let Some(items) = pos_items_map.get(ref_name) {
                    items.get(&pos).cloned().unwrap_or_else(|| {
                        let (pi, _) = itemcreate(ref_name.clone(), pos, ref_base);
                        pi
                    })
                } else {
                    let (pi, _) = itemcreate(ref_name.clone(), pos, ref_base);
                    pi
                };
                
                // 获取负链数据
                let neg_item = if let Some(items) = neg_items_map.get(ref_name) {
                    items.get(&pos).cloned().unwrap_or_else(|| {
                        let (_, ni) = itemcreate(ref_name.clone(), pos, ref_base);
                        ni
                    })
                } else {
                    let (_, ni) = itemcreate(ref_name.clone(), pos, ref_base);
                    ni
                };
                
                // 默认过滤模式：只输出有深度的链（与传统方法保持一致）
                if pos_item.depth > 0 {
                    write_buffer.push_str(&format!("{},{},{},{},{},{},{},{},{},{},{},{},{}\n",
                        ref_name, pos_item.strand, pos_item.position + 1, pos_item.refbase, 
                        pos_item.rfc, pos_item.pipc, pos_item.depth, pos_item.ins, pos_item.del, 
                        pos_item.a, pos_item.c, pos_item.g, pos_item.t));
                }
                
                if neg_item.depth > 0 {
                    write_buffer.push_str(&format!("{},{},{},{},{},{},{},{},{},{},{},{},{}\n",
                        ref_name, neg_item.strand, neg_item.position + 1, neg_item.refbase, 
                        neg_item.rfc, neg_item.pipc, neg_item.depth, neg_item.ins, neg_item.del, 
                        neg_item.a, neg_item.c, neg_item.g, neg_item.t));
                }
                
                // 当缓冲区达到一定大小时，批量写入文件
                if write_buffer.len() > 512 * 1024 { // 512KB
                    output.write_all(write_buffer.as_bytes())?;
                    write_buffer.clear();
                }
            }
        }
    }
    
    // 最后刷新一次，确保显示达到100%
    write_progress.finish()?;
    
    // 写入剩余的缓冲区内容
    if !write_buffer.is_empty() {
        output.write_all(write_buffer.as_bytes())?;
    }
    
    // 确保所有数据都写入磁盘
    output.flush()?;
    // 确保输出落盘后再完整打印换行，避免终端悬停
    {
        use std::io::Write as _;
        std::io::stdout().flush()?;
    }
    
    let duration = start.elapsed();
    println!("\r[Output] Pileup: {}", output_path);
    println!("{}", crate::progress::format_time_used(duration));
    {
        use std::io::Write as _;
        std::io::stdout().flush()?;
    }
    Ok(())
}

/// 超高效的区域块处理，使用优化的I/O策略
fn process_region_chunk_ultra_fast(
    bam_path: &str,
    reference_sequences: &HashMap<String, Vec<char>>,
    windows: Vec<WindowTask>
) -> Result<(HashMap<String, HashMap<u32, Item>>, HashMap<String, HashMap<u32, Item>>, HashMap<String, Vec<u32>>), Box<dyn Error + Send + Sync>> {
    
    use rust_htslib::bam::IndexedReader;
    
    // 使用优化的BAM读取器
    let mut bam = IndexedReader::from_path(bam_path)?;
    
    // 设置BAM读取器参数以优化性能
    bam.set_threads(1); // 避免内部多线程冲突
    
    let header = bam.header().to_owned();
    
    // 预分配内存以减少动态分配
    let mut pos_items_map: HashMap<String, HashMap<u32, Item>> = HashMap::new();
    let mut neg_items_map: HashMap<String, HashMap<u32, Item>> = HashMap::new();
    let mut region_positions: HashMap<String, Vec<u32>> = HashMap::new();
    
    // 对每个窗口进行查询处理
    for w in &windows {
        if !reference_sequences.contains_key(&w.region) { continue; }
        let ref_seq = &reference_sequences[&w.region];
        let ref_len = ref_seq.len() as u32;

        // 初始化该region容器（如果不存在）
        if !pos_items_map.contains_key(&w.region) {
            pos_items_map.insert(w.region.clone(), HashMap::new());
        }
        if !neg_items_map.contains_key(&w.region) {
            neg_items_map.insert(w.region.clone(), HashMap::new());
        }
        if !region_positions.contains_key(&w.region) {
            region_positions.insert(w.region.clone(), Vec::new());
        }

        // 为每个窗口创建独立的BAM读取器
        let mut window_bam = IndexedReader::from_path(bam_path)?;
        window_bam.set_threads(1);

        // 使用窗口查询，只处理当前窗口范围的数据
        let start1 = (w.start + 1) as u32; // 1-based inclusive
        let end1 = w.end as u32;           // 1-based inclusive end
        let region_query = format!("{}:{}-{}", w.region, start1, end1);

        // 查询该窗口的所有reads
        if let Ok(()) = window_bam.fetch(&region_query) {
            // 获取header
            let window_header = window_bam.header().to_owned();
            
            // 创建该区域的pileup，使用优化的参数
            let mut pileups = window_bam.pileup();
            pileups.set_max_depth(100_000);
            
            // 预分配临时数据结构
            let mut temp_positions = Vec::with_capacity(1000);
            let mut temp_pos_items = HashMap::new();
            let mut temp_neg_items = HashMap::new();
            
            // 处理该区域的pileup数据
            for pileup in pileups {
                let pileup = pileup?;
                let ref_name = String::from_utf8_lossy(window_header.tid2name(pileup.tid())).to_string();
                
                // 确保只处理当前区域
                if ref_name != w.region {
                    continue;
                }
                
                let pos = pileup.pos();
                
                // 检查位置是否在有效范围内
                if pos >= ref_len {
                    continue;
                }
                // 过滤到当前窗口范围（使用闭区间 [start, end]，允许stop跨窗写入）
                if pos < w.start || pos > w.end { continue; }
                
                let ref_base = ref_seq[pos as usize];
                
                // 记录当前region的位点（去重）
                if !temp_positions.contains(&pos) {
                    temp_positions.push(pos);
                }
                
                // 确保当前位点的Item存在
                if !temp_pos_items.contains_key(&pos) {
                    let (pi, ni) = itemcreate(w.region.clone(), pos, ref_base);
                    temp_pos_items.insert(pos, pi);
                    temp_neg_items.insert(pos, ni);
                }
                
                // 处理每个alignment
                for alignment in pileup.alignments() {
                    let record = alignment.record();
                    let id = alignment.indel();
                    
                    if let Some(qpos) = alignment.qpos() {
                        let base = findchar(qpos, record.seq());
                        
                        // 更新计数的闭包
                        let update_counts = |pi_or_ni: &mut Item, base| {
                            match base {
                                1 => pi_or_ni.a += 1,
                                2 => pi_or_ni.c += 1,
                                4 => pi_or_ni.g += 1,
                                8 => pi_or_ni.t += 1,
                                _ => {},
                            }
                            match id {
                                rust_htslib::bam::pileup::Indel::Ins(_) => {
                                    pi_or_ni.ins += 1;
                                }
                                rust_htslib::bam::pileup::Indel::Del(_) => {
                                    pi_or_ni.del += 1;
                                    pi_or_ni.depth += 1;
                                }
                                _ => {}
                            }
                            pi_or_ni.depth += 1;
                            let numb = chartonum(ref_base);
                            if numb != base {
                                pi_or_ni.rfc += 1;
                            }
                        };
                        
                        // stop信号累加
                        if alignment.is_head() {
                            let stop_signal_pos = if pos > 0 { pos - 1 } else { 0 };
                            
                            // 确保stop信号位点的Item存在
                            if !temp_pos_items.contains_key(&stop_signal_pos) {
                                let ref_base_stop = if stop_signal_pos < ref_len {
                                    ref_seq[stop_signal_pos as usize]
                                } else {
                                    'N'
                                };
                                let (pi, ni) = itemcreate(w.region.clone(), stop_signal_pos, ref_base_stop);
                                temp_pos_items.insert(stop_signal_pos, pi);
                                temp_neg_items.insert(stop_signal_pos, ni);
                            }
                            
                            // 记录stop信号位点
                            if !temp_positions.contains(&stop_signal_pos) {
                                temp_positions.push(stop_signal_pos);
                            }
                            
                            // 在正确的位点累加PIPC
                            if record.is_reverse() {
                                temp_neg_items.get_mut(&stop_signal_pos).unwrap().pipc += 1;
                            } else {
                                temp_pos_items.get_mut(&stop_signal_pos).unwrap().pipc += 1;
                            }
                        }
                        
                        // 更新当前位点的其他统计信息
                        if record.is_reverse() {
                            update_counts(temp_neg_items.get_mut(&pos).unwrap(), base);
                        } else {
                            update_counts(temp_pos_items.get_mut(&pos).unwrap(), base);
                        }
                    }
                }
            }
            
            // 批量更新主数据结构（合并到region级别），不做窗口级过滤，避免边界位点被忽略
            let pos_map = pos_items_map.get_mut(&w.region).unwrap();
            let neg_map = neg_items_map.get_mut(&w.region).unwrap();
            for (k, v) in temp_pos_items { pos_map.insert(k, v); }
            for (k, v) in temp_neg_items { neg_map.insert(k, v); }
            
            // 对位置进行排序并更新（无条件合并窗口内观测到的位点）
            let mut filtered_positions = temp_positions;
            filtered_positions.sort();
            let pos_vec = region_positions.get_mut(&w.region).unwrap();
            
            // 调试信息已移除
            
            // 合并并去重（简单方式，后续全局排序）
            for p in filtered_positions { if !pos_vec.contains(&p) { pos_vec.push(p); } }
        }
    }
    
    Ok((pos_items_map, neg_items_map, region_positions))
}