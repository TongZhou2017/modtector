// Main.rs content copied from stone project

// Version information constants
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

use rust_htslib::bam::Read;

use bio::io::fasta;
use rust_htslib::bam::record::Seq;
use std::collections::{HashMap, HashSet};
use std::error::Error;
use std::io::{BufWriter, Write};
use std::time::Instant;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;

use clap::{Args, Parser, Subcommand};
use rayon::prelude::*;
use std::path::Path;
use std::fs;
use glob::glob;
use ndarray::{Array1, Array2};

mod compare;
mod convert;
mod correct;
mod duet;
mod evaluate;
mod extract;
mod norm;
mod plot;
mod progress;
mod reactivity;

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

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand, Debug)]
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
    /// Duet analysis: infer alternative RNA conformations from normalized reactivity
    Duet(duet::DuetArgs),
    /// Extract gene regions from count results using GTF annotation
    Extract(ExtractArgs),
    /// Convert bamreadcount format to modtector pileup CSV format
    Convert(convert::ConvertArgs),
    /// Apply PCR bias correction to pileup CSV file
    Correct(correct::CorrectArgs),
}

#[derive(Args, Debug)]
struct CountArgs {
    /// BAM file path (or glob pattern for batch/single-cell mode)
    /// Batch mode: use glob pattern with wildcards (*, ?)
    /// Single-cell mode: use glob pattern, unified processing with cell labels
    /// Example: /path/to/bam/*sort.bam
    #[arg(short = 'b', long = "bam")]
    pub bam: String,
    /// Reference FASTA file path
    #[arg(short = 'f', long = "fasta")]
    pub fasta: String,
    /// Output CSV path (or directory for batch/single-cell mode)
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
    /// Strand filter: '+' for forward strand, '-' for reverse strand, '+/-' for both
    #[arg(short = 's', long = "strand", default_value = "+/-")]
    pub strand: String,
    /// Enable batch mode: process multiple BAM files sequentially
    /// Example: --batch /path/to/bam/*sort.bam
    #[arg(long = "batch", default_value_t = false)]
    pub batch: bool,
    /// Enable single-cell unified mode: unified processing with cell labels
    /// Reads are labeled by cell ID (extracted from filename), then processed together
    /// Example: --single-cell /path/to/bam/*sort.bam
    #[arg(long = "single-cell", default_value_t = false)]
    pub single_cell: bool,
    /// Minimum base quality score (Phred score) to count a mutation
    /// Only mutations with quality >= this threshold will be counted
    /// Set to 0 to disable base quality filtering (default: 0)
    /// Recommended value: 20 (same as RNAFramework rf-count)
    /// Note: Quality scores are decoded Phred scores (0-93), not ASCII
    #[arg(long = "min-base-qual", default_value_t = 0)]
    pub min_base_qual: u8,
    /// Enable PCR bias correction for depth calculation
    /// This corrects for PCR amplification bias where high depths cause mutation rate dilution
    /// Similar to shapemapper2's effective_depth calculation
    /// When enabled, effective_depth will be calculated and used instead of raw depth
    #[arg(long = "pcr-bias-correction", default_value_t = false)]
    pub pcr_bias_correction: bool,
    /// Weight for increasing depth correction (default: 1.0)
    /// Used when depth is too low and needs to be increased
    /// Matches Python script default
    #[arg(long = "weight-increase", default_value_t = 1.0)]
    pub weight_increase: f64,
    /// Weight for decreasing depth correction (default: 0.5)
    /// Used when depth is too high and needs to be decreased
    /// Matches Python script default
    #[arg(long = "weight-decrease", default_value_t = 0.5)]
    pub weight_decrease: f64,
}

#[derive(Args, Debug)]
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

#[derive(Args, Debug)]
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
    /// Reactive bases for AUC/ROC (e.g. AC for DMS; default ACGT = all)
    #[arg(long = "reactive-bases", default_value = "ACGT")]
    pub reactive_bases: String,

    // Performance configuration
    /// Use optimized version (default false)
    #[arg(short = 'O', long = "optimized", default_value_t = false)]
    pub optimized: bool,
}

#[derive(Args, Debug)]
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
    /// Circle fill type: true for filled circles, false for hollow rings
    #[arg(long = "svg-circle-filled", default_value_t = false)]
    pub svg_circle_filled: bool,
    /// Font color for text elements (e.g., "black", "red", "#000000")
    #[arg(long = "svg-font-color")]
    pub svg_font_color: Option<String>,
    /// Legend item width for horizontal layout (default: 30.0)
    #[arg(long = "svg-legend-width", default_value_t = 30.0)]
    pub svg_legend_width: f64,
    /// Legend item height for vertical layout (default: 15.0)
    #[arg(long = "svg-legend-height", default_value_t = 15.0)]
    pub svg_legend_height: f64,
    /// Generate interactive HTML visualization instead of static SVG
    #[arg(long = "interactive", default_value_t = false)]
    pub interactive: bool,

    // Performance configuration
    /// Number of parallel threads
    #[arg(short = 't', long = "threads")]
    pub threads: Option<usize>,

    // Log configuration
    /// Log file path (optional)
    #[arg(short = 'l', long = "log")]
    pub log: Option<String>,
}

#[derive(Args, Debug)]
struct ExtractArgs {
    /// Count result CSV file or bam-readcount output file
    #[arg(short = 'i', long = "input")]
    pub input: String,
    /// GTF annotation file
    #[arg(short = 'g', long = "gtf")]
    pub gtf: String,
    /// Target gene name or ID (e.g., "18S", "RN18S1", "rRNA")
    #[arg(short = 't', long = "target-gene")]
    pub target_gene: String,
    /// Output file prefix (each gene region will be saved as: {prefix}_{gene_name}_depth{avg_depth}.csv)
    #[arg(short = 'o', long = "output")]
    pub output: String,
    /// Input format: "csv" (default) or "bam-readcount"
    #[arg(long = "input-format", default_value = "csv")]
    pub input_format: String,
    /// Number of parallel threads for processing
    #[arg(long = "threads", default_value = "1")]
    pub threads: usize,
    /// Use relative position (1-based, relative to gene start) instead of absolute genomic position
    #[arg(long = "relative-position", default_value_t = false)]
    pub relative_position: bool,
    /// Log file path (optional)
    #[arg(short = 'l', long = "log")]
    pub log: Option<String>,
}

#[derive(Debug, Clone)]
struct Item {
    tid: String,
    strand: char,
    position: u32,
    refbase: char,
    rfc: usize,
    pipc: usize,
    depth: usize,
    effective_depth: usize,  // Quality-filtered depth (similar to shapemapper2)
    ins: usize,
    del: usize,
    a: u32,
    c: u32,
    g: u32,
    t: u32,
}

/// Chi-Square distribution PDF for df=2
/// y = a * chi2.pdf(x / scale, df=2)
fn chi_square_distribution_simple(x: f64, a: f64, scale: f64) -> f64 {
    let df = 2.0;
    let x_safe = if x == 0.0 { 1e-10 } else { x };
    let x_scaled = x_safe / scale;
    let x_scaled = x_scaled.max(1e-10);
    
    // chi2.pdf(x, df=2) = exp(-x/2) / 2
    let chi2_pdf = (-x_scaled / 2.0).exp() / 2.0;
    let result = a * chi2_pdf;
    result.max(0.0)
}

/// Chi-Square PDF for df=2 (inline version for optimization)
/// pdf(z) = 0.5 * exp(-z/2) where z = x/scale
#[inline]
fn chi2_pdf_df2(x_over_scale: f64) -> f64 {
    (-x_over_scale / 2.0).exp() * 0.5
}

/// Given scale, compute closed-form optimal a and SSE
/// Returns (sse, a) where a is projected to bounds [a_lb, a_ub]
fn sse_with_closed_form_a(
    binned: &[(f64, f64)],        // (depth, mutation_rate)
    scale: f64,
    a_lb: f64,
    a_ub: f64,
) -> (f64, f64) {
    let mut s1 = 0.0; // sum(y*g)
    let mut s2 = 0.0; // sum(g^2)

    for &(x, y) in binned {
        let z = x / scale;
        let g = chi2_pdf_df2(z);
        s1 += y * g;
        s2 += g * g;
    }

    // a*(scale) = (sum(y*g)) / (sum(g^2))
    let mut a = if s2 > 0.0 { s1 / s2 } else { a_lb };
    if !a.is_finite() {
        a = a_lb;
    }
    // project to bounds
    if a < a_lb { a = a_lb; }
    if a > a_ub { a = a_ub; }

    // SSE
    let mut sse = 0.0;
    for &(x, y) in binned {
        let g = chi2_pdf_df2(x / scale);
        let pred = a * g;
        let r = pred - y;
        sse += r * r;
    }
    (sse, a)
}

/// Brent's method for 1D minimization (simplified but stable)
/// Finds minimum of f in [a, c] with initial guess b (where a < b < c and f(b) <= f(a), f(c))
fn brent_minimize(
    f: &dyn Fn(f64) -> f64,
    mut a: f64, // left
    mut b: f64, // mid (near min)
    mut c: f64, // right
    tol: f64,
    max_iter: usize,
) -> (f64, f64) {
    // Ensure ordering: a < b < c
    if a > c {
        std::mem::swap(&mut a, &mut c);
    }
    if b < a { b = a; }
    if b > c { b = c; }

    let mut x = b;
    let mut w = b;
    let mut v = b;
    let mut fx = f(x);
    let mut fw = fx;
    let mut fv = fx;

    let mut d: f64 = 0.0;
    let mut e: f64 = 0.0;

    for _ in 0..max_iter {
        let xm = 0.5 * (a + c);
        let tol1 = tol * x.abs() + 1e-12;
        let tol2 = 2.0 * tol1;

        // Convergence check
        if (x - xm).abs() <= (tol2 - 0.5 * (c - a)).abs() {
            return (x, fx);
        }

        let mut p: f64 = 0.0;
        let mut q: f64 = 0.0;
        let mut r: f64 = 0.0;
        let mut use_parabolic = false;

        if e.abs() > tol1 {
            // Parabolic interpolation
            r = (x - w) * (fx - fv);
            q = (x - v) * (fx - fw);
            p = (x - v) * q - (x - w) * r;
            q = 2.0 * (q - r);
            if q > 0.0 { p = -p; }
            q = q.abs();

            let etemp = e;
            e = d;

            if p.abs() < 0.5 * q * etemp && p > q * (a - x) && p < q * (c - x) {
                // Accept parabolic step
                d = p / q;
                let u = x + d;
                if (u - a) < tol2 || (c - u) < tol2 {
                    d = if xm - x >= 0.0 { tol1 } else { -tol1 };
                }
                use_parabolic = true;
            }
        }

        if !use_parabolic {
            // Golden section step
            e = if x >= xm { a - x } else { c - x };
            d = 0.3819660112501051 * e; // 1 - 1/phi
        }

        let u = if d.abs() >= tol1 { x + d } else { x + if d > 0.0 { tol1 } else { -tol1 } };
        let fu = f(u);

        // Update bracket
        if fu <= fx {
            if u >= x { a = x; } else { c = x; }
            v = w; fv = fw;
            w = x; fw = fx;
            x = u; fx = fu;
        } else {
            if u < x { a = u; } else { c = u; }
            if fu <= fw || (w - x).abs() < 1e-18 {
                v = w; fv = fw;
                w = u; fw = fu;
            } else if fu <= fv || (v - x).abs() < 1e-18 || (v - w).abs() < 1e-18 {
                v = u; fv = fu;
            }
        }
    }

    (x, fx)
}

/// Minimize SSE(scale) in log-space using coarse grid + Brent refinement
/// Returns (best_scale, best_a, best_sse)
fn minimize_scale_logspace(
    binned: &[(f64, f64)],
    scale_lb: f64,
    scale_ub: f64,
    a_lb: f64,
    a_ub: f64,
) -> (f64, f64, f64) {
    let log_lb = scale_lb.ln();
    let log_ub = scale_ub.ln();

    // 1) Coarse grid search to find approximate minimum
    let grid_n = 120;
    let mut best_i = 0usize;
    let mut best_val = f64::INFINITY;
    let mut best_scale_grid = 0.0;

    let mut vals: Vec<f64> = Vec::with_capacity(grid_n);
    for i in 0..grid_n {
        let t = i as f64 / (grid_n as f64 - 1.0);
        let s = log_lb + t * (log_ub - log_lb);
        let scale = s.exp();
        let (sse, _a) = sse_with_closed_form_a(binned, scale, a_lb, a_ub);
        vals.push(sse);
        if sse < best_val {
            best_val = sse;
            best_i = i;
            best_scale_grid = scale;
        }
        // Print every 20th point for debugging
        if i % 20 == 0 || i == grid_n - 1 {
            eprintln!("[DEBUG] Grid i={}: scale={:.2}, sse={:.6e}", i, scale, sse);
        }
    }
    
    eprintln!("[DEBUG] Grid search: best_i={}, best_scale={:.2}, best_sse={:.6e}", best_i, best_scale_grid, best_val);

    // 2) Construct bracket: a < b < c where f(b) <= f(a), f(c)
    let i0 = best_i.saturating_sub(1);
    let i1 = best_i;
    let i2 = (best_i + 1).min(grid_n - 1);

    // Expand bracket if at boundary
    let left = if i0 == i1 { best_i.saturating_sub(2) } else { i0 };
    let right = if i2 == i1 { (best_i + 2).min(grid_n - 1) } else { i2 };

    let s_a = log_lb + (left as f64 / (grid_n as f64 - 1.0)) * (log_ub - log_lb);
    let s_b = log_lb + (best_i as f64 / (grid_n as f64 - 1.0)) * (log_ub - log_lb);
    let s_c = log_lb + (right as f64 / (grid_n as f64 - 1.0)) * (log_ub - log_lb);

    // Ensure ordering
    let s_a = s_a.min(s_c);
    let s_c = s_a.max(s_c);
    let s_b = s_b.max(s_a).min(s_c);

    let f = |s: f64| -> f64 {
        let scale = s.exp();
        let (sse, _a) = sse_with_closed_form_a(binned, scale, a_lb, a_ub);
        sse
    };

    // 3) Brent refinement
    eprintln!("[DEBUG] Brent bracket: s_a={:.6e} (scale={:.2}), s_b={:.6e} (scale={:.2}), s_c={:.6e} (scale={:.2})", 
              s_a, s_a.exp(), s_b, s_b.exp(), s_c, s_c.exp());
    
    // Print SSE at bracket points for comparison
    let (sse_a, _) = sse_with_closed_form_a(binned, s_a.exp(), a_lb, a_ub);
    let (sse_b, _) = sse_with_closed_form_a(binned, s_b.exp(), a_lb, a_ub);
    let (sse_c, _) = sse_with_closed_form_a(binned, s_c.exp(), a_lb, a_ub);
    eprintln!("[DEBUG] Brent bracket SSE: f(s_a)={:.6e}, f(s_b)={:.6e}, f(s_c)={:.6e}", sse_a, sse_b, sse_c);
    
    let (s_best, sse_best) = brent_minimize(&f, s_a, s_b, s_c, 1e-10, 200);
    let scale_best = s_best.exp();
    let (sse_final, a_final) = sse_with_closed_form_a(binned, scale_best, a_lb, a_ub);
    
    eprintln!("[DEBUG] Brent result: s_best={:.6e}, scale_best={:.2}, sse_best={:.6e}, sse_final={:.6e}, a_final={:.6e}", 
              s_best, scale_best, sse_best, sse_final, a_final);
    
    // Compare with grid search result
    eprintln!("[DEBUG] Comparison: grid_scale={:.2} (sse={:.6e}) vs brent_scale={:.2} (sse={:.6e})", 
              best_scale_grid, best_val, scale_best, sse_final);

    (scale_best, a_final, sse_final.min(sse_best))
}

/// Calculate minimum area line (median of fitted curve)
fn calculate_min_area_line(y_fit: &[f64]) -> f64 {
    let mut sorted = y_fit.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let n = sorted.len();
    if n == 0 {
        return 0.0;
    } else if n % 2 == 0 {
        (sorted[n / 2 - 1] + sorted[n / 2]) / 2.0
    } else {
        sorted[n / 2]
    }
}

/// Fit Chi-Square distribution to data (simplified version matching Python script)
/// Returns (a, scale, r2) parameters
fn fit_chi_square_to_data(
    depths: &[f64],
    mutation_rates: &[f64],
) -> Option<(f64, f64, f64)> {
    if depths.len() < 10 {
        return None;
    }
    
    // Bin data by depth using value-based bins (matching Python's pd.cut)
    let mut data: Vec<(f64, f64)> = depths.iter().zip(mutation_rates.iter())
        .map(|(d, mr)| (*d, *mr))
        .collect();
    data.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
    
    if data.is_empty() {
        return None;
    }
    
    let min_depth = data[0].0;
    let max_depth = data[data.len() - 1].0;
    // Always use 50 bins (matching Python script exactly)
    let n_bins = 50;
    
    if (max_depth - min_depth) < 1e-6 {
        return None;
    }
    
    // Create value-based bins (matching pd.cut exactly: right=True, i.e., (left, right])
    // pd.cut creates bins from min to max, with right-closed intervals
    let bin_width = (max_depth - min_depth) / n_bins as f64;
    let mut binned: Vec<(f64, f64)> = Vec::new();
    
    for i in 0..n_bins {
        // pd.cut uses right-closed intervals: (left, right]
        // First bin: (min_depth - epsilon, min_depth + bin_width]
        // Last bin: (min_depth + (n_bins-1)*bin_width, max_depth + epsilon]
        let bin_left = if i == 0 {
            min_depth - 1e-6 // Slightly below min to include first point
        } else {
            min_depth + (i as f64 * bin_width)
        };
        let bin_right = if i == n_bins - 1 {
            max_depth + 1e-6 // Slightly above max to include last point
        } else {
            min_depth + ((i + 1) as f64 * bin_width)
        };
        
        // Collect all points in this bin: left < d <= right (right-closed)
        let mut depths_in_bin = Vec::new();
        let mut rates_in_bin = Vec::new();
        
        for &(d, mr) in &data {
            if d > bin_left && d <= bin_right {
                depths_in_bin.push(d);
                rates_in_bin.push(mr);
            }
        }
        
        if !rates_in_bin.is_empty() {
            // Calculate median
            depths_in_bin.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
            rates_in_bin.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
            
            let median_depth = depths_in_bin[depths_in_bin.len() / 2];
            let median_rate = rates_in_bin[rates_in_bin.len() / 2];
            
            if median_rate > 0.0 {
                binned.push((median_depth, median_rate));
            }
        }
    }
    
    if binned.len() < 5 {
        return None;
    }
    
    let max_mutation = binned.iter().map(|(_, mr)| *mr).fold(0.0, f64::max);
    // Calculate median_depth matching Python's np.median (average of two middle values for even length)
    let median_depth = {
        let mut depths: Vec<f64> = binned.iter().map(|(d, _)| *d).collect();
        depths.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        let n = depths.len();
        if n == 0 {
            0.0
        } else if n % 2 == 0 {
            // For even length, average the two middle values (matching np.median)
            (depths[n / 2 - 1] + depths[n / 2]) / 2.0
        } else {
            depths[n / 2]
        }
    };
    
    // Try multiple initial parameters (matching Python script)
    let initial_params = vec![
        (max_mutation * 2.0, median_depth * 0.1),
        (max_mutation * 5.0, median_depth * 0.1),
        (max_mutation * 10.0, median_depth * 0.1),
        (max_mutation * 20.0, median_depth * 0.1),
        (max_mutation * 2.0, median_depth * 0.2),
        (max_mutation * 5.0, median_depth * 0.2),
    ];
    
    let mut best_r2 = f64::NEG_INFINITY;
    let mut best_params: Option<(f64, f64)> = None;
    
    // Use 1D optimization (eliminate linear parameter a) instead of 2D Trust Region
    // This matches scipy's behavior and avoids the a-scale degeneracy issue
    let bounds_a = (1e-6, 100.0);
    let bounds_scale = (1e-3, median_depth * 0.2);
    
    // Direct 1D global minimization in log-space
    // Debug: print bounds and binned data info
    eprintln!("[DEBUG] minimize_scale_logspace: bounds_scale=[{:.6e}, {:.2}], bounds_a=[{:.6e}, {:.1}], binned.len()={}", 
              bounds_scale.0, bounds_scale.1, bounds_a.0, bounds_a.1, binned.len());
    
    // Print first few binned data points for comparison with Python
    eprintln!("[DEBUG] First 5 binned points: {:?}", &binned[..binned.len().min(5)]);
    
    let (scale, a, sse) = minimize_scale_logspace(
        &binned,
        bounds_scale.0,
        bounds_scale.1,
        bounds_a.0,
        bounds_a.1,
    );
    
    eprintln!("[DEBUG] minimize_scale_logspace result: scale={:.2}, a={:.6e}, sse={:.6e}", scale, a, sse);
    
    // Calculate R² (matching existing logic)
    let mut y_obs: Vec<f64> = Vec::with_capacity(binned.len());
    let mut y_pred: Vec<f64> = Vec::with_capacity(binned.len());
    for &(d, mr_obs) in &binned {
        y_obs.push(mr_obs);
        y_pred.push(chi_square_distribution_simple(d, a, scale));
    }
    let mean_y = y_obs.iter().sum::<f64>() / (y_obs.len() as f64);
    let ss_res = sse;
    let ss_tot = y_obs.iter().map(|yy| (yy - mean_y) * (yy - mean_y)).sum::<f64>();
    let r2 = if ss_tot > 0.0 { 1.0 - ss_res / ss_tot } else { f64::NEG_INFINITY };
    
    if r2 > best_r2 {
        best_r2 = r2;
        best_params = Some((a, scale));
    }
    
    // Old 2D Trust Region code (replaced by 1D optimization above)
    // This is kept commented for reference but is no longer used
    /*
    for (_a_init, _scale_init) in initial_params {
        // Use Trust Region Reflective optimization (matching scipy's 'trf' method)
        // Ensure strict feasibility: adjust initial parameters if they're on bounds
        let mut a = a_init.max(1e-6).min(100.0);
        let mut scale = scale_init.max(1e-3).min(median_depth * 0.2);
        let bounds_a = (1e-6, 100.0);
        let bounds_scale = (1e-3, median_depth * 0.2);
        
        // Make strictly feasible (matching scipy's make_strictly_feasible)
        // If parameter is at or very close to bound, move it slightly inward
        if (a - bounds_a.0).abs() < 1e-8 {
            a = bounds_a.0 + 1e-6;
        } else if (a - bounds_a.1).abs() < 1e-8 {
            a = bounds_a.1 - 1e-6;
        }
        if (scale - bounds_scale.0).abs() < 1e-8 {
            scale = bounds_scale.0 + 1e-6;
        } else if (scale - bounds_scale.1).abs() < 1e-8 {
            scale = bounds_scale.1 - 1e-6;
        }
        
        // Trust region radius initialization (matching scipy's TRF)
        // scipy uses initial_trust_radius = 1.0, but adapts based on step quality
        // For bounded problems, scipy allows very large steps (step_norm up to 9+)
        // We need a large initial radius to allow jumping from scale=52.8 to scale=7.22
        // scipy's step_norm starts at ~9.16, so we need initial_radius >= 10
        let param_scale = (a.abs() + scale.abs()).sqrt();
        // Use much larger initial radius to match scipy's behavior
        // scipy allows step_norm up to ~9, so initial radius should be at least 10
        let mut trust_radius = param_scale.max(50.0).min(1000.0);  // Start with larger radius
        // Increase max_trust_radius significantly to allow very large steps
        let max_trust_radius = 1000.0;
        let min_trust_radius = 1e-8;
        
        let mut best_error = f64::INFINITY;
        let mut best_a = a;
        let mut best_scale = scale;
        
        // Convergence tolerances (matching scipy defaults: ftol=1e-8, xtol=1e-8, gtol=1e-8)
        let ftol = 1e-8;
        let xtol = 1e-8;
        let gtol = 1e-8;
        let mut prev_error = f64::INFINITY;
        let mut prev_a = a;
        let mut prev_scale = scale;
        
        // Maximum function evaluations (increased to allow better convergence)
        // Python script uses maxfev=50, but we increase it to allow more iterations
        // This helps when starting from scale=52.8 and needing to reach scale=7.22
        let max_nfev = 200;
        let mut nfev = 0;
        
        // Trust region parameters (matching scipy's TRF more closely)
        // scipy uses eta around 0.125 for step acceptance
        let eta1 = 0.125;  // Lower threshold for step acceptance (matching scipy default)
        let eta2 = 0.75;  // Upper threshold for step acceptance
        let gamma1 = 0.5;  // Trust radius reduction factor
        let gamma2 = 2.0;  // Trust radius increase factor
        
        for iteration in 0..max_nfev {
            // Calculate current error (residual sum of squares)
            let mut current_error = 0.0;
            for &(d, mr_obs) in &binned {
                let mr_pred = chi_square_distribution_simple(d, a, scale);
                let residual = mr_pred - mr_obs;
                current_error += residual * residual;
            }
            nfev += 1;
            
            if current_error < best_error {
                best_error = current_error;
                best_a = a;
                best_scale = scale;
            }
            
            prev_error = current_error;
            prev_a = a;
            prev_scale = scale;
            
            // Calculate Jacobian (gradients of residuals)
            let mut jacobian = Vec::new();
            let eps = 1e-8; // Smaller epsilon for better gradient accuracy
            
            for &(d, mr_obs) in &binned {
                let mr_pred = chi_square_distribution_simple(d, a, scale);
                let residual = mr_pred - mr_obs;
                
                // Numerical gradients
                let mr_a_eps = chi_square_distribution_simple(d, a + eps, scale);
                let mr_s_eps = chi_square_distribution_simple(d, a, scale + eps);
                let grad_a = (mr_a_eps - mr_pred) / eps;
                let grad_scale = (mr_s_eps - mr_pred) / eps;
                
                jacobian.push((grad_a, grad_scale, residual));
            }
            
            // Calculate gradient of cost function: J^T * r
            let mut grad_cost_a = 0.0;
            let mut grad_cost_scale = 0.0;
            for &(jac_a, jac_scale, residual) in &jacobian {
                grad_cost_a += jac_a * residual;
                grad_cost_scale += jac_scale * residual;
            }
            
            // Calculate approximate Hessian: J^T * J
            let mut hess_aa = 0.0;
            let mut hess_ss = 0.0;
            let mut hess_as = 0.0;
            for &(jac_a, jac_scale, _) in &jacobian {
                hess_aa += jac_a * jac_a;
                hess_ss += jac_scale * jac_scale;
                hess_as += jac_a * jac_scale;
            }
            
            // Check convergence (matching scipy's termination conditions)
            // Now we have grad_cost_a and grad_cost_scale, so we can check gtol
            if iteration > 0 {
                // ftol: relative change in cost function
                let error_change = (prev_error - current_error).abs() / (prev_error + 1e-10);
                // xtol: relative change in parameters
                let param_change_a = (a - prev_a).abs() / (prev_a.abs() + 1e-10);
                let param_change_scale = (scale - prev_scale).abs() / (prev_scale.abs() + 1e-10);
                let max_param_change = param_change_a.max(param_change_scale);
                
                // gtol: gradient norm
                let grad_norm = (grad_cost_a * grad_cost_a + grad_cost_scale * grad_cost_scale).sqrt();
                let grad_norm_scaled = grad_norm / (current_error.sqrt() + 1e-10);
                
                if error_change < ftol && max_param_change < xtol && grad_norm_scaled < gtol {
                    break; // Converged
                }
            }
            
            // Trust Region Reflective step (improved implementation matching scipy's TRF)
            // Solve trust region subproblem: minimize ||J*delta + r||^2 subject to ||delta|| <= trust_radius
            // Using Levenberg-Marquardt approach: (J^T*J + lambda*I) * delta = -J^T*r
            // where lambda is chosen such that ||delta|| = trust_radius
            
            // Initialize lambda based on Hessian eigenvalues (matching scipy's approach)
            // scipy uses a more sophisticated initialization based on the problem scale
            // For large steps (step_norm ~9), lambda needs to be very small
            let hess_trace = hess_aa + hess_ss;
            let hess_det = hess_aa * hess_ss - hess_as * hess_as;
            // Use much smaller initial lambda to allow very large steps
            // scipy's step_norm starts at ~9, which requires very small lambda
            let lambda_init = if hess_trace > 0.0 {
                (hess_trace / 100.0).max(1e-10)  // Much smaller initial lambda
            } else {
                1e-5  // Much smaller initial lambda
            };
            let mut lambda = lambda_init;
            let mut step_accepted = false;
            let max_trust_iter = 50;  // Further increased iterations for better convergence
            
            for trust_iter in 0..max_trust_iter {
                let det = (hess_aa + lambda) * (hess_ss + lambda) - hess_as * hess_as;
                
                if det.abs() > 1e-12 {
                    // Solve for step direction
                    let delta_a = -((hess_ss + lambda) * grad_cost_a - hess_as * grad_cost_scale) / det;
                    let delta_scale = -((hess_aa + lambda) * grad_cost_scale - hess_as * grad_cost_a) / det;
                    
                    // Check if step is within trust region
                    let step_norm = (delta_a * delta_a + delta_scale * delta_scale).sqrt();
                    
                    // Apply reflective step for bounds (matching scipy's trf method)
                    // This is a simplified version - scipy uses more sophisticated reflection
                    let mut new_a = a + delta_a;
                    let mut new_scale = scale + delta_scale;
                    let mut reflection_count = 0;
                    const MAX_REFLECTIONS: usize = 10;
                    
                    // Apply reflections until parameters are within bounds
                    while reflection_count < MAX_REFLECTIONS {
                        let mut reflected = false;
                        
                        if new_a < bounds_a.0 {
                            new_a = 2.0 * bounds_a.0 - new_a;
                            reflected = true;
                        } else if new_a > bounds_a.1 {
                            new_a = 2.0 * bounds_a.1 - new_a;
                            reflected = true;
                        }
                        
                        if new_scale < bounds_scale.0 {
                            new_scale = 2.0 * bounds_scale.0 - new_scale;
                            reflected = true;
                        } else if new_scale > bounds_scale.1 {
                            new_scale = 2.0 * bounds_scale.1 - new_scale;
                            reflected = true;
                        }
                        
                        if !reflected {
                            break;  // Parameters are within bounds
                        }
                        reflection_count += 1;
                    }
                    
                    // Ensure strict feasibility after reflection
                    new_a = new_a.max(bounds_a.0 + 1e-10).min(bounds_a.1 - 1e-10);
                    new_scale = new_scale.max(bounds_scale.0 + 1e-10).min(bounds_scale.1 - 1e-10);
                    
                    // Check if step (after reflection) is within trust region
                    let actual_step_a = new_a - a;
                    let actual_step_scale = new_scale - scale;
                    let actual_step_norm = (actual_step_a * actual_step_a + actual_step_scale * actual_step_scale).sqrt();
                    
                    if actual_step_norm <= trust_radius * 1.1 || trust_iter == max_trust_iter - 1 {
                        // Calculate new error
                        let mut new_error = 0.0;
                        for &(d, mr_obs) in &binned {
                            let mr_pred = chi_square_distribution_simple(d, new_a, new_scale);
                            let residual = mr_pred - mr_obs;
                            new_error += residual * residual;
                        }
                        nfev += 1;
                        
                        // Calculate actual vs predicted reduction
                        // Predicted reduction using quadratic model
                        let predicted_reduction = 0.5 * (
                            hess_aa * delta_a * delta_a +
                            2.0 * hess_as * delta_a * delta_scale +
                            hess_ss * delta_scale * delta_scale
                        ) + (grad_cost_a * delta_a + grad_cost_scale * delta_scale);
                        let actual_reduction = current_error - new_error;
                        
                        // Calculate rho (actual/predicted reduction ratio)
                        let rho = if predicted_reduction.abs() > 1e-12 {
                            actual_reduction / predicted_reduction
                        } else {
                            if actual_reduction > 0.0 {
                                1.0  // Good step even if prediction is poor
                            } else {
                                -1.0  // Bad step
                            }
                        };
                        
                        // Update trust region based on rho (matching scipy's TRF)
                        // scipy uses a more nuanced strategy for updating trust radius
                        if rho > eta2 {
                            // Very good step (rho > 0.75)
                            // Increase trust radius aggressively
                            trust_radius = (trust_radius * gamma2).min(max_trust_radius);
                            a = new_a;
                            scale = new_scale;
                            step_accepted = true;
                            // Reset lambda to allow larger steps in next iteration
                            lambda = lambda_init * 0.1;  // Reset to smaller value
                            break;
                        } else if rho > eta1 {
                            // Acceptable step (0.125 < rho <= 0.75)
                            // Keep or slightly increase trust radius
                            if rho > 0.5 {
                                // Good step, increase trust radius
                                trust_radius = (trust_radius * 1.5).min(max_trust_radius);
                            } else {
                                // Acceptable but not great, keep trust radius
                                trust_radius = trust_radius.min(max_trust_radius);
                            }
                            a = new_a;
                            scale = new_scale;
                            step_accepted = true;
                            // Reduce lambda slightly to allow larger steps
                            lambda = lambda.max(1e-8) * 0.8;
                            break;
                        } else {
                            // Poor step (rho <= 0.125)
                            // Reduce trust radius
                            trust_radius = (trust_radius * gamma1).max(min_trust_radius);
                            if trust_radius < min_trust_radius {
                                // Trust region too small, accept step anyway if it improves
                                if actual_reduction > 0.0 {
                                    a = new_a;
                                    scale = new_scale;
                                    step_accepted = true;
                                }
                                break;
                            }
                            // Increase lambda to reduce step size
                            lambda = lambda.max(1e-8) * 2.0;
                        }
                    } else {
                        // Step too large, increase lambda to reduce step size
                        // But also try to increase trust radius if possible
                        if trust_radius < max_trust_radius {
                            trust_radius = (trust_radius * 1.1).min(max_trust_radius);
                        }
                        lambda = lambda.max(1e-8) * 2.0;
                    }
                } else {
                    // Determinant too small, increase lambda
                    lambda = lambda.max(1e-8) * 2.0;
                }
            }
            
            if !step_accepted {
                // If no step accepted, use gradient descent with adaptive learning rate
                // Use larger learning rate initially to allow bigger steps
                let base_learning_rate = 0.5;  // Increased from 0.1
                let learning_rate = base_learning_rate / (1.0 + iteration as f64 * 0.05);  // Slower decay
                
                // Apply gradient descent step
                let mut new_a = a - learning_rate * grad_cost_a;
                let mut new_scale = scale - learning_rate * grad_cost_scale;
                
                // Apply bounds with reflection
                if new_a < bounds_a.0 {
                    new_a = 2.0 * bounds_a.0 - new_a;
                } else if new_a > bounds_a.1 {
                    new_a = 2.0 * bounds_a.1 - new_a;
                }
                if new_scale < bounds_scale.0 {
                    new_scale = 2.0 * bounds_scale.0 - new_scale;
                } else if new_scale > bounds_scale.1 {
                    new_scale = 2.0 * bounds_scale.1 - new_scale;
                }
                
                // Ensure strict feasibility
                new_a = new_a.max(bounds_a.0 + 1e-10).min(bounds_a.1 - 1e-10);
                new_scale = new_scale.max(bounds_scale.0 + 1e-10).min(bounds_scale.1 - 1e-10);
                
                let mut new_error = 0.0;
                for &(d, mr_obs) in &binned {
                    let mr_pred = chi_square_distribution_simple(d, new_a, new_scale);
                    let residual = mr_pred - mr_obs;
                    new_error += residual * residual;
                }
                nfev += 1;
                
                if new_error < current_error {
                    a = new_a;
                    scale = new_scale;
                    // If gradient descent step is accepted, slightly increase trust radius
                    trust_radius = (trust_radius * 1.1).min(max_trust_radius);
                } else {
                    // If gradient descent also fails, reduce trust radius
                    trust_radius = (trust_radius * gamma1).max(min_trust_radius);
                }
            }
            
            if nfev >= max_nfev {
                break;
            }
        }
        
        // Use best parameters found during optimization
        a = best_a;
        scale = best_scale;
        
        // Calculate R² (matching Python script: sum of squared residuals)
        let mut ss_res = 0.0;
        let mut ss_tot = 0.0;
        let mean_rate = binned.iter().map(|(_, mr)| *mr).sum::<f64>() / binned.len() as f64;
        
        // Calculate residuals: y_pred - y_obs (matching Python's curve_fit)
        for &(d, mr_obs) in &binned {
            let mr_pred = chi_square_distribution_simple(d, a, scale);
            let residual = mr_pred - mr_obs;  // Python: y_pred - y_obs
            ss_res += residual * residual;
            ss_tot += (mr_obs - mean_rate).powi(2);
        }
        
        let r2 = if ss_tot > 0.0 {
            1.0 - (ss_res / ss_tot)
        } else {
            f64::NEG_INFINITY
        };
        
        if r2 > best_r2 {
            best_r2 = r2;
            best_params = Some((a, scale));
        }
    }
    */
    
    if let Some((a, scale)) = best_params {
        Some((a, scale, best_r2))
    } else {
        None
    }
}

/// Calculate correction factor for a single position (matches Python script exactly)
fn calculate_correction_factor(
    depth: f64,
    mutation_rate: f64,
    chi2_params: Option<(f64, f64)>,
    target_y: f64,
    weight_increase: f64,
    weight_decrease: f64,
) -> f64 {
    if chi2_params.is_none() || target_y <= 0.0 {
        return 1.0;
    }
    
    let (a, scale) = chi2_params.unwrap();
    
    // Calculate expected mutation rate at this depth (red line)
    let y_at_depth = chi_square_distribution_simple(depth, a, scale);
    
    // Calculate correction factor
    if y_at_depth > target_y {
        // Red line above gray line: depth too low, need to increase
        let distance = (y_at_depth - target_y) / target_y;
        let sqrt_distance = distance.sqrt();
        1.0 + sqrt_distance * weight_increase
    } else {
        // Red line below gray line: depth too high, need to decrease
        let base_ratio_sigmoid = 0.88 + 0.12 / (1.0 + (depth / 500.0).powi(2));
        
        let base_ratio = if depth > 1000.0 {
            let log_factor = (depth / 1000.0).log10();
            let decay_factor = 0.08;
            let ratio = base_ratio_sigmoid - log_factor * decay_factor;
            ratio.max(0.5)
        } else {
            base_ratio_sigmoid
        };
        
        let base_ratio = base_ratio.max(0.5).min(1.0);
        
        if y_at_depth > 0.0 && target_y > 0.0 {
            let normalized_y = y_at_depth / target_y;
            let log2_normalized = (1.0 + normalized_y).log2() / 2.0_f64.log2();
            base_ratio + (1.0 - base_ratio) * log2_normalized
        } else {
            base_ratio
        }
    }
}

/// Calculate effective depth with PCR bias correction using Chi-Square distribution model
/// This function accounts for PCR amplification bias where high sequencing depths
/// cause mutation rate dilution due to preferential amplification of unmodified fragments.
/// 
/// Note: This function is deprecated in favor of global Chi-Square fitting approach.
/// It uses simplified fixed parameters and is only used when global fitting is not available.
/// 
/// The model uses:
/// - Sigmoid-like function for medium depths: base_ratio = 0.88 + 0.12 / (1 + (depth/500)^2)
/// - Logarithmic decay for high depths (>1000): ratio decreases by 0.08 per 10x depth increase
/// - Minimum ratio: 0.5 (maximum 50% reduction)
fn calculate_effective_depth(
    read_depth: usize, 
    mutation_count: usize,
    weight_increase: f64,
    weight_decrease: f64,
) -> usize {
    // Apply correction to all depths
    // This matches the Python script behavior which fits Chi-Square distribution globally
    
    // Calculate mutation rate for this position
    let mutation_rate = if read_depth > 0 {
        mutation_count as f64 / read_depth as f64
    } else {
        0.0
    };
    
    // Skip correction if mutation_rate is 0 or invalid (matches Python script filtering)
    if mutation_rate <= 0.0 || mutation_rate > 0.1 {
        return read_depth;
    }
    
    // Estimate target mutation rate (simplified: use a fixed target based on typical values)
    // In the full Chi-Square fitting approach, this would be calculated from the fitted curve
    // Here we use a conservative estimate: median mutation rate is typically around 0.01-0.02
    let target_mutation_rate = 0.015; // Typical median mutation rate
    
    // Calculate expected mutation rate at this depth using Chi-Square-like model
    // Simplified Chi-Square distribution: y = a * chi2.pdf(x/scale, df=2)
    // For estimation, we use a simplified model based on depth
    let depth_f64 = read_depth as f64;
    
    // Simplified Chi-Square model parameters (estimated from typical data)
    // These would normally come from fitting, but here we use fixed values
    // Note: This function is deprecated in favor of global fitting approach
    // Scale is typically around 20% of typical depth (1000)
    let scale = 1000.0 * 0.2; // Typical scale value (200)
    let a = target_mutation_rate * 10.0; // Amplitude factor
    
    // Calculate expected mutation rate at this depth
    // Simplified: chi2.pdf(x/scale, df=2) = exp(-x/(2*scale)) / 2
    let x_scaled = depth_f64 / scale;
    let expected_mutation_rate = if x_scaled > 1e-10 {
        // Simplified Chi-Square PDF approximation for df=2
        // chi2.pdf(x, df=2) = exp(-x/2) / 2
        let chi2_pdf = (-x_scaled / 2.0).exp() / 2.0;
        a * chi2_pdf.max(0.0)
    } else {
        target_mutation_rate
    };
    
    // Calculate correction factor based on expected vs target mutation rate
    let correction_factor = if expected_mutation_rate > target_mutation_rate {
        // Expected rate is above target: depth is too low, need to increase
        // Use square root of distance for gentler correction
        let distance = (expected_mutation_rate - target_mutation_rate) / target_mutation_rate;
        let sqrt_distance = distance.sqrt();
        1.0 + sqrt_distance * weight_increase
    } else {
        // Expected rate is below target: depth is too high, need to decrease
        // Use sigmoid-like function with logarithmic decay for high depths
        let base_ratio_sigmoid = 0.88 + 0.12 / (1.0 + (depth_f64 / 500.0).powi(2));
        
        // Apply logarithmic decay for depths > 1000
        let base_ratio = if depth_f64 > 1000.0 {
            let log_factor = (depth_f64 / 1000.0).log10();
            let decay_factor = 0.08;
            let ratio = base_ratio_sigmoid - log_factor * decay_factor;
            ratio.max(0.5) // Minimum ratio: 0.5
        } else {
            base_ratio_sigmoid
        };
        
        // Blend with mutation rate information
        let normalized_y = if target_mutation_rate > 0.0 {
            expected_mutation_rate / target_mutation_rate
        } else {
            1.0
        };
        
        // Use log2 transformation for smoothing
        let log2_normalized = (1.0 + normalized_y).log2() / 2.0_f64.log2();
        
        // Final correction factor: blend base_ratio with log2_normalized
        base_ratio + (1.0 - base_ratio) * log2_normalized * weight_decrease
    };
    
    // Ensure correction factor is reasonable
    let final_correction = correction_factor.max(0.5).min(10.0);
    
    // Calculate effective depth
    let effective_depth = (read_depth as f64 * final_correction) as usize;
    
    // Ensure effective depth is at least 1
    effective_depth.max(1)
}

// Read FASTA file, return a hash map: TID -> sequence
fn load_reference_sequences(
    fasta_path: &str,
) -> Result<HashMap<String, Vec<char>>, Box<dyn Error>> {
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

/// Find matching BAM files using glob pattern
fn find_matching_bam_files(glob_pattern: &str) -> Result<Vec<String>, Box<dyn Error>> {
    let mut matching_files = Vec::new();
    
    // Use glob to match files
    for entry in glob(glob_pattern)? {
        match entry {
            Ok(path) => {
                // Only process .bam files
                if let Some(ext) = path.extension() {
                    if ext == "bam" {
                        if let Some(full_path) = path.to_str() {
                            matching_files.push(full_path.to_string());
                        }
                    }
                }
            }
            Err(e) => {
                // Ignore files that cannot be accessed
                eprintln!("Warning: Unable to access file: {}", e);
            }
        }
    }
    
    // Sort by filename
    matching_files.sort();
    
    if matching_files.is_empty() {
        return Err(format!("Error: No matching BAM files found (pattern: {})", glob_pattern).into());
    }
    
    Ok(matching_files)
}

/// Check if BAM file has an index, return error if not
fn check_bam_index(bam_path: &str, logger: &mut Option<&mut Logger>) -> Result<(), Box<dyn Error>> {
    let index_path = format!("{}.bai", bam_path);
    
    // Check if index file exists
    if Path::new(&index_path).exists() {
        if let Some(ref mut logger) = logger {
            logger.log(&format!("BAM index file already exists: {}", index_path))?;
        }
        return Ok(());
    }
    
    // If not exists, return error message
    let error_msg = format!(
        "Error: BAM index file does not exist: {}\nPlease use the following command to create index:\nsamtools index -b {} -o {}",
        index_path, bam_path, index_path
    );
    
    if let Some(ref mut logger) = logger {
        logger.log(&error_msg)?;
    }
    
    Err(error_msg.into())
}

/// Scan BAM file data distribution to identify regions with data
/// Optimized for single-cell data: uses more efficient scanning strategy
fn scan_bam_distribution(
    bam_path: &str,
    reference_sequences: &HashMap<String, Vec<char>>,
    num_threads: usize,
    sample_size: usize, // Number of reads to sample, 0 means scan all
) -> Result<HashMap<String, Vec<(u32, u32)>>, Box<dyn Error>> {
    use rust_htslib::bam::IndexedReader;
    
    let mut bam = IndexedReader::from_path(bam_path)?;
    bam.set_threads(1);
    
    let header = bam.header().to_owned();
    
    // Store intervals with data for each region (start, end)
    let mut region_intervals: HashMap<String, Vec<(u32, u32)>> = HashMap::new();
    
    // Initialize all regions
    for ref_name in reference_sequences.keys() {
        region_intervals.insert(ref_name.clone(), Vec::new());
    }
    
    // Use more efficient method: directly record intervals to avoid large memory usage
    // Scan reads and record coverage intervals
    let mut read_count = 0;
    let mut skipped_ref_count = 0;
    let mut invalid_tid_count = 0;
    let mut records = bam.records();
    
    // Use sparse interval recording to avoid large memory usage
    let mut temp_intervals: HashMap<String, Vec<(u32, u32)>> = HashMap::new();
    
    // Record encountered reference sequence names (for debugging)
    let mut encountered_refs: HashSet<String> = HashSet::new();
    let mut expected_refs: HashSet<String> = region_intervals.keys().cloned().collect();
    
    while let Some(Ok(record)) = records.next() {
        read_count += 1;
        
        // If sample size is limited, only process first sample_size reads
        if sample_size > 0 && read_count > sample_size {
            break;
        }
        
        let tid = record.tid();
        if tid < 0 {
            invalid_tid_count += 1;
            continue;
        }
        
        let ref_name = String::from_utf8_lossy(header.tid2name(tid as u32)).to_string();
        encountered_refs.insert(ref_name.clone());
        
        if !region_intervals.contains_key(&ref_name) {
            skipped_ref_count += 1;
            continue;
        }
        
        let start = record.pos() as u32;
        let read_len = record.seq().len() as u32;
        let end = start + read_len;
        
        // Add to temporary interval list
        let intervals = temp_intervals.entry(ref_name).or_insert_with(Vec::new);
        intervals.push((start, end));
    }
    
    // Debug information: output reference sequence matching status
    if read_count > 0 {
        eprintln!("Data distribution scan statistics ({}):", bam_path);
        eprintln!("  Total reads: {}", read_count);
        eprintln!("  Invalid tid count: {}", invalid_tid_count);
        eprintln!("  Skipped reference sequences: {}", skipped_ref_count);
        eprintln!("  Found reference sequences: {}", temp_intervals.len());
        eprintln!("  Reference sequences in FASTA: {}", expected_refs.len());
        eprintln!("  Reference sequences encountered in BAM: {}", encountered_refs.len());
        
        if skipped_ref_count > read_count / 2 && read_count > 100 {
            eprintln!("Warning: More than 50% of reads have reference sequence names not in FASTA file");
            eprintln!("  First 5 reference sequences in BAM: {:?}", encountered_refs.iter().take(5).collect::<Vec<_>>());
            eprintln!("  First 5 reference sequences in FASTA: {:?}", expected_refs.iter().take(5).collect::<Vec<_>>());
        }
    }
    
    // Merge intervals for each region
    for (ref_name, intervals) in temp_intervals {
        if intervals.is_empty() {
            continue;
        }
        
        // Sort by start position
        let mut sorted_intervals = intervals;
        sorted_intervals.sort_by_key(|(s, _)| *s);
        
        // Merge overlapping or close intervals (single-cell data is usually concentrated)
        let mut merged: Vec<(u32, u32)> = Vec::new();
        let mut current = sorted_intervals[0];
        
        for &(start, end) in sorted_intervals.iter().skip(1) {
            if start <= current.1 + 5000 {
                // Merge intervals within 5kb (single-cell data characteristic)
                current.1 = current.1.max(end);
            } else {
                merged.push(current);
                current = (start, end);
            }
        }
        merged.push(current);
        
        region_intervals.insert(ref_name, merged);
    }
    
    Ok(region_intervals)
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
        'T' | 'U' => return 8,  // U and T are equivalent in RNA/DNA
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
        effective_depth: 0,  // Quality-filtered depth
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
        effective_depth: 0,  // Quality-filtered depth
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
        if len <= head_skip + tail_skip + 10 {
            continue;
        }
        let start = head_skip;
        let end = len - tail_skip;
        let mut region: Vec<f64> = bd_vec[start..end].iter().cloned().collect();
        region.sort_by(|a, b| b.partial_cmp(a).unwrap()); // Descending order
        let n = region.len();
        let s = (n as f64 * norm_percentile_start).round() as usize;
        let e = (n as f64 * norm_percentile_end).round() as usize;
        if s >= e || e > n {
            continue;
        }
        let slice = &region[s..e];
        let scaling = if use_median {
            let mut v = slice.to_vec();
            v.sort_by(|a, b| a.partial_cmp(b).unwrap());
            let m = v.len();
            if m == 0 {
                1.0
            } else if m % 2 == 0 {
                (v[m / 2 - 1] + v[m / 2]) / 2.0
            } else {
                v[m / 2]
            }
        } else {
            let sum: f64 = slice.iter().sum();
            if slice.len() == 0 {
                1.0
            } else {
                sum / slice.len() as f64
            }
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
    eprintln!("=== PROGRAM STARTED ===");
    let cli = Cli::parse();
    eprintln!("=== MAIN FUNCTION ENTERED ===");
    eprintln!("DEBUG: Cli = {:?}", cli);
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
            logger.log(&format!(
                "Runtime: {}",
                chrono::Utc::now().format("%Y-%m-%d %H:%M:%S")
            ))?;
            logger.log(&format!("BAM File: {}", args.bam))?;
            logger.log(&format!("FASTA File: {}", args.fasta))?;
            let strand_filter = normalize_strand_filter(&args.strand);
            logger.log(&format!("Output File: {}", args.output))?;
            logger.log(&format!("Threads: {}", num_threads))?;
            logger.log(&format!("Strand Filter: {}", strand_filter))?;
            
            // Handle single-cell unified processing mode
            if args.single_cell {
                logger.log("=== Single-cell unified mode enabled ===")?;
                logger.log(&format!("Glob pattern: {}", args.bam))?;
                
                // Initialize thread pool before batch processing (initialize only once)
                rayon::ThreadPoolBuilder::new()
                    .num_threads(num_threads)
                    .build_global()
                    .ok(); // Ignore if already initialized
                
                // Find matching BAM files using glob pattern
                logger.log("Finding matching BAM files...")?;
                let bam_files = find_matching_bam_files(&args.bam)?;
                
                logger.log(&format!("Found {} matching BAM file(s)", bam_files.len()))?;
                
                // Validate all BAM files
                logger.log("Validating BAM files and indexes...")?;
                for (idx, bam_file) in bam_files.iter().enumerate() {
                    logger.log(&format!("  [{}/{}] Checking: {}", idx + 1, bam_files.len(), bam_file))?;
                    let mut logger_opt = Some(&mut logger);
                    check_bam_index(bam_file, &mut logger_opt)?;
                }
                
                // Check output directory
                let output_dir = Path::new(&args.output);
                if !output_dir.exists() {
                    fs::create_dir_all(output_dir)?;
                    logger.log(&format!("Created output directory: {}", args.output))?;
                } else if !output_dir.is_dir() {
                    return Err(format!("Error: Output path is not a directory: {}", args.output).into());
                }
                
                // Execute unified processing strategy
                let mut logger_opt = Some(&mut logger);
                let result =                 // Note: single-cell mode PCR bias correction is handled in process_single_cell_unified
                // which needs to be updated to accept and use these parameters
                process_single_cell_unified(
                    &bam_files,
                    &args.fasta,
                    output_dir,
                    num_threads,
                    args.window,
                    strand_filter,
                    args.min_base_qual,
                    args.pcr_bias_correction,
                    args.weight_increase,
                    args.weight_decrease,
                    &mut logger_opt,
                );
                
                match &result {
                    Ok(_) => logger.log("Single-cell unified processing completed")?,
                    Err(e) => logger.log(&format!("Single-cell unified processing failed: {}", e))?,
                }
                
                result
            } else if args.batch {
                // Handle batch mode (process files one by one)
                logger.log("=== Batch mode enabled ===")?;
                logger.log(&format!("Glob pattern: {}", args.bam))?;
                
                // Initialize thread pool before batch processing (initialize only once)
                rayon::ThreadPoolBuilder::new()
                    .num_threads(num_threads)
                    .build_global()
                    .ok(); // Ignore if already initialized
                
                // Find matching BAM files using glob pattern
                logger.log("Finding matching BAM files...")?;
                let bam_files = find_matching_bam_files(&args.bam)?;
                
                logger.log(&format!("Found {} matching BAM file(s)", bam_files.len()))?;
                
                // Validate all BAM files
                logger.log("Validating BAM files and indexes...")?;
                for (idx, bam_file) in bam_files.iter().enumerate() {
                    logger.log(&format!("  [{}/{}] Checking: {}", idx + 1, bam_files.len(), bam_file))?;
                    let mut logger_opt = Some(&mut logger);
                    check_bam_index(bam_file, &mut logger_opt)?;
                }
                
                // Check output directory
                let output_dir = Path::new(&args.output);
                if !output_dir.exists() {
                    fs::create_dir_all(output_dir)?;
                    logger.log(&format!("Created output directory: {}", args.output))?;
                } else if !output_dir.is_dir() {
                    return Err(format!("Error: Output path is not a directory: {}", args.output).into());
                }
                
                // Batch process all matching BAM files
                logger.log("Starting batch processing...")?;
                let mut success_count = 0;
                let mut fail_count = 0;
                
                for (idx, bam_file) in bam_files.iter().enumerate() {
                    logger.log(&format!("=== Processing [{}/{}] ===", idx + 1, bam_files.len()))?;
                    logger.log(&format!("BAM File: {}", bam_file))?;
                    
                    // Generate output filename
                    let bam_name = Path::new(bam_file)
                        .file_stem()
                        .and_then(|s| s.to_str())
                        .unwrap_or("unknown");
                    let output_file = output_dir.join(format!("{}.csv", bam_name));
                    let output_str = output_file.to_str().unwrap();
                    
                    // Generate log filename
                    let log_file_path = if let Some(log_path) = &args.log {
                        let log_dir = Path::new(log_path).parent().unwrap_or(Path::new("."));
                        if !log_dir.exists() {
                            fs::create_dir_all(log_dir)?;
                        }
                        log_dir.join(format!("{}.log", bam_name)).to_str().unwrap().to_string()
                    } else {
                        format!("{}.log", bam_name)
                    };
                    
                    logger.log(&format!("Output: {}", output_str))?;
                    logger.log(&format!("Log: {}", log_file_path))?;
                    
                    // Process single file
                    let mut file_logger_opt = Some(&mut logger);
                    match process_sample_pileup_ultra_fast(
                        bam_file,
                        &args.fasta,
                        output_str,
                        num_threads,
                        args.window,
                        strand_filter,
                        args.min_base_qual,
                        args.pcr_bias_correction,
                        args.weight_increase,
                        args.weight_decrease,
                        &mut file_logger_opt,
                    ) {
                        Ok(_) => {
                            logger.log(&format!("✓ Successfully processed: {}", bam_file))?;
                            success_count += 1;
                        }
                        Err(e) => {
                            logger.log(&format!("✗ Failed to process {}: {}", bam_file, e))?;
                            fail_count += 1;
                        }
                    }
                }
                
                logger.log("=== Batch processing completed ===")?;
                logger.log(&format!("Success: {}, Failed: {}", success_count, fail_count))?;
                
                if fail_count > 0 {
                    return Err(format!("Processing completed, but {} file(s) failed", fail_count).into());
                }
                
                Ok(())
            } else {
                // Standard mode: process single BAM file
                logger.log("=== Standard mode ===")?;
                
                // Check if BAM index exists
                logger.log("Checking BAM index...")?;
                let mut logger_opt = Some(&mut logger);
                check_bam_index(&args.bam, &mut logger_opt)?;
                
                logger.log("Starting BAM file processing...")?;

                let mut logger_opt = Some(&mut logger);
                let result = process_sample_pileup_ultra_fast(
                    &args.bam,
                    &args.fasta,
                    &args.output,
                    num_threads,
                    args.window,
                    strand_filter,
                    args.min_base_qual,
                    args.pcr_bias_correction,
                    args.weight_increase,
                    args.weight_decrease,
                    &mut logger_opt,
                );
                
                match &result {
                    Ok(_) => logger.log("BAM file processing completed")?,
                    Err(e) => logger.log(&format!("BAM file processing failed: {}", e))?,
                }

                result
            }
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
                    &args.reactive_bases,
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
                // Full evaluation using base matching and auto-shift correction
                evaluate::evaluate_reactivity_accuracy_with_base_matching_main(
                    &args.reactivity_file,
                    &args.structure_file,
                    &args.output_dir,
                    &args.signal_type,
                    &args.gene_id,
                    &args.strand,
                    &args.reactive_bases,
                    &mut logger,
                )
            } else if args.auto_shift {
                // Only use auto-shift correction
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
                // Basic evaluation (without base matching and auto-shift correction)
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
            // Validate compare command parameters
            validate_compare_args(&args)?;

            let log_file = if let Some(log_path) = &args.log {
                std::fs::File::create(log_path)?
            } else {
                std::fs::File::create("compare.log")?
            };
            let mut logger = Logger::new(log_file);

            // Parse statistical test type
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
            eprintln!("=== PLOT COMMAND ENTERED ===");
            eprintln!("DEBUG: PlotArgs = {:?}", args);
            // Validate plot command parameters
            validate_plot_args(&args)?;

            let log_file = if let Some(log_path) = &args.log {
                std::fs::File::create(log_path)?
            } else {
                std::fs::File::create("plot.log")?
            };
            let mut logger = Logger::new(log_file);

            // Check if only SVG plotting is desired
            eprintln!("DEBUG: svg_template = {:?}", args.svg_template);
            eprintln!("DEBUG: reactivity_file = {:?}", args.reactivity_file);
            let svg_only_mode = args.svg_template.is_some() && args.reactivity_file.is_some();
            eprintln!("SVG-only mode: {}", svg_only_mode);
            eprintln!("This is a test message to verify the code is running");
            eprintln!("About to start SVG processing");

            if !svg_only_mode {
                // Execute regular plot functionality
                let mod_csv = args.mod_csv.as_ref().unwrap(); // Already validated in validate_plot_args
                let unmod_csv = args.unmod_csv.as_ref().unwrap(); // Already validated in validate_plot_args

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

            // If SVG template is provided, execute SVG plotting
            eprintln!("=== CHECKING SVG TEMPLATE CONDITION ===");
            eprintln!("svg_template: {:?}", args.svg_template);
            if let Some(svg_template) = &args.svg_template {
                eprintln!("=== SVG TEMPLATE FOUND ===");
                eprintln!("reactivity_file: {:?}", args.reactivity_file);
                if let Some(reactivity_file) = &args.reactivity_file {
                    eprintln!("=== REACTIVITY FILE FOUND ===");
                    logger.log("Starting SVG reactivity plotting...")?;

                    // Ensure output directory exists
                    std::fs::create_dir_all(&args.output)?;

                    // Check if interactive mode is requested
                    if args.interactive {
                        logger.log("Generating interactive HTML visualization...")?;
                        
                        if args.svg_signal == "all" {
                            // For "all" signal type, generate interactive HTML for each detected signal
                            logger.log("Auto-detecting signal types for interactive visualization...")?;
                            // First generate static SVG to get the structure, then convert to interactive HTML
                            let temp_svg_dir = format!("{}/temp_svg", args.output);
                            std::fs::create_dir_all(&temp_svg_dir)?;
                            
                            plot::plot_multiple_signals_to_svg(
                                reactivity_file,
                                svg_template,
                                &temp_svg_dir,
                                &args.svg_bases,
                                &args.svg_strand,
                                args.svg_ref.as_deref(),
                                args.svg_max_shift,
                                None,
                            )?;
                            
                            // Generate interactive HTML for each SVG file
                            let svg_files = std::fs::read_dir(&temp_svg_dir)?
                                .filter_map(|e| e.ok())
                                .filter(|e| e.path().extension().map(|s| s == "svg").unwrap_or(false))
                                .collect::<Vec<_>>();
                            
                            for svg_file in svg_files {
                                let svg_path = svg_file.path();
                                let signal_name = svg_path.file_stem()
                                    .and_then(|s| s.to_str())
                                    .and_then(|s| s.strip_prefix("rna_structure_colored_"))
                                    .unwrap_or("unknown");
                                
                                let html_output = format!(
                                    "{}/rna_structure_interactive_{}.html",
                                    args.output, signal_name
                                );
                                
                                plot::generate_interactive_visualization(
                                    reactivity_file,
                                    svg_path.to_str().unwrap(),
                                    &html_output,
                                    &args.svg_bases,
                                    signal_name,
                                    &args.svg_strand,
                                    None, // Use default color ranges
                                    args.svg_circle_filled,
                                )?;
                                
                                logger.log(&format!("Interactive HTML created: {}", html_output))?;
                            }
                            
                            // Clean up temporary SVG files
                            std::fs::remove_dir_all(&temp_svg_dir)?;
                        } else {
                            // Generate interactive HTML for specific signal type
                            let html_output = format!(
                                "{}/rna_structure_interactive_{}.html",
                                args.output, args.svg_signal
                            );
                            
                            // First generate a temporary SVG to get the structure with colors
                            let temp_svg = format!("{}/temp_{}.svg", args.output, args.svg_signal);
                            let mut config = plot::SvgDrawConfig::default();
                            config.circle_filled = args.svg_circle_filled;
                            config.font_color = args.svg_font_color.clone();
                            config.legend_item_width = args.svg_legend_width;
                            config.legend_item_height = args.svg_legend_height;
                            
                            plot::draw_svg_with_custom_layers(
                                reactivity_file,
                                svg_template,
                                &temp_svg,
                                &args.svg_bases,
                                &args.svg_signal,
                                &args.svg_strand,
                                None,
                                config,
                                args.svg_ref.as_deref(),
                                args.svg_max_shift,
                            )?;
                            
                            // Generate interactive HTML from the SVG
                            plot::generate_interactive_visualization(
                                reactivity_file,
                                &temp_svg,
                                &html_output,
                                &args.svg_bases,
                                &args.svg_signal,
                                &args.svg_strand,
                                None, // Use default color ranges
                                args.svg_circle_filled,
                            )?;
                            
                            // Clean up temporary SVG
                            std::fs::remove_file(&temp_svg)?;
                            
                            logger.log(&format!("Interactive HTML visualization created: {}", html_output))?;
                        }
                    } else {
                        // Standard SVG plotting mode
                    // Decide plotting method based on signal type
                    if args.svg_signal == "all" {
                        // Plot all signal types
                        logger.log("Auto-detecting signal types and plotting all...")?;
                        plot::plot_multiple_signals_to_svg(
                            reactivity_file,
                            svg_template,
                            &args.output,
                            &args.svg_bases,
                            &args.svg_strand,
                            args.svg_ref.as_deref(),
                            args.svg_max_shift,
                            None, // Use default color range
                        )?;
                        logger.log("All signal types plotted successfully")?;
                    } else {
                        // Plot specified signal type
                        eprintln!("=== ENTERING SPECIFIC SIGNAL TYPE BRANCH ===");
                        logger.log(&format!("Plotting signal type: {}", args.svg_signal))?;
                        let svg_output = format!(
                            "{}/rna_structure_colored_{}.svg",
                            args.output, args.svg_signal
                        );

                        // Force use new architecture for testing
                        eprintln!("=== FORCING NEW ARCHITECTURE BRANCH ===");
                        eprintln!(
                            "Calling NEW ARCHITECTURE draw_svg_with_custom_layers with signal: {}",
                            args.svg_signal
                        );
                        let mut config = plot::SvgDrawConfig::default();
                        config.circle_filled = args.svg_circle_filled; // If flag is specified, use filled circles
                        config.font_color = args.svg_font_color.clone(); // Set font color
                        config.legend_item_width = args.svg_legend_width;
                        config.legend_item_height = args.svg_legend_height;
                        let result = plot::draw_svg_with_custom_layers(
                            reactivity_file,
                            svg_template,
                            &svg_output,
                            &args.svg_bases,
                            &args.svg_signal,
                            &args.svg_strand,
                            None, // Use default color range
                            config,
                            args.svg_ref.as_deref(),
                            args.svg_max_shift,
                        );
                        eprintln!("draw_svg_with_custom_layers returned: {:?}", result);
                        result?;
                        logger.log(&format!(
                            "SVG plot created with new architecture: {}",
                            svg_output
                        ))?;
                        }
                    }
                } else {
                    logger.log("Warning: SVG template provided but no reactivity file specified. Skipping SVG plotting.")?;
                }
            }

            Ok(())
        }
        Commands::Duet(args) => {
            duet::validate_args(&args)?;
            let log_file = if let Some(log_path) = &args.log {
                std::fs::File::create(log_path)?
            } else {
                std::fs::File::create("duet.log")?
            };
            let mut logger = Logger::new(log_file);
            duet::analyze_reactivity(&args, &mut logger)
        }
        Commands::Extract(args) => {
            let log_file = if let Some(log_path) = &args.log {
                std::fs::File::create(log_path)?
            } else {
                std::fs::File::create("extract.log")?
            };
            let mut logger = Logger::new(log_file);
            
            logger.log("=== ModDetector Extract Function Log ===")?;
            logger.log(&format!("Software Version: v{}", VERSION))?;
            logger.log(&format!(
                "Runtime: {}",
                chrono::Utc::now().format("%Y-%m-%d %H:%M:%S")
            ))?;
            logger.log(&format!("Input File: {}", args.input))?;
            logger.log(&format!("Input Format: {}", args.input_format))?;
            logger.log(&format!("GTF File: {}", args.gtf))?;
            logger.log(&format!("Target Gene: {}", args.target_gene))?;
            logger.log(&format!("Output Prefix: {}", args.output))?;
            logger.log(&format!("Relative Position: {}", args.relative_position))?;
            
            // Set rayon thread pool size
            if args.threads > 1 {
                rayon::ThreadPoolBuilder::new()
                    .num_threads(args.threads)
                    .build_global()
                    .unwrap_or_default();
            }
            
            if args.input_format == "bam-readcount" {
                extract::extract_gene_from_bamreadcount(
                    &args.input,
                    &args.gtf,
                    &args.target_gene,
                    &args.output,
                    args.threads,
                    args.relative_position,
                )?;
            } else {
                extract::extract_gene_from_count(
                    &args.input,
                    &args.gtf,
                    &args.target_gene,
                    &args.output,
                    args.threads,
                    args.relative_position,
                )?;
            }
            
            logger.log("Gene extraction completed successfully - each gene region saved to separate files")?;
            Ok(())
        }
        Commands::Correct(args) => {
            let log_file = std::fs::File::create("correct.log")?;
            let mut logger = Logger::new(log_file);
            logger.log("=== PROGRAM STARTED ===")?;
            logger.log("=== MAIN FUNCTION ENTERED ===")?;
            correct::apply_pcr_bias_correction(&args, &mut logger)?;
            Ok(())
        }
        Commands::Convert(args) => {
            let log_file = if let Some(log_path) = &args.log {
                std::fs::File::create(log_path)?
            } else {
                std::fs::File::create("convert.log")?
            };
            let mut logger = Logger::new(log_file);
            
            convert::convert_bamreadcount_to_pileup(&args, &mut logger)?;
            Ok(())
        }
    }
}

/// Validate count command parameters
fn validate_count_args(args: &CountArgs) -> Result<(), Box<dyn Error>> {
    // Validate BAM file or single-cell mode input
    if args.bam.trim().is_empty() {
        return Err("Error: BAM file path cannot be empty".into());
    }
    
    if args.single_cell || args.batch {
        // Single-cell unified mode or batch mode: validate glob pattern
        // Check if it contains wildcards
        if !args.bam.contains('*') && !args.bam.contains('?') {
            return Err("Error: Batch/single-cell mode requires glob pattern (containing * or ?), e.g.: /path/to/bam/*sort.bam".into());
        }
        
        // Validate if glob pattern is valid (try to find files)
        match find_matching_bam_files(&args.bam) {
            Ok(files) => {
                if files.is_empty() {
                    return Err(format!("Error: No matching BAM files found (pattern: {})", args.bam).into());
                }
            }
            Err(e) => {
                return Err(format!("Error: Glob pattern is invalid or inaccessible: {}", e).into());
            }
        }
        
        // Validate output is a directory
        if Path::new(&args.output).exists() && !Path::new(&args.output).is_dir() {
            return Err(format!("Error: Output path must be a directory in batch/single-cell mode: {}", args.output).into());
        }
    } else {
        // Standard mode: validate single BAM file
        if !Path::new(&args.bam).exists() {
            return Err(format!("Error: BAM file does not exist: {}", args.bam).into());
        }
        if !args.bam.ends_with(".bam") {
            return Err(format!("Error: BAM file path must end with .bam: {}", args.bam).into());
        }
    }

    // Validate FASTA file
    if args.fasta.trim().is_empty() {
        return Err("Error: FASTA file path cannot be empty".into());
    }
    if !Path::new(&args.fasta).exists() {
        return Err(format!("Error: FASTA file does not exist: {}", args.fasta).into());
    }
    if !args.fasta.ends_with(".fa") && !args.fasta.ends_with(".fasta") {
        return Err(format!("Error: FASTA file path must end with .fa or .fasta: {}", args.fasta).into());
    }

    // Validate output path
    if args.output.trim().is_empty() {
        return Err("Error: Output path cannot be empty".into());
    }
    // In batch/single-cell mode output must be a directory, in standard mode output must be a CSV file
    if args.single_cell || args.batch {
        // Batch/single-cell mode: output must be a directory (if exists)
        if Path::new(&args.output).exists() && !Path::new(&args.output).is_dir() {
            return Err(format!("Error: Output path must be a directory in batch/single-cell mode: {}", args.output).into());
        }
    } else {
        // Standard mode: output must be a CSV file
        if !args.output.ends_with(".csv") {
            return Err(format!("Error: Output file path must end with .csv: {}", args.output).into());
        }
    }

    // Validate thread count
    if let Some(threads) = args.threads {
        if threads == 0 {
            return Err("Error: Thread count cannot be 0".into());
        }
        // Remove hard limit, only provide suggestion
        if threads > 128 {
            eprintln!("Warning: Thread count {} may exceed optimal range, recommend using 64-128 threads for best performance", threads);
        }
    }
    // Validate window size (optional)
    if let Some(w) = args.window {
        if w == 0 {
            return Err("Error: Window size cannot be 0".into());
        }
        if w > 10_000_000 {
            return Err(format!("Error: Window size cannot exceed 10000000 (current: {})", w).into());
        }
    }

    // Validate strand direction parameter
    match args.strand.trim() {
        "+" | "-" | "+/-" | "both" | "Both" => {}
        other => {
            return Err(format!("Error: Invalid strand direction parameter: {}. Supported values: +, -, +/-", other).into())
        }
    }

    Ok(())
}

fn normalize_strand_filter(filter: &str) -> &'static str {
    match filter.trim() {
        "+" | "plus" | "Plus" => "+",
        "-" | "minus" | "Minus" => "-",
        "+/-" | "+/−" | "both" | "Both" | "BOTH" => "both",
        _ => "both",
    }
}

/// Validate compare command parameters
fn validate_compare_args(args: &CompareArgs) -> Result<(), Box<dyn Error>> {
    // Validate comparison mode
    match args.mode.as_str() {
        "mod-vs-unmod" => {
            // Validate modified sample file
            let mod_csv = args.mod_csv.as_ref().ok_or("Error: mod-vs-unmod mode requires modified sample file to be specified")?;
            if mod_csv.trim().is_empty() {
                return Err("Error: Modified sample file path cannot be empty".into());
            }
            if !Path::new(mod_csv).exists() {
                return Err(format!("Error: Modified sample file does not exist: {}", mod_csv).into());
            }
            if !mod_csv.ends_with(".csv") {
                return Err(format!("Error: Modified sample file path must end with .csv: {}", mod_csv).into());
            }

            // Validate unmodified sample file
            let unmod_csv = args.unmod_csv.as_ref().ok_or("Error: mod-vs-unmod mode requires unmodified sample file to be specified")?;
            if unmod_csv.trim().is_empty() {
                return Err("Error: Unmodified sample file path cannot be empty".into());
            }
            if !Path::new(unmod_csv).exists() {
                return Err(format!("Error: Unmodified sample file does not exist: {}", unmod_csv).into());
            }
            if !unmod_csv.ends_with(".csv") {
                return Err(format!("Error: Unmodified sample file path must end with .csv: {}", unmod_csv).into());
            }
        },
        "biological-replicates" => {
            // Validate group1 files
            let group1_files = args.group1_files.as_ref().ok_or("Error: biological-replicates mode requires group1 files to be specified")?;
            if group1_files.trim().is_empty() {
                return Err("Error: Group1 file list cannot be empty".into());
            }
            for file in group1_files.split(',') {
                let file = file.trim();
                if file.is_empty() {
                    continue;
                }
                if !Path::new(file).exists() {
                    return Err(format!("Error: Group1 file does not exist: {}", file).into());
                }
                if !file.ends_with(".csv") {
                    return Err(format!("Error: Group1 file path must end with .csv: {}", file).into());
                }
            }

            // Validate group2 files
            let group2_files = args.group2_files.as_ref().ok_or("Error: biological-replicates mode requires group2 files to be specified")?;
            if group2_files.trim().is_empty() {
                return Err("Error: Group2 file list cannot be empty".into());
            }
            for file in group2_files.split(',') {
                let file = file.trim();
                if file.is_empty() {
                    continue;
                }
                if !Path::new(file).exists() {
                    return Err(format!("Error: Group2 file does not exist: {}", file).into());
                }
                if !file.ends_with(".csv") {
                    return Err(format!("Error: Group2 file path must end with .csv: {}", file).into());
                }
            }
        },
        "reactivity-groups" => {
            // Validate reactivity group1 files
            let group1_files = args.reactivity_group1.as_ref().ok_or("Error: reactivity-groups mode requires group1 files to be specified")?;
            if group1_files.trim().is_empty() {
                return Err("Error: Reactivity group1 file list cannot be empty".into());
            }
            for file in group1_files.split(',') {
                let file = file.trim();
                if file.is_empty() {
                    continue;
                }
                if !Path::new(file).exists() {
                    return Err(format!("Error: Reactivity group1 file does not exist: {}", file).into());
                }
                if !file.ends_with(".csv") {
                    return Err(format!("Error: Reactivity group1 file path must end with .csv: {}", file).into());
                }
            }

            // Validate reactivity group2 files
            let group2_files = args.reactivity_group2.as_ref().ok_or("Error: reactivity-groups mode requires group2 files to be specified")?;
            if group2_files.trim().is_empty() {
                return Err("Error: Reactivity group2 file list cannot be empty".into());
            }
            for file in group2_files.split(',') {
                let file = file.trim();
                if file.is_empty() {
                    continue;
                }
                if !Path::new(file).exists() {
                    return Err(format!("Error: Reactivity group2 file does not exist: {}", file).into());
                }
                if !file.ends_with(".csv") {
                    return Err(format!("Error: Reactivity group2 file path must end with .csv: {}", file).into());
                }
            }
        },
        "reactivity-results" => {
            // Validate first reactivity file
            let reactivity1 = args.reactivity1.as_ref().ok_or("Error: reactivity-results mode requires first reactivity file to be specified")?;
            if reactivity1.trim().is_empty() {
                return Err("Error: First reactivity file path cannot be empty".into());
            }
            if !Path::new(reactivity1).exists() {
                return Err(format!("Error: First reactivity file does not exist: {}", reactivity1).into());
            }
            if !reactivity1.ends_with(".csv") {
                return Err(format!("Error: First reactivity file path must end with .csv: {}", reactivity1).into());
            }

            // Validate second reactivity file
            let reactivity2 = args.reactivity2.as_ref().ok_or("Error: reactivity-results mode requires second reactivity file to be specified")?;
            if reactivity2.trim().is_empty() {
                return Err("Error: Second reactivity file path cannot be empty".into());
            }
            if !Path::new(reactivity2).exists() {
                return Err(format!("Error: Second reactivity file does not exist: {}", reactivity2).into());
            }
            if !reactivity2.ends_with(".csv") {
                return Err(format!("Error: Second reactivity file path must end with .csv: {}", reactivity2).into());
            }
        },
        _ => return Err("Error: Invalid comparison mode. Supported modes: mod-vs-unmod, biological-replicates, reactivity-groups, reactivity-results".into()),
    }

    // Validate output path
    if args.output.trim().is_empty() {
        return Err("Error: Output file path cannot be empty".into());
    }
    if !args.output.ends_with(".csv") {
        return Err(format!("Error: Output file path must end with .csv: {}", args.output).into());
    }

    // Validate minimum depth
    if args.min_depth == 0 {
        return Err("Error: Minimum depth cannot be 0".into());
    }
    if args.min_depth > 10000 {
        return Err(format!("Error: Minimum depth cannot exceed 10000 (current: {})", args.min_depth).into());
    }

    // Validate minimum fold change
    if args.min_fold <= 0.0 {
        return Err("Error: Minimum fold change must be greater than 0".into());
    }
    if args.min_fold > 1000.0 {
        return Err(format!(
            "Error: Minimum fold change cannot exceed 1000 (current: {})",
            args.min_fold
        )
        .into());
    }

    // Validate statistical test type
    match args.test_type.as_str() {
        "t-test" | "mann-whitney" | "wilcoxon" | "chi-square" | "continuity" | "diffscan" | "deltashape" => {},
        _ => return Err("Error: Invalid statistical test type. Supported tests: t-test, mann-whitney, wilcoxon, chi-square, continuity, diffscan, deltashape".into()),
    }

    // Validate significance level
    if args.alpha <= 0.0 || args.alpha >= 1.0 {
        return Err("Error: Significance level must be between 0 and 1".into());
    }

    Ok(())
}

/// Validate evaluate command parameters
fn validate_evaluate_args(args: &EvaluateArgs) -> Result<(), Box<dyn Error>> {
    // Validate reactivity file
    if args.reactivity_file.trim().is_empty() {
        return Err("Error: Reactivity file path cannot be empty".into());
    }
    if !Path::new(&args.reactivity_file).exists() {
        return Err(format!("Error: Reactivity file does not exist: {}", args.reactivity_file).into());
    }
    if !args.reactivity_file.ends_with(".csv") {
        return Err(format!(
            "Error: Reactivity file path must end with .csv: {}",
            args.reactivity_file
        )
        .into());
    }

    // Validate secondary structure file
    if args.structure_file.trim().is_empty() {
        return Err("Error: Secondary structure file path cannot be empty".into());
    }
    if !Path::new(&args.structure_file).exists() {
        return Err(format!("Error: Secondary structure file does not exist: {}", args.structure_file).into());
    }
    if !args.structure_file.ends_with(".dp") && !args.structure_file.ends_with(".ct") {
        return Err(format!(
            "Error: Secondary structure file path must end with .dp or .ct: {}",
            args.structure_file
        )
        .into());
    }

    // Validate output directory
    if args.output_dir.trim().is_empty() {
        return Err("Error: Output directory path cannot be empty".into());
    }

    // Validate signal type
    if args.signal_type != "stop" && args.signal_type != "mutation" && args.signal_type != "both" {
        return Err(format!(
            "Error: Signal type must be 'stop', 'mutation', or 'both', current: {}",
            args.signal_type
        )
        .into());
    }

    // Validate reactive-bases: only A/C/G/T/U allowed
    if !args.reactive_bases.chars().all(|c| "ACGTU".contains(c)) {
        return Err(format!(
            "Error: reactive-bases must contain only A,C,G,T,U; current: {}",
            args.reactive_bases
        )
        .into());
    }
    if args.reactive_bases.is_empty() {
        return Err("Error: reactive-bases cannot be empty".into());
    }

    // Validate gene ID
    if args.gene_id.trim().is_empty() {
        return Err("Error: Gene ID cannot be empty".into());
    }

    // Validate strand information
    if args.strand != "+" && args.strand != "-" {
        return Err(format!("Error: Strand must be '+' or '-', current: {}", args.strand).into());
    }

    // Validate parameter combinations
    if args.optimized && (args.use_base_matching || args.auto_shift) {
        return Err("Error: Optimized mode cannot be used with base-matching or auto-shift parameters".into());
    }

    if !args.use_base_matching && args.auto_shift {
        return Err("Error: auto-shift parameter requires base-matching parameter to be used together".into());
    }

    Ok(())
}

/// Validate plot command parameters
fn validate_plot_args(args: &PlotArgs) -> Result<(), Box<dyn Error>> {
    // Check if only SVG plotting is desired
    eprintln!("VALIDATE: svg_template = {:?}", args.svg_template);
    eprintln!("VALIDATE: reactivity_file = {:?}", args.reactivity_file);
    let svg_only_mode = args.svg_template.is_some() && args.reactivity_file.is_some();
    eprintln!("VALIDATE: svg_only_mode = {}", svg_only_mode);

    // If not SVG-only mode, validate modified sample file
    if !svg_only_mode {
        let mod_csv = args
            .mod_csv
            .as_ref()
            .ok_or("Error: Modified sample file path cannot be empty")?;
        if mod_csv.trim().is_empty() {
            return Err("Error: Modified sample file path cannot be empty".into());
        }
        if !Path::new(mod_csv).exists() {
            return Err(format!("Error: Modified sample file does not exist: {}", mod_csv).into());
        }
        if !mod_csv.ends_with(".csv") {
            return Err(format!("Error: Modified sample file path must end with .csv: {}", mod_csv).into());
        }

        // Validate unmodified sample file
        let unmod_csv = args
            .unmod_csv
            .as_ref()
            .ok_or("Error: Unmodified sample file path cannot be empty (except in SVG-only mode)")?;
        if unmod_csv.trim().is_empty() {
            return Err("Error: Unmodified sample file path cannot be empty".into());
        }
        if !Path::new(unmod_csv).exists() {
            return Err(format!("Error: Unmodified sample file does not exist: {}", unmod_csv).into());
        }
        if !unmod_csv.ends_with(".csv") {
            return Err(format!("Error: Unmodified sample file path must end with .csv: {}", unmod_csv).into());
        }
    }

    // Validate output directory
    if args.output.trim().is_empty() {
        return Err("Error: Output directory path cannot be empty".into());
    }

    // Validate coverage threshold
    if let Some(coverage) = args.coverage_threshold {
        if coverage < 0.0 || coverage > 1.0 {
            return Err(format!("Error: Coverage threshold must be between 0.0 and 1.0, current: {}", coverage).into());
        }
    }

    // Validate depth threshold
    if let Some(depth) = args.depth_threshold {
        if depth <= 0.0 {
            return Err("Error: Depth threshold must be greater than 0".into());
        }
        if depth > 10000.0 {
            return Err(format!("Error: Depth threshold cannot exceed 10000, current: {}", depth).into());
        }
    }

    // Validate thread count
    if let Some(threads) = args.threads {
        if threads == 0 {
            return Err("Error: Thread count cannot be 0".into());
        }
        // Remove hard limit, only provide suggestion
        if threads > 128 {
            eprintln!("Warning: Thread count {} may exceed optimal range, recommend using 64-128 threads for best performance", threads);
        }
    }

    // Validate optional reactivity file
    if let Some(reactivity_file) = &args.reactivity_file {
        if !Path::new(reactivity_file).exists() {
            return Err(format!("Error: Reactivity file does not exist: {}", reactivity_file).into());
        }
        if !reactivity_file.ends_with(".csv") {
            return Err(format!(
                "Error: Reactivity file path must end with .csv: {}",
                reactivity_file
            )
            .into());
        }
    }

    // Validate GFF file
    if let Some(gff) = &args.gff {
        if !Path::new(gff).exists() {
            return Err(format!("Error: GFF file does not exist: {}", gff).into());
        }
        if !gff.ends_with(".gff") && !gff.ends_with(".gff3") && !gff.ends_with(".gtf") {
            return Err(format!("Error: GFF file path must end with .gff, .gff3, or .gtf: {}", gff).into());
        }
    }

    // Validate SVG template file
    if let Some(svg_template) = &args.svg_template {
        if !Path::new(svg_template).exists() {
            return Err(format!("Error: SVG template file does not exist: {}", svg_template).into());
        }
        if !svg_template.ends_with(".svg") {
            return Err(format!("Error: SVG template file path must end with .svg: {}", svg_template).into());
        }
    }

    // Validate SVG reference sequence file
    if let Some(svg_ref) = &args.svg_ref {
        if !Path::new(svg_ref).exists() {
            return Err(format!("Error: SVG reference sequence file does not exist: {}", svg_ref).into());
        }
    }

    // Validate SVG signal type parameter
    if !matches!(
        args.svg_signal.as_str(),
        "all" | "stop" | "mutation" | "score"
    ) {
        return Err(format!(
            "Error: Invalid SVG signal type: {}. Supported types: all, stop, mutation, score",
            args.svg_signal
        )
        .into());
    }

    // Validate SVG strand type parameter
    if !matches!(args.svg_strand.as_str(), "+" | "-" | "both") {
        return Err(format!(
            "Error: Invalid SVG strand type: {}. Supported types: +, -, both",
            args.svg_strand
        )
        .into());
    }

    // Validate SVG base filter string
    if args.svg_bases.is_empty() {
        return Err("Error: SVG base filter string cannot be empty".into());
    }
    for base in args.svg_bases.chars() {
        if !"ACGT".contains(base.to_ascii_uppercase()) {
            return Err(format!("Error: SVG base filter string contains invalid base: {}", base).into());
        }
    }

    // Validate SVG maximum shift
    if args.svg_max_shift > 100 {
        return Err(format!(
            "Error: SVG maximum shift cannot exceed 100, current: {}",
            args.svg_max_shift
        )
        .into());
    }

    Ok(())
}

/// Ultra-efficient multi-threaded pileup processing using optimized I/O strategy
#[derive(Clone)]
struct WindowTask {
    region: String,
    start: u32,
    end: u32,
}

/// Function to merge Item
fn merge_item_into(target: &mut Item, src: &Item) {
    target.a += src.a;
    target.c += src.c;
    target.g += src.g;
    target.t += src.t;
    target.depth += src.depth;
    target.effective_depth += src.effective_depth;  // Merge quality-filtered depth
    target.ins += src.ins;
    target.del += src.del;
    target.rfc += src.rfc;
    target.pipc += src.pipc;
}

/// Extract cell label from BAM file name
/// Example: in_vivo_single_cell_Mut_transfected_RNAs_in_HEK293T_DMSO_RHX672.sort.bam -> RHX672
fn extract_cell_label(bam_path: &str) -> String {
    let path = Path::new(bam_path);
    if let Some(file_name) = path.file_stem() {
        if let Some(name_str) = file_name.to_str() {
            // Try to extract cell ID from file name (usually the part after the last underscore, or the part starting with RHX)
            // Match pattern: *RHX[0-9]+ or the part after the last underscore
            if let Some(rhx_pos) = name_str.find("RHX") {
                let after_rhx = &name_str[rhx_pos..];
                if let Some(end_pos) = after_rhx.find(|c: char| !c.is_ascii_digit() && c != 'X' && c != 'H' && c != 'R') {
                    return after_rhx[..end_pos].to_string();
                }
                return after_rhx.to_string();
            }
            // If RHX pattern is not found, use the part after the last underscore
            if let Some(last_underscore) = name_str.rfind('_') {
                return name_str[last_underscore + 1..].to_string();
            }
            // If neither, use the file name (without extension)
            return name_str.to_string();
        }
    }
    // If extraction fails, use file name hash
    format!("cell_{}", bam_path.len())
}

/// Calculate optimal CPU count for single-cell mode
/// Considerations:
/// - Number of BAM files (each file requires I/O)
/// - I/O bandwidth limitations
/// - Memory bandwidth limitations
/// - Parallel processing efficiency
fn calculate_optimal_cpus(bam_file_count: usize, current_cpus: usize) -> usize {
    // Base recommendation: 1-2 times the number of BAM files
    let base_recommendation = (bam_file_count * 2).min(64);
    
    // If there are many BAM files, consider I/O bottlenecks, recommend using more CPUs
    // But don't exceed 128, as there will be diminishing returns
    let recommended = if bam_file_count <= 10 {
        // Small number of files: use 2-4 times the file count
        (bam_file_count * 4).min(32)
    } else if bam_file_count <= 50 {
        // Medium number: use 1.5-2 times the file count
        (bam_file_count * 2).min(64)
    } else {
        // Large number of files: use 1-1.5 times the file count, but not exceeding 128
        (bam_file_count * 3 / 2).min(128)
    };
    
    // Ensure at least using current CPU count, but not exceeding system limits
    recommended.max(current_cpus).min(256)
}

/// Single-cell unified processing: read all BAM files at once, add cell labels, process uniformly, then split output by cell
fn process_single_cell_unified(
    bam_files: &[String],
    fasta_path: &str,
    output_dir: &Path,
    num_threads: usize,
    window: Option<usize>,
    strand_filter: &str,
    min_base_qual: u8,
    pcr_bias_correction: bool,
    weight_increase: f64,
    weight_decrease: f64,
    logger: &mut Option<&mut Logger>,
) -> Result<(), Box<dyn Error>> {
    if let Some(ref mut logger) = logger {
        logger.log("=== Starting single-cell unified processing ===")?;
        logger.log(&format!("Total BAM files: {}", bam_files.len()))?;
    }
    
    let start = Instant::now();
    
    // 1. Load reference sequences (keep in memory, do not release)
    if let Some(ref mut logger) = logger {
        logger.log("Loading reference sequences...")?;
    }
    let reference_sequences = load_reference_sequences(fasta_path)?;
    if let Some(ref mut logger) = logger {
        logger.log(&format!("Loaded {} reference sequences", reference_sequences.len()))?;
    }
    
    // 2. Single-cell mode: directly process all reference sequences, no data distribution scan needed
    if let Some(ref mut logger) = logger {
        logger.log("Single-cell mode: processing all reference sequences directly (no distribution scan)...")?;
        logger.log(&format!("Total reference sequences: {}", reference_sequences.len()))?;
    }
    
    // 3. Build window tasks (based on all reference sequences)
    let (num_chunks, progress_total, chunk_iter) = if let Some(wsize) = window {
        let mut windows: Vec<WindowTask> = Vec::new();
        let overlap_size = 3u32;
        
        // Directly process all regions of all reference sequences
        for (ref_name, ref_seq) in &reference_sequences {
            let ref_len = ref_seq.len() as u32;
            let step = wsize as u32;
            let mut s = 0u32;
            
            while s < ref_len {
                let e = (s + step).min(ref_len);
                windows.push(WindowTask {
                    region: ref_name.clone(),
                    start: s,
                    end: e,
                });
                if e >= ref_len {
                    break;
                }
                s = e.saturating_sub(overlap_size);
            }
        }
        
        // Shuffle windows to improve load balancing
        use std::collections::hash_map::DefaultHasher;
        use std::hash::{Hash, Hasher};
        let mut keyed: Vec<(u64, WindowTask)> = windows
            .into_iter()
            .map(|w| {
                let mut hasher = DefaultHasher::new();
                w.region.hash(&mut hasher);
                w.start.hash(&mut hasher);
                (hasher.finish(), w)
            })
            .collect();
        keyed.sort_by_key(|(k, _)| *k);
        let windows_shuffled: Vec<WindowTask> = keyed.into_iter().map(|(_, w)| w).collect();
        
        let target_multiplicity: usize = 32;
        let mut batch_size = (windows_shuffled.len() / (num_threads.saturating_mul(target_multiplicity))).max(1);
        if batch_size > 4096 {
            batch_size = 4096;
        }
        
        let chunks: Vec<Vec<WindowTask>> = windows_shuffled
            .chunks(batch_size)
            .map(|c| c.to_vec())
            .collect();
        
        (chunks.len(), chunks.iter().map(|c| c.len()).sum(), chunks)
    } else {
        // No window mode: directly process all reference sequences
        let mut windows: Vec<WindowTask> = Vec::new();
        for (ref_name, ref_seq) in &reference_sequences {
            let ref_len = ref_seq.len() as u32;
            // Directly process the entire reference sequence
            windows.push(WindowTask {
                region: ref_name.clone(),
                start: 0,
                end: ref_len,
            });
        }
        
        let chunk_size = (windows.len() + num_threads - 1) / num_threads.max(1);
        let chunks: Vec<Vec<WindowTask>> = windows
            .chunks(chunk_size.max(1))
            .map(|chunk| chunk.to_vec())
            .collect();
        (chunks.len(), windows.len(), chunks)
    };
    
    if let Some(ref mut logger) = logger {
        logger.log(&format!("Created {} chunks with {} total windows", num_chunks, progress_total))?;
    }
    
    if num_chunks == 0 {
        if let Some(ref mut logger) = logger {
            logger.log("Error: No processing chunks created!")?;
            logger.log(&format!("  Reference sequence count: {}", reference_sequences.len()))?;
            logger.log(&format!("  Window size: {:?}", window))?;
        }
        return Err("No processing chunks created, please check reference sequence loading results".into());
    }
    
    // 4. Unified processing of all BAM files (parallel processing of windows, each window processes all files)
    use std::sync::atomic::{AtomicUsize, Ordering};
    use std::sync::Arc;
    
    let completed = Arc::new(AtomicUsize::new(0));
    
    // Store statistical results for each cell: HashMap<cell_id, HashMap<region, HashMap<position, Item>>>
    let mut cell_results: HashMap<String, (HashMap<String, HashMap<u32, Item>>, HashMap<String, HashMap<u32, Item>>, HashMap<String, Vec<u32>>)> = HashMap::new();
    
    // Initialize result structure for each cell
    for bam_file in bam_files {
        let cell_id = extract_cell_label(bam_file);
        cell_results.insert(cell_id, (HashMap::new(), HashMap::new(), HashMap::new()));
    }
    
    if let Some(ref mut logger) = logger {
        logger.log("Processing windows with unified reads...")?;
        logger.log(&format!("Total windows: {}, Total chunks: {}, BAM files: {}", progress_total, num_chunks, bam_files.len()))?;
        
        // CPU usage recommendation
        let recommended_cpus = calculate_optimal_cpus(bam_files.len(), num_threads);
        if num_threads < recommended_cpus {
            logger.log(&format!("Tip: Detected {} BAM files, recommend using {} CPUs for better performance (current: {})", 
                bam_files.len(), recommended_cpus, num_threads))?;
        }
    }
    
    // Use atomic counter to track progress (these types are already imported at the top of the file)
    let completed = Arc::new(AtomicUsize::new(0));
    let completed_clone = completed.clone();
    let start_time = Arc::new(std::sync::Mutex::new(Instant::now()));
    
    // Process each window in parallel
    let chunk_results: Vec<Result<HashMap<String, (HashMap<String, HashMap<u32, Item>>, HashMap<String, HashMap<u32, Item>>, HashMap<String, Vec<u32>>)>, Box<dyn Error + Send + Sync>>> = chunk_iter
        .into_par_iter()
        .map(|window_chunk| {
            let result = process_region_chunk_unified(
                bam_files,
                &reference_sequences,
                window_chunk,
                strand_filter,
            );
            
            // Update progress
            let completed_count = completed_clone.fetch_add(1, Ordering::Relaxed) + 1;
            
            // Update progress every 10 chunks or every 5 seconds
            if completed_count % 10 == 0 || completed_count == num_chunks {
                let elapsed = {
                    let start = start_time.lock().unwrap();
                    start.elapsed()
                };
                let percentage = (completed_count as f64 * 100.0) / num_chunks as f64;
                let elapsed_secs = elapsed.as_secs_f64();
                let rate = if elapsed_secs > 0.0 {
                    completed_count as f64 / elapsed_secs
                } else {
                    0.0
                };
                let remaining = if rate > 0.0 {
                    (num_chunks - completed_count) as f64 / rate
                } else {
                    0.0
                };
                
                print!("\r[Running] Processing {}/{} chunks ({:.1}%) | Windows: {} | Speed: {:.1} chunks/s | ETA: {:.0}s", 
                    completed_count, num_chunks, percentage, progress_total, rate, remaining);
                std::io::Write::flush(&mut std::io::stdout()).ok();
            }
            
            result
        })
        .collect();
    
    let total_elapsed = start_time.lock().unwrap().elapsed();
    println!("\r[Done] Processed {} chunks ({} windows) in {:.1}s", 
        num_chunks, progress_total, total_elapsed.as_secs_f64());
    
    // 4.5. Use intermediate files to avoid memory accumulation
    // Create temporary files for each cell, write immediately after processing chunk
    if let Some(ref mut logger) = logger {
        logger.log("Using intermediate files to avoid memory accumulation...")?;
    }
    
    // Create temporary directory to store intermediate files
    let temp_dir = output_dir.join(".tmp_intermediate");
    fs::create_dir_all(&temp_dir)?;
    
    // Create temporary file writer for each cell
    let mut cell_temp_writers: HashMap<String, std::io::BufWriter<std::fs::File>> = HashMap::new();
    
    for bam_file in bam_files {
        let cell_id = extract_cell_label(bam_file);
        if !cell_temp_writers.contains_key(&cell_id) {
            let temp_file_path = temp_dir.join(format!("{}.tmp.csv", cell_id));
            let temp_file = std::fs::File::create(&temp_file_path)?;
            let writer = std::io::BufWriter::with_capacity(2 * 1024 * 1024, temp_file); // 2MB buffer
            cell_temp_writers.insert(cell_id, writer);
        }
    }
    
    // Process each chunk's results, write to intermediate files immediately
    let mut merge_count = 0;
    for chunk_result in chunk_results {
        let cell_chunk_results = chunk_result.map_err(|e| format!("{}", e))?;
        
        // Collect current chunk's data for each cell and write immediately
        for (cell_id, (pos_items, neg_items, _positions)) in cell_chunk_results {
            if let Some(writer) = cell_temp_writers.get_mut(&cell_id) {
                use std::io::Write;
                
                // Write positive strand data
                for (region, items) in pos_items {
                    for (pos, item) in items {
                        if item.depth > 0 {
                            writeln!(writer, "{},{},{},{},{},{},{},{},{},{},{},{},{}",
                                region, "+", pos + 1, item.refbase,
                                item.rfc, item.pipc, item.depth,
                                item.ins, item.del,
                                item.a, item.c, item.g, item.t)?;
                        }
                    }
                }
                
                // Write negative strand data
                for (region, items) in neg_items {
                    for (pos, item) in items {
                        if item.depth > 0 {
                            writeln!(writer, "{},{},{},{},{},{},{},{},{},{},{},{},{}",
                                region, "-", pos + 1, item.refbase,
                                item.rfc, item.pipc, item.depth,
                                item.ins, item.del,
                                item.a, item.c, item.g, item.t)?;
                        }
                    }
                }
            }
        }
        
        merge_count += 1;
        if merge_count % 100 == 0 {
            // Periodically flush all writers
            for writer in cell_temp_writers.values_mut() {
                writer.flush()?;
            }
            if let Some(ref mut logger) = logger {
                logger.log(&format!("Flushed intermediate files after {} chunks", merge_count))?;
            }
        }
    }
    
    // Flush and close all writers
    for writer in cell_temp_writers.values_mut() {
        writer.flush()?;
    }
    drop(cell_temp_writers); // Ensure all files are closed
    
    if let Some(ref mut logger) = logger {
        logger.log("Intermediate files written, merging and deduplicating to final outputs...")?;
    }
    
    // 5. Read intermediate files, merge and deduplicate by position, output final CSV
    // Read intermediate files for each cell, merge Items by position, then output
    let collect_pos = strand_filter != "-";
    let collect_neg = strand_filter != "+";
    
    for bam_file in bam_files {
        let cell_id = extract_cell_label(bam_file);
        let temp_file_path = temp_dir.join(format!("{}.tmp.csv", cell_id));
        
        if !temp_file_path.exists() {
            continue; // Skip cells with no data
        }
        
        if let Some(ref mut logger) = logger {
            logger.log(&format!("Merging and writing output for cell: {}", cell_id))?;
        }
        
        // Read intermediate file, merge Items by position
        let mut cell_results: HashMap<String, HashMap<u32, (Item, Item)>> = HashMap::new(); // (region, pos) -> (pos_item, neg_item)
        
        let temp_file = std::fs::File::open(&temp_file_path)?;
        let reader = std::io::BufReader::new(temp_file);
        use std::io::BufRead;
        
        for line in reader.lines() {
            let line = line?;
            if line.starts_with("ChrID") {
                continue; // Skip header (if present)
            }
            
            let parts: Vec<&str> = line.split(',').collect();
            if parts.len() < 13 {
                continue;
            }
            
            let region = parts[0].to_string();
            let strand = parts[1];
            let pos_str = parts[2].trim();
            let pos = pos_str.parse::<u32>().unwrap_or(0);
            if pos == 0 {
                continue;
            }
            let pos = pos - 1; // Convert back to 0-based
            
            let ref_base = parts[3].chars().next().unwrap_or('N');
            let rfc = parts[4].parse::<usize>().unwrap_or(0);
            let pipc = parts[5].parse::<usize>().unwrap_or(0);
            let depth = parts[6].parse::<usize>().unwrap_or(0);
            // Read effective_depth if available (column 7), otherwise use depth
            let effective_depth = if parts.len() > 13 {
                parts[7].parse::<usize>().unwrap_or(depth)
            } else {
                depth  // Fallback to depth if effective_depth column not present
            };
            let ins = if parts.len() > 13 { parts[8].parse::<usize>().unwrap_or(0) } else { parts[7].parse::<usize>().unwrap_or(0) };
            let del = if parts.len() > 13 { parts[9].parse::<usize>().unwrap_or(0) } else { parts[8].parse::<usize>().unwrap_or(0) };
            let a = if parts.len() > 13 { parts[10].parse::<u32>().unwrap_or(0) } else { parts[9].parse::<u32>().unwrap_or(0) };
            let c = if parts.len() > 13 { parts[11].parse::<u32>().unwrap_or(0) } else { parts[10].parse::<u32>().unwrap_or(0) };
            let g = if parts.len() > 13 { parts[12].parse::<u32>().unwrap_or(0) } else { parts[11].parse::<u32>().unwrap_or(0) };
            let t = if parts.len() > 13 { parts[13].parse::<u32>().unwrap_or(0) } else { parts[12].parse::<u32>().unwrap_or(0) };
            
            let item = Item {
                tid: region.clone(),
                strand: if strand == "+" { '+' } else { '-' },
                position: pos,
                refbase: ref_base,
                rfc,
                pipc,
                depth,
                effective_depth,
                ins,
                del,
                a,
                c,
                g,
                t,
            };
            
            let region_map = cell_results.entry(region.clone()).or_insert_with(HashMap::new);
            let (pos_item, neg_item) = region_map.entry(pos).or_insert_with(|| {
                let (pi, ni) = itemcreate(region.clone(), pos, ref_base);
                (pi, ni)
            });
            
            if strand == "+" && collect_pos {
                merge_item_into(pos_item, &item);
            } else if strand == "-" && collect_neg {
                merge_item_into(neg_item, &item);
            }
        }
        
        // Output final CSV file
        let output_file = output_dir.join(format!("{}.csv", cell_id));
        let output_str = output_file.to_str().unwrap();
        let output_file = std::fs::File::create(output_str)?;
        let mut output = std::io::BufWriter::with_capacity(2 * 1024 * 1024, output_file);
        
        use std::io::Write;
        // Write header with or without effective_depth column
        // effective_depth is output if quality filtering is enabled OR PCR bias correction is enabled
        if pcr_bias_correction || min_base_qual > 0 {
            writeln!(output, "ChrID,pipe_truncation_Strand,pipe_truncation_ChrPos,rf_mutation_Base,rf_mutation_Count,pipe_truncation_count,depth,effective_depth,rf_mutation_ins,rf_mutation_del,base_A,base_C,base_G,base_T")?;
        } else {
            writeln!(output, "ChrID,pipe_truncation_Strand,pipe_truncation_ChrPos,rf_mutation_Base,rf_mutation_Count,pipe_truncation_count,depth,rf_mutation_ins,rf_mutation_del,base_A,base_C,base_G,base_T")?;
        }
        
        // Output in reference sequence order
        for (ref_name, ref_seq) in &reference_sequences {
            if let Some(region_map) = cell_results.get(ref_name) {
                // Collect all positions and sort
                let mut positions: Vec<u32> = region_map.keys().cloned().collect();
                positions.sort();
                
                for pos in positions {
                    let (pos_item, neg_item) = &region_map[&pos];
                    
                    if collect_pos && pos_item.depth > 0 {
                        if pcr_bias_correction || min_base_qual > 0 {
                            let effective_depth = if pcr_bias_correction {
                                // Apply PCR bias correction to quality-filtered depth
                                let base_depth = if min_base_qual > 0 {
                                    pos_item.effective_depth
                                } else {
                                    pos_item.depth
                                };
                                calculate_effective_depth(
                                    base_depth, 
                                    pos_item.rfc,
                                    weight_increase,
                                    weight_decrease,
                                )
                            } else {
                                // Only quality filtering, no PCR bias correction
                                pos_item.effective_depth
                            };
                            writeln!(output, "{},{},{},{},{},{},{},{},{},{},{},{},{},{}",
                                ref_name, pos_item.strand, pos + 1, pos_item.refbase,
                                pos_item.rfc, pos_item.pipc, pos_item.depth, effective_depth,
                                pos_item.ins, pos_item.del,
                                pos_item.a, pos_item.c, pos_item.g, pos_item.t)?;
                        } else {
                            writeln!(output, "{},{},{},{},{},{},{},{},{},{},{},{},{}",
                                ref_name, pos_item.strand, pos + 1, pos_item.refbase,
                                pos_item.rfc, pos_item.pipc, pos_item.depth,
                                pos_item.ins, pos_item.del,
                                pos_item.a, pos_item.c, pos_item.g, pos_item.t)?;
                        }
                    }
                    
                    if collect_neg && neg_item.depth > 0 {
                        if pcr_bias_correction || min_base_qual > 0 {
                            let effective_depth = if pcr_bias_correction {
                                // Apply PCR bias correction to quality-filtered depth
                                let base_depth = if min_base_qual > 0 {
                                    neg_item.effective_depth
                                } else {
                                    neg_item.depth
                                };
                                calculate_effective_depth(
                                    base_depth, 
                                    neg_item.rfc,
                                    weight_increase,
                                    weight_decrease,
                                )
                            } else {
                                // Only quality filtering, no PCR bias correction
                                neg_item.effective_depth
                            };
                            writeln!(output, "{},{},{},{},{},{},{},{},{},{},{},{},{},{}",
                                ref_name, neg_item.strand, pos + 1, neg_item.refbase,
                                neg_item.rfc, neg_item.pipc, neg_item.depth, effective_depth,
                                neg_item.ins, neg_item.del,
                                neg_item.a, neg_item.c, neg_item.g, neg_item.t)?;
                        } else {
                            writeln!(output, "{},{},{},{},{},{},{},{},{},{},{},{},{}",
                                ref_name, neg_item.strand, pos + 1, neg_item.refbase,
                                neg_item.rfc, neg_item.pipc, neg_item.depth,
                                neg_item.ins, neg_item.del,
                                neg_item.a, neg_item.c, neg_item.g, neg_item.t)?;
                        }
                    }
                }
            }
        }
        
        output.flush()?;
        
        // Delete temporary file
        fs::remove_file(&temp_file_path)?;
    }
    
    // Delete temporary directory
    fs::remove_dir(&temp_dir)?;
    
    let duration = start.elapsed();
    if let Some(ref mut logger) = logger {
        logger.log(&format!("Single-cell unified processing completed in {:?}", duration))?;
    }
    
    Ok(())
}

/// Unified processing of region chunks: read all BAM file reads at once, process uniformly, group by cell label
fn process_region_chunk_unified(
    bam_files: &[String],
    reference_sequences: &HashMap<String, Vec<char>>,
    windows: Vec<WindowTask>,
    strand_filter: &str,
) -> Result<HashMap<String, (HashMap<String, HashMap<u32, Item>>, HashMap<String, HashMap<u32, Item>>, HashMap<String, Vec<u32>>)>, Box<dyn Error + Send + Sync>> {
    use rust_htslib::bam::{IndexedReader, Read};
    
    // Store results for each cell
    let mut cell_results: HashMap<String, (HashMap<String, HashMap<u32, Item>>, HashMap<String, HashMap<u32, Item>>, HashMap<String, Vec<u32>>)> = HashMap::new();
    
    // Initialize result structure for each cell
    for bam_file in bam_files {
        let cell_id = extract_cell_label(bam_file);
        cell_results.insert(cell_id, (HashMap::new(), HashMap::new(), HashMap::new()));
    }
    
    let strand_mode = normalize_strand_filter(strand_filter);
    let collect_pos = strand_mode != "-";
    let collect_neg = strand_mode != "+";
    
    // Process each window
    for w in &windows {
        if !reference_sequences.contains_key(&w.region) {
            continue;
        }
        
        let ref_seq = &reference_sequences[&w.region];
        let ref_len = ref_seq.len() as u32;
        
        // Unified processing: collect all reads for this window from all BAM files at once
        // First collect all reads, then process uniformly by position
        let mut all_reads: Vec<(String, rust_htslib::bam::Record)> = Vec::new();
        
        // Create reader for each BAM file and extract cell label, collect reads
        for bam_file in bam_files {
            let cell_id = extract_cell_label(bam_file);
            let mut bam = IndexedReader::from_path(bam_file)?;
            bam.set_threads(1);
            
            // Query this window, collect all reads
            let start1 = (w.start + 1) as u32;
            let end1 = w.end as u32;
            let region_query = format!("{}:{}-{}", w.region, start1, end1);
            
            if bam.fetch(&region_query).is_ok() {
                let mut record = rust_htslib::bam::Record::new();
                loop {
                    match bam.read(&mut record) {
                        Some(Ok(())) => {
                            let read_start = record.pos() as u32;
                            if read_start < ref_len && read_start + record.seq().len() as u32 >= w.start {
                                // Collect all reads covering this window, mark with cell label
                                all_reads.push((cell_id.clone(), record.clone()));
                            }
                        }
                        Some(Err(_)) | None => break,
                    }
                }
            }
        }
        
        // Unified processing: iterate by position, process all reads covering each position uniformly
        let mut processed_positions = std::collections::HashSet::new();
        
        let end_pos = w.end.min(ref_len);
        for pos in w.start..end_pos {
            if processed_positions.contains(&pos) {
                continue;
            }
            processed_positions.insert(pos);
            if pos >= ref_len || pos < w.start || pos > w.end {
                continue;
            }
            
            let ref_base = ref_seq[pos as usize];
            
            // Initialize Item for this position for each cell
            let mut cell_pos_items: HashMap<String, Item> = HashMap::new();
            let mut cell_neg_items: HashMap<String, Item> = HashMap::new();
            
            // Process all reads covering this position uniformly
            for (cell_id, record) in &all_reads {
                let read_start = record.pos() as u32;
                let read_end = read_start + record.seq().len() as u32;
                
                // Check if this read covers current position
                if pos < read_start || pos >= read_end {
                    continue;
                }
                
                // Initialize Item for this cell (if not already)
                if !cell_pos_items.contains_key(cell_id) && collect_pos {
                    let (pi, _) = itemcreate(w.region.clone(), pos, ref_base);
                    cell_pos_items.insert(cell_id.clone(), pi);
                }
                if !cell_neg_items.contains_key(cell_id) && collect_neg {
                    let (_, ni) = itemcreate(w.region.clone(), pos, ref_base);
                    cell_neg_items.insert(cell_id.clone(), ni);
                }
            }
            
            // Process all reads covering this position uniformly
            for (cell_id, record) in &all_reads {
                let read_start = record.pos() as u32;
                let read_end = read_start + record.seq().len() as u32;
                
                // Check if this read covers current position
                if pos < read_start || pos >= read_end {
                    continue;
                }
                let is_reverse = record.is_reverse();
                
                // Process stop signal (read start position)
                // Theoretical correction: positive strand start_site - 1, negative strand start_site + 1
                if read_start == pos {
                    let stop_signal_pos = if is_reverse {
                        // Negative strand: correct to start_site + 1
                        pos + 1
                    } else {
                        // Positive strand: correct to start_site - 1
                        if pos > 0 { pos - 1 } else { 0 }
                    };
                    
                    // Check if corrected position is within valid range
                    if stop_signal_pos < ref_len && stop_signal_pos >= w.start {
                        let ref_base_stop = ref_seq[stop_signal_pos as usize];
                        
                        if is_reverse && collect_neg {
                            if let Some(neg_item) = cell_neg_items.get_mut(cell_id) {
                                neg_item.pipc += 1;
                            } else {
                                let (_, ni) = itemcreate(w.region.clone(), stop_signal_pos, ref_base_stop);
                                cell_neg_items.insert(cell_id.clone(), ni);
                                cell_neg_items.get_mut(cell_id).unwrap().pipc += 1;
                            }
                            processed_positions.insert(stop_signal_pos);
                        } else if !is_reverse && collect_pos {
                            if let Some(pos_item) = cell_pos_items.get_mut(cell_id) {
                                pos_item.pipc += 1;
                            } else {
                                let (pi, _) = itemcreate(w.region.clone(), stop_signal_pos, ref_base_stop);
                                cell_pos_items.insert(cell_id.clone(), pi);
                                cell_pos_items.get_mut(cell_id).unwrap().pipc += 1;
                            }
                            processed_positions.insert(stop_signal_pos);
                        }
                    }
                }
                
                // Process mutation signal at current position
                let qpos = (pos - read_start) as usize;
                if qpos < record.seq().len() {
                    let base = findchar(qpos, record.seq());
                    let ref_numb = chartonum(ref_base);
                    
                    let update_item = |item: &mut Item| {
                        match base {
                            1 => item.a += 1,
                            2 => item.c += 1,
                            4 => item.g += 1,
                            8 => item.t += 1,
                            _ => {}
                        }
                        item.depth += 1;
                        if ref_numb != base {
                            item.rfc += 1;
                        }
                    };
                    
                    if is_reverse && collect_neg {
                        if let Some(neg_item) = cell_neg_items.get_mut(cell_id) {
                            update_item(neg_item);
                        }
                    } else if !is_reverse && collect_pos {
                        if let Some(pos_item) = cell_pos_items.get_mut(cell_id) {
                            update_item(pos_item);
                        }
                    }
                }
            }
            
            // Merge results into cell results
            for (cell_id, pos_item) in cell_pos_items {
                let (pos_items_map, _, region_positions) = cell_results.get_mut(&cell_id).unwrap();
                let pos_map = pos_items_map.entry(w.region.clone()).or_insert_with(HashMap::new);
                if let Some(existing) = pos_map.get_mut(&pos) {
                    merge_item_into(existing, &pos_item);
                } else {
                    pos_map.insert(pos, pos_item);
                }
                let pos_vec = region_positions.entry(w.region.clone()).or_insert_with(Vec::new);
                if !pos_vec.contains(&pos) {
                    pos_vec.push(pos);
                }
            }
            
            for (cell_id, neg_item) in cell_neg_items {
                let (_, neg_items_map, region_positions) = cell_results.get_mut(&cell_id).unwrap();
                let neg_map = neg_items_map.entry(w.region.clone()).or_insert_with(HashMap::new);
                if let Some(existing) = neg_map.get_mut(&pos) {
                    merge_item_into(existing, &neg_item);
                } else {
                    neg_map.insert(pos, neg_item);
                }
                let pos_vec = region_positions.entry(w.region.clone()).or_insert_with(Vec::new);
                if !pos_vec.contains(&pos) {
                    pos_vec.push(pos);
                }
            }
            
        }
    }
    
    Ok(cell_results)
}

fn process_sample_pileup_ultra_fast(
    bam_path: &str,
    fasta_path: &str,
    output_path: &str,
    num_threads: usize,
    window: Option<usize>,
    strand_filter: &str,
    min_base_qual: u8,
    pcr_bias_correction: bool,
    weight_increase: f64,
    weight_decrease: f64,
    logger: &mut Option<&mut Logger>,
) -> Result<(), Box<dyn Error>> {
    if !Path::new(bam_path).exists() {
        return Err(format!("BAM file does not exist: {}", bam_path).into());
    }
    if !Path::new(fasta_path).exists() {
        return Err(format!("FASTA file does not exist: {}", fasta_path).into());
    }

    if let Some(w) = window {
        println!("[loading data] Threads={}, Window={}bp", num_threads, w);
    } else {
        println!("[loading data] Threads={} (no window)", num_threads);
    }
    println!("    BAM={}", bam_path);
    println!("    FASTA={}", fasta_path);
    println!("    Strand Filter={}", strand_filter);
    if min_base_qual > 0 {
        println!("    Base Quality Filter={} (Phred score)", min_base_qual);
    } else {
        println!("    Base Quality Filter=disabled");
    }
    println!();

    let start = Instant::now();

    // Load reference sequences (keep in memory, do not release)
    if let Some(ref mut logger) = logger {
        logger.log("Loading reference sequences...")?;
    }
    let reference_sequences = load_reference_sequences(fasta_path)?;
    if let Some(ref mut logger) = logger {
        logger.log(&format!("Loaded {} reference sequences", reference_sequences.len()))?;
    }
    
    // For single-cell data: first scan data distribution
    if let Some(ref mut logger) = logger {
        logger.log("Scanning BAM data distribution (single-cell optimization)...")?;
    }
    let distribution = scan_bam_distribution(bam_path, &reference_sequences, num_threads, 100000)?;
    
    // Count regions with data
    let mut total_intervals = 0;
    let mut total_covered_bases = 0u64;
    for (ref_name, intervals) in &distribution {
        total_intervals += intervals.len();
        for &(start, end) in intervals {
            total_covered_bases += (end - start) as u64;
        }
    }
    
    if let Some(ref mut logger) = logger {
        logger.log(&format!(
            "Data distribution scan completed: {} intervals covering {} bases",
            total_intervals, total_covered_bases
        ))?;
    }

    // Set up rayon thread pool (ignore error if already initialized)
    // When batch processing in single-cell mode, thread pool may already be initialized
    rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build_global()
        .ok(); // Use ok() instead of unwrap() to avoid duplicate initialization error

    // Build tasks: windowed or regionalized
    // For single-cell data: optimize window division based on data distribution, only process regions with data
    let (num_chunks, progress_total, chunk_iter): (usize, usize, Vec<Vec<WindowTask>>) =
        if let Some(wsize) = window {
            let mut windows: Vec<WindowTask> = Vec::new();
            let overlap_size = 3u32; // 3bp overlapping windows to ensure boundary sites are not lost

            // Build windows based on data distribution: only process regions with data
            for (ref_name, ref_seq) in &reference_sequences {
                let ref_len = ref_seq.len() as u32;
                
                // Get data distribution intervals for this region
                let intervals = distribution.get(ref_name).cloned().unwrap_or_default();
                
                if intervals.is_empty() {
                    // If no data distribution information, use traditional method: entire region
                    let step = wsize as u32;
                    if step >= ref_len {
                        windows.push(WindowTask {
                            region: ref_name.clone(),
                            start: 0,
                            end: ref_len,
                        });
                    } else {
                        let mut s: u32 = 0;
                        while s < ref_len {
                            let e = (s + step).min(ref_len);
                            windows.push(WindowTask {
                                region: ref_name.clone(),
                                start: s,
                                end: e,
                            });
                            if e >= ref_len {
                                break;
                            }
                            s = e.saturating_sub(overlap_size);
                        }
                    }
                } else {
                    // Based on data distribution: only process intervals with data
                    for (interval_start, interval_end) in intervals {
                        let step = wsize as u32;
                        let mut s = interval_start;
                        let end = interval_end.min(ref_len);
                        
                        while s < end {
                            let e = (s + step).min(end);
                            windows.push(WindowTask {
                                region: ref_name.clone(),
                                start: s,
                                end: e,
                            });
                            
                            if e >= end {
                                break;
                            }
                            s = e.saturating_sub(overlap_size);
                        }
                    }
                }
            }

            // Approximate coverage-weighted load balancing:
            // 1) Hash shuffle different chromosomes/positions to reduce heavy region clustering
            // 2) Create many more small batch tasks than thread count for dynamic parallelization
            use std::collections::hash_map::DefaultHasher;
            use std::hash::{Hash, Hasher};

            let mut keyed: Vec<(u64, WindowTask)> = windows
                .into_iter()
                .map(|w| {
                    let mut hasher = DefaultHasher::new();
                    // Use region and start as key for hashing to achieve cross-chromosome interleaving and shuffling
                    w.region.hash(&mut hasher);
                    w.start.hash(&mut hasher);
                    let key = hasher.finish();
                    (key, w)
                })
                .collect();

            // Sort by hash key to produce stable shuffle order
            keyed.sort_by_key(|(k, _)| *k);
            let windows_shuffled: Vec<WindowTask> = keyed.into_iter().map(|(_, w)| w).collect();

            // Calculate small batch task size: target is 20-50 times the number of threads
            let target_multiplicity: usize = 32; // Empirical value, can be exposed as parameter later
            let mut batch_size =
                (windows_shuffled.len() / (num_threads.saturating_mul(target_multiplicity))).max(1);
            // Avoid batches being too small causing excessive scheduling/merging overhead
            if batch_size > 4096 {
                batch_size = 4096;
            }

            let chunks: Vec<Vec<WindowTask>> = windows_shuffled
                .chunks(batch_size)
                .map(|c| c.to_vec())
                .collect();

            (chunks.len(), chunks.iter().map(|c| c.len()).sum(), chunks)
        } else {
            // No window: based on data distribution, one window per interval with data
            let mut windows: Vec<WindowTask> = Vec::new();
            for (ref_name, ref_seq) in &reference_sequences {
                let ref_len = ref_seq.len() as u32;
                
                // Get data distribution intervals for this region
                let intervals = distribution.get(ref_name).cloned().unwrap_or_default();
                
                if intervals.is_empty() {
                    // If no data distribution information, use entire region
                    windows.push(WindowTask {
                        region: ref_name.clone(),
                        start: 0,
                        end: ref_len,
                    });
                } else {
                    // Only process intervals with data
                    for (interval_start, interval_end) in intervals {
                        windows.push(WindowTask {
                            region: ref_name.clone(),
                            start: interval_start,
                            end: interval_end.min(ref_len),
                        });
                    }
                }
            }
            let chunk_size = (windows.len() + num_threads - 1) / num_threads.max(1);
            let chunks: Vec<Vec<WindowTask>> = windows
                .chunks(chunk_size.max(1))
                .map(|chunk| chunk.to_vec())
                .collect();
            (chunks.len(), windows.len(), chunks)
        };

    // Use atomic counter to track progress
    use std::sync::atomic::{AtomicUsize, Ordering};
    use std::sync::Arc;

    let completed = Arc::new(AtomicUsize::new(0));

    // Use dynamic progress display
    print!(
        "\r[Running] Dividing {} tasks into {} chunks",
        progress_total, num_chunks
    );
    std::io::stdout().flush()?;

    // Process each region chunk in parallel
    let chunk_results: Vec<
        Result<
            (
                HashMap<String, HashMap<u32, Item>>,
                HashMap<String, HashMap<u32, Item>>,
                HashMap<String, Vec<u32>>,
            ),
            Box<dyn Error + Send + Sync>,
        >,
    > = chunk_iter
        .into_par_iter()
        .map(|window_chunk| {
            let result = process_region_chunk_ultra_fast(
                bam_path,
                &reference_sequences,
                window_chunk,
                strand_filter,
                min_base_qual,
            );

            // Update progress
            let completed_count = completed.fetch_add(1, Ordering::Relaxed) + 1;
            let percentage = (completed_count as f64 * 100.0) / num_chunks as f64;
            print!(
                "\r[Running] Processing {}/{} chunks ({:.1}%)",
                completed_count, num_chunks, percentage
            );
            std::io::stdout().flush().ok();

            result
        })
        .collect();

    // Overwrite display to show completion status
    println!("\r[Done] Processing {} chunks", num_chunks);
    let mut pos_items_map: HashMap<String, HashMap<u32, Item>> = HashMap::new();
    let mut neg_items_map: HashMap<String, HashMap<u32, Item>> = HashMap::new();
    let mut region_positions: HashMap<String, Vec<u32>> = HashMap::new();

    for (chunk_idx, chunk_result) in chunk_results.into_iter().enumerate() {
        let (pos_items, neg_items, positions) = match chunk_result {
            Ok(result) => result,
            Err(e) => return Err(format!("Parallel processing error (chunk {}): {}", chunk_idx, e).into()),
        };

        // Merge positive strand data (accumulate)
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

        // Merge negative strand data (accumulate)
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

        // Merge position information
        for (region, mut positions_vec) in positions {
            let target_vec = region_positions.entry(region).or_insert_with(Vec::new);
            target_vec.append(&mut positions_vec);
        }
    }

    // Sort and deduplicate position information (handle duplicate sites from overlapping windows)
    for (_, positions_vec) in region_positions.iter_mut() {
        positions_vec.sort();
        positions_vec.dedup();
    }

    // Use buffered writing to optimize output performance
    let total_regions = reference_sequences.len();
    let collect_pos = strand_filter != "-";
    let collect_neg = strand_filter != "+";
    let mut write_progress = crate::progress::SimpleProgress::new(total_regions);

    // For PCR bias correction with global fitting: collect all data first
    let mut all_items: Vec<(String, char, u32, Item)> = Vec::new();
    
    if pcr_bias_correction {
        print!("\r[Collecting] Collecting data for global fitting...");
        std::io::stdout().flush()?;
        
        for (ref_name, ref_seq) in &reference_sequences {
            let empty_vec = Vec::new();
            let actual_positions = region_positions.get(ref_name).unwrap_or(&empty_vec);
            
            for &pos in actual_positions {
                let ref_base = if pos < ref_seq.len() as u32 {
                    ref_seq[pos as usize]
                } else {
                    'N'
                };
                
                if collect_pos {
                    if let Some(items) = pos_items_map.get(ref_name) {
                        if let Some(item) = items.get(&pos) {
                            if item.depth > 0 {
                                all_items.push((ref_name.clone(), '+', pos, item.clone()));
                            }
                        }
                    }
                }
                
                if collect_neg {
                    if let Some(items) = neg_items_map.get(ref_name) {
                        if let Some(item) = items.get(&pos) {
                            if item.depth > 0 {
                                all_items.push((ref_name.clone(), '-', pos, item.clone()));
                            }
                        }
                    }
                }
            }
        }
        
        println!("\r[Fitting] Fitting Chi-Square distribution to {} positions...", all_items.len());
    }

    print!("\r[Writing] Writing result file...");
    std::io::stdout().flush()?;
    let output_file = std::fs::File::create(output_path)?;
    let mut output = std::io::BufWriter::with_capacity(1024 * 1024, output_file); // 1MB buffer

    use std::io::Write;
    // Write header with or without effective_depth column
    // effective_depth is output if quality filtering is enabled OR PCR bias correction is enabled
    if pcr_bias_correction || min_base_qual > 0 {
        writeln!(output, "ChrID,pipe_truncation_Strand,pipe_truncation_ChrPos,rf_mutation_Base,rf_mutation_Count,pipe_truncation_count,depth,effective_depth,rf_mutation_ins,rf_mutation_del,base_A,base_C,base_G,base_T")?;
    } else {
        writeln!(output, "ChrID,pipe_truncation_Strand,pipe_truncation_ChrPos,rf_mutation_Base,rf_mutation_Count,pipe_truncation_count,depth,rf_mutation_ins,rf_mutation_del,base_A,base_C,base_G,base_T")?;
    }

    // Batch write data for all regions
    let mut write_buffer = String::with_capacity(1024 * 1024); // 1MB write buffer

    // If PCR bias correction is enabled, perform global fitting first
    let effective_depth_map: HashMap<usize, usize> = if pcr_bias_correction && !all_items.is_empty() {
        // Collect depth and mutation_count data
        // Use quality-filtered effective_depth if quality filtering is enabled
        let mut depths = Vec::new();
        let mut mutation_rates = Vec::new();
        
        for (_, _, _, item) in &all_items {
            // Use quality-filtered depth for fitting if quality filtering is enabled
            let base_depth = if min_base_qual > 0 {
                item.effective_depth
            } else {
                item.depth
            };
            
            if base_depth > 0 {
                let mutation_rate = item.rfc as f64 / base_depth as f64;
                if mutation_rate > 0.0 && mutation_rate <= 0.1 {
                    depths.push(base_depth as f64);
                    mutation_rates.push(mutation_rate);
                }
            }
        }
        
        if depths.len() >= 10 {
            // Fit Chi-Square distribution
            if let Some((a, scale, r2)) = fit_chi_square_to_data(&depths, &mutation_rates) {
                if let Some(ref mut logger) = logger {
                    logger.log(&format!("Chi-Square fit: a={:.6e}, scale={:.2}, R²={:.4}", a, scale, r2))?;
                }
                
                // Calculate minimum area line (target_y)
                let max_depth = depths.iter().fold(0.0_f64, |a, &b| a.max(b));
                let x_fit: Vec<f64> = (0..200)
                    .map(|i| (i as f64 / 199.0) * max_depth * 1.2)
                    .collect();
                let y_fit: Vec<f64> = x_fit.iter().map(|&x| chi_square_distribution_simple(x, a, scale)).collect();
                let mut target_y = calculate_min_area_line(&y_fit);
                
                if target_y <= 1e-6 {
                    // Fallback to median mutation rate
                    let mut rates_sorted = mutation_rates.clone();
                    rates_sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
                    target_y = rates_sorted[rates_sorted.len() / 2];
                    if let Some(ref mut logger) = logger {
                        logger.log(&format!("Warning: Target line too small, using data median: y={:.6}", target_y))?;
                    }
                } else {
                    if let Some(ref mut logger) = logger {
                        logger.log(&format!("Target line (min area): y={:.6}", target_y))?;
                    }
                }
                
                // Calculate effective depth for each item
                // Apply PCR bias correction to quality-filtered depth if quality filtering is enabled
                let mut result_map = HashMap::new();
                for (idx, (_, _, _, item)) in all_items.iter().enumerate() {
                    // Use quality-filtered depth as base for PCR bias correction
                    let base_depth = if min_base_qual > 0 {
                        item.effective_depth
                    } else {
                        item.depth
                    };
                    
                    let mutation_rate = if base_depth > 0 {
                        item.rfc as f64 / base_depth as f64
                    } else {
                        0.0
                    };
                    
                    if mutation_rate <= 0.0 || mutation_rate > 0.1 {
                        result_map.insert(idx, base_depth);
                    } else {
                        let correction_factor = calculate_correction_factor(
                            base_depth as f64,
                            mutation_rate,
                            Some((a, scale)),
                            target_y,
                            weight_increase,
                            weight_decrease,
                        );
                        
                        let final_correction = correction_factor.max(0.5).min(10.0);
                        let effective_depth = ((base_depth as f64 * final_correction).round() as usize).max(1);
                        result_map.insert(idx, effective_depth);
                    }
                }
                
                result_map
            } else {
                // Fit failed, use quality-filtered depths (or raw depths if quality filtering disabled)
                if let Some(ref mut logger) = logger {
                    logger.log("Warning: Failed to fit Chi-Square distribution, using quality-filtered depths")?;
                }
                let mut result_map = HashMap::new();
                for (idx, (_, _, _, item)) in all_items.iter().enumerate() {
                    let base_depth = if min_base_qual > 0 {
                        item.effective_depth
                    } else {
                        item.depth
                    };
                    result_map.insert(idx, base_depth);
                }
                result_map
            }
        } else {
            // Not enough data, use quality-filtered depths (or raw depths if quality filtering disabled)
            let mut result_map = HashMap::new();
            for (idx, (_, _, _, item)) in all_items.iter().enumerate() {
                let base_depth = if min_base_qual > 0 {
                    item.effective_depth
                } else {
                    item.depth
                };
                result_map.insert(idx, base_depth);
            }
            result_map
        }
    } else {
        HashMap::new()
    };

    // Write data
    if pcr_bias_correction && !all_items.is_empty() {
        // Write using pre-calculated effective depths (PCR bias correction applied to quality-filtered depth)
        // Only output filtered data (mutation_rate > 0 and <= 0.1) to match Python script
        for (idx, (ref_name, strand, pos, item)) in all_items.iter().enumerate() {
            // Use quality-filtered effective_depth as base for PCR bias correction
            let base_depth = if min_base_qual > 0 {
                item.effective_depth  // Use quality-filtered depth
            } else {
                item.depth  // Use raw depth if quality filtering disabled
            };
            
            let mutation_rate = if base_depth > 0 {
                item.rfc as f64 / base_depth as f64
            } else {
                0.0
            };
            
            // Filter: only output data with valid mutation_rate (matches Python script)
            if mutation_rate > 0.0 && mutation_rate <= 0.1 {
                // Apply PCR bias correction to quality-filtered depth
                let effective_depth = if min_base_qual > 0 {
                    // PCR bias correction on quality-filtered depth
                    effective_depth_map.get(&idx).copied().unwrap_or(item.effective_depth)
                } else {
                    // PCR bias correction on raw depth
                    effective_depth_map.get(&idx).copied().unwrap_or(item.depth)
                };
                
                write_buffer.push_str(&format!(
                    "{},{},{},{},{},{},{},{},{},{},{},{},{},{}\n",
                    ref_name,
                    strand,
                    pos + 1,
                    item.refbase,
                    item.rfc,
                    item.pipc,
                    item.depth,
                    effective_depth,
                    item.ins,
                    item.del,
                    item.a,
                    item.c,
                    item.g,
                    item.t
                ));
                
                if write_buffer.len() > 512 * 1024 {
                    output.write_all(write_buffer.as_bytes())?;
                    write_buffer.clear();
                }
            }
        }
    } else if min_base_qual > 0 {
        // Quality filtering enabled but no PCR bias correction - output quality-filtered effective_depth
        for (i, (ref_name, ref_seq)) in reference_sequences.iter().enumerate() {
            write_progress.update(i)?;
            let empty_vec = Vec::new();
            let actual_positions = region_positions.get(ref_name).unwrap_or(&empty_vec);

            if actual_positions.is_empty() {
                continue;
            } else {
                for &pos in actual_positions {
                    let ref_base = if pos < ref_seq.len() as u32 {
                        ref_seq[pos as usize]
                    } else {
                        'N'
                    };

                    if collect_pos {
                        let pos_item = if let Some(items) = pos_items_map.get(ref_name) {
                            items.get(&pos).cloned().unwrap_or_else(|| {
                                let (pi, _) = itemcreate(ref_name.clone(), pos, ref_base);
                                pi
                            })
                        } else {
                            let (pi, _) = itemcreate(ref_name.clone(), pos, ref_base);
                            pi
                        };

                        if pos_item.depth > 0 {
                            write_buffer.push_str(&format!(
                                "{},{},{},{},{},{},{},{},{},{},{},{},{},{}\n",
                                ref_name,
                                pos_item.strand,
                                pos_item.position + 1,
                                pos_item.refbase,
                                pos_item.rfc,
                                pos_item.pipc,
                                pos_item.depth,
                                pos_item.effective_depth,  // Quality-filtered depth
                                pos_item.ins,
                                pos_item.del,
                                pos_item.a,
                                pos_item.c,
                                pos_item.g,
                                pos_item.t
                            ));
                        }
                    }

                    if collect_neg {
                        let neg_item = if let Some(items) = neg_items_map.get(ref_name) {
                            items.get(&pos).cloned().unwrap_or_else(|| {
                                let (_, ni) = itemcreate(ref_name.clone(), pos, ref_base);
                                ni
                            })
                        } else {
                            let (_, ni) = itemcreate(ref_name.clone(), pos, ref_base);
                            ni
                        };

                        if neg_item.depth > 0 {
                            write_buffer.push_str(&format!(
                                "{},{},{},{},{},{},{},{},{},{},{},{},{},{}\n",
                                ref_name,
                                neg_item.strand,
                                neg_item.position + 1,
                                neg_item.refbase,
                                neg_item.rfc,
                                neg_item.pipc,
                                neg_item.depth,
                                neg_item.effective_depth,  // Quality-filtered depth
                                neg_item.ins,
                                neg_item.del,
                                neg_item.a,
                                neg_item.c,
                                neg_item.g,
                                neg_item.t
                            ));
                        }
                    }
                    
                    if write_buffer.len() > 512 * 1024 {
                        output.write_all(write_buffer.as_bytes())?;
                        write_buffer.clear();
                    }
                }
            }
        }
    } else {
        // Original writing logic (no PCR bias correction or not enough data)
        for (i, (ref_name, ref_seq)) in reference_sequences.iter().enumerate() {
            write_progress.update(i)?;
            let empty_vec = Vec::new();
            let actual_positions = region_positions.get(ref_name).unwrap_or(&empty_vec);

            if actual_positions.is_empty() {
                continue;
            } else {
                for &pos in actual_positions {
                    let ref_base = if pos < ref_seq.len() as u32 {
                        ref_seq[pos as usize]
                    } else {
                        'N'
                    };

                    if collect_pos {
                        let pos_item = if let Some(items) = pos_items_map.get(ref_name) {
                            items.get(&pos).cloned().unwrap_or_else(|| {
                                let (pi, _) = itemcreate(ref_name.clone(), pos, ref_base);
                                pi
                            })
                        } else {
                            let (pi, _) = itemcreate(ref_name.clone(), pos, ref_base);
                            pi
                        };

                        if pos_item.depth > 0 {
                            write_buffer.push_str(&format!(
                                "{},{},{},{},{},{},{},{},{},{},{},{},{}\n",
                                ref_name,
                                pos_item.strand,
                                pos_item.position + 1,
                                pos_item.refbase,
                                pos_item.rfc,
                                pos_item.pipc,
                                pos_item.depth,
                                pos_item.ins,
                                pos_item.del,
                                pos_item.a,
                                pos_item.c,
                                pos_item.g,
                                pos_item.t
                            ));
                        }
                    }

                    if collect_neg {
                        let neg_item = if let Some(items) = neg_items_map.get(ref_name) {
                            items.get(&pos).cloned().unwrap_or_else(|| {
                                let (_, ni) = itemcreate(ref_name.clone(), pos, ref_base);
                                ni
                            })
                        } else {
                            let (_, ni) = itemcreate(ref_name.clone(), pos, ref_base);
                            ni
                        };

                        if neg_item.depth > 0 {
                            write_buffer.push_str(&format!(
                                "{},{},{},{},{},{},{},{},{},{},{},{},{}\n",
                                ref_name,
                                neg_item.strand,
                                neg_item.position + 1,
                                neg_item.refbase,
                                neg_item.rfc,
                                neg_item.pipc,
                                neg_item.depth,
                                neg_item.ins,
                                neg_item.del,
                                neg_item.a,
                                neg_item.c,
                                neg_item.g,
                                neg_item.t
                            ));
                        }
                    }

                    if write_buffer.len() > 512 * 1024 {
                        output.write_all(write_buffer.as_bytes())?;
                        write_buffer.clear();
                    }
                }
            }
        }
    }

    // Final flush to ensure display reaches 100%
    write_progress.finish()?;

    // Write remaining buffer content
    if !write_buffer.is_empty() {
        output.write_all(write_buffer.as_bytes())?;
    }

    // Ensure all data is written to disk
    output.flush()?;
    // Ensure output is flushed to disk before printing newline completely, avoid terminal hanging
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

/// Ultra-efficient region chunk processing using optimized I/O strategy
fn process_region_chunk_ultra_fast(
    bam_path: &str,
    reference_sequences: &HashMap<String, Vec<char>>,
    windows: Vec<WindowTask>,
    strand_filter: &str,
    min_base_qual: u8,
) -> Result<
    (
        HashMap<String, HashMap<u32, Item>>,
        HashMap<String, HashMap<u32, Item>>,
        HashMap<String, Vec<u32>>,
    ),
    Box<dyn Error + Send + Sync>,
> {
    use rust_htslib::bam::IndexedReader;

    // Use optimized BAM reader
    let mut bam = IndexedReader::from_path(bam_path)?;

    // Set BAM reader parameters to optimize performance
    bam.set_threads(1); // Avoid internal multi-threading conflicts

    let header = bam.header().to_owned();

    // Pre-allocate memory to reduce dynamic allocation
    let mut pos_items_map: HashMap<String, HashMap<u32, Item>> = HashMap::new();
    let mut neg_items_map: HashMap<String, HashMap<u32, Item>> = HashMap::new();
    let mut region_positions: HashMap<String, Vec<u32>> = HashMap::new();
    let strand_mode = normalize_strand_filter(strand_filter);
    let collect_pos = strand_mode != "-";
    let collect_neg = strand_mode != "+";

    // Query and process each window
    for w in &windows {
        if !reference_sequences.contains_key(&w.region) {
            continue;
        }
        let ref_seq = &reference_sequences[&w.region];
        let ref_len = ref_seq.len() as u32;

        // Initialize container for this region (if not exists)
        if !pos_items_map.contains_key(&w.region) {
            pos_items_map.insert(w.region.clone(), HashMap::new());
        }
        if !neg_items_map.contains_key(&w.region) {
            neg_items_map.insert(w.region.clone(), HashMap::new());
        }
        if !region_positions.contains_key(&w.region) {
            region_positions.insert(w.region.clone(), Vec::new());
        }

        // Create independent BAM reader for each window
        let mut window_bam = IndexedReader::from_path(bam_path)?;
        window_bam.set_threads(1);

        // Use window query, only process data in current window range
        let start1 = (w.start + 1) as u32; // 1-based inclusive
        let end1 = w.end as u32; // 1-based inclusive end
        let region_query = format!("{}:{}-{}", w.region, start1, end1);

        // Query all reads for this window
        if let Ok(()) = window_bam.fetch(&region_query) {
            // Get header
            let window_header = window_bam.header().to_owned();

            // Create pileup for this region, using optimized parameters
            let mut pileups = window_bam.pileup();
            pileups.set_max_depth(1_000_000);

            // Pre-allocate temporary data structures
            let mut temp_positions = Vec::with_capacity(1000);
            let mut temp_pos_items = HashMap::new();
            let mut temp_neg_items = HashMap::new();

            // Process pileup data for this region
            for pileup in pileups {
                let pileup = pileup?;
                let ref_name =
                    String::from_utf8_lossy(window_header.tid2name(pileup.tid())).to_string();

                // Ensure only processing current region
                if ref_name != w.region {
                    continue;
                }

                let pos = pileup.pos();

                // Check if position is within valid range
                if pos >= ref_len {
                    continue;
                }
                // Filter to current window range (use closed interval [start, end], allow stop to cross window boundaries)
                if pos < w.start || pos > w.end {
                    continue;
                }

                let ref_base = ref_seq[pos as usize];

                // Record sites for current region (deduplicate)
                if !temp_positions.contains(&pos) {
                    temp_positions.push(pos);
                }

                if collect_pos && !temp_pos_items.contains_key(&pos) {
                    let (pi, _) = itemcreate(w.region.clone(), pos, ref_base);
                    temp_pos_items.insert(pos, pi);
                }
                if collect_neg && !temp_neg_items.contains_key(&pos) {
                    let (_, ni) = itemcreate(w.region.clone(), pos, ref_base);
                    temp_neg_items.insert(pos, ni);
                }

                // Process each alignment
                for alignment in pileup.alignments() {
                    let record = alignment.record();
                    let id = alignment.indel();

                    if let Some(qpos) = alignment.qpos() {
                        let base = findchar(qpos, record.seq());
                        let numb = chartonum(ref_base);
                        let is_mutation = numb != base;
                        
                        // Check base quality for depth calculation (similar to shapemapper2)
                        // Quality filtering affects both mutation counting and effective_depth
                        let quals = record.qual();
                        let has_good_quality = if min_base_qual > 0 {
                            // Check current position quality
                            if qpos >= 0 && (qpos as usize) < quals.len() {
                                let base_qual = quals[qpos as usize];
                                // rust_htslib returns decoded Phred scores (0-93), not ASCII
                                // So we compare directly without adding 33
                                (base_qual as u8) >= min_base_qual
                            } else {
                                false  // No quality information for this position or invalid qpos
                            }
                        } else {
                            true  // Quality filtering disabled
                        };
                        
                        // For mutations, check quality (matching shapemapper2's filterQscoresCountDepths)
                        let should_count_mutation = if is_mutation {
                            has_good_quality  // Mutation quality check uses same logic as depth
                        } else {
                            true  // Not a mutation
                        };
                        
                        // Note: update_item closure removed - we update counts directly below

                        // Accumulate stop signal
                        // Theoretical correction: positive strand start_site - 1, negative strand start_site + 1
                        if alignment.is_head() {
                            let stop_signal_pos = if record.is_reverse() {
                                // Negative strand: correct to start_site + 1
                                pos + 1
                            } else {
                                // Positive strand: correct to start_site - 1
                                if pos > 0 { pos - 1 } else { 0 }
                            };
                            
                            // Check if corrected position is within valid range
                            if stop_signal_pos < ref_len {
                                let ref_base_stop = ref_seq[stop_signal_pos as usize];

                                if record.is_reverse() {
                                    if collect_neg {
                                        let entry = temp_neg_items
                                            .entry(stop_signal_pos)
                                            .or_insert_with(|| {
                                                let (_, ni) = itemcreate(
                                                    w.region.clone(),
                                                    stop_signal_pos,
                                                    ref_base_stop,
                                                );
                                                ni
                                            });
                                        entry.pipc += 1;
                                        if !temp_positions.contains(&stop_signal_pos) {
                                            temp_positions.push(stop_signal_pos);
                                        }
                                    }
                                } else if collect_pos {
                                    let entry =
                                        temp_pos_items.entry(stop_signal_pos).or_insert_with(|| {
                                            let (pi, _) = itemcreate(
                                                w.region.clone(),
                                                stop_signal_pos,
                                                ref_base_stop,
                                            );
                                            pi
                                        });
                                    entry.pipc += 1;
                                    if !temp_positions.contains(&stop_signal_pos) {
                                        temp_positions.push(stop_signal_pos);
                                    }
                                }
                            }
                        }

                        // Update other statistics for current position
                        // Update directly without closure to ensure correct mutation counting
                        if record.is_reverse() {
                            if collect_neg {
                                let entry = temp_neg_items.entry(pos).or_insert_with(|| {
                                    let (_, ni) = itemcreate(w.region.clone(), pos, ref_base);
                                    ni
                                });
                                match base {
                                    1 => entry.a += 1,
                                    2 => entry.c += 1,
                                    4 => entry.g += 1,
                                    8 => entry.t += 1,
                                    _ => {}
                                }
                                match id {
                                    rust_htslib::bam::pileup::Indel::Ins(_) => {
                                        entry.ins += 1;
                                    }
                                    rust_htslib::bam::pileup::Indel::Del(_) => {
                                        entry.del += 1;
                                        entry.depth += 1;
                                    }
                                    _ => {}
                                }
                                entry.depth += 1;
                                // Count effective_depth only if quality check passed (similar to shapemapper2)
                                if has_good_quality {
                                    entry.effective_depth += 1;
                                }
                                // Only count mutation if quality check passed
                                if is_mutation && should_count_mutation {
                                    entry.rfc += 1;
                                }
                            }
                        } else if collect_pos {
                            let entry = temp_pos_items.entry(pos).or_insert_with(|| {
                                let (pi, _) = itemcreate(w.region.clone(), pos, ref_base);
                                pi
                            });
                            match base {
                                1 => entry.a += 1,
                                2 => entry.c += 1,
                                4 => entry.g += 1,
                                8 => entry.t += 1,
                                _ => {}
                            }
                            match id {
                                rust_htslib::bam::pileup::Indel::Ins(_) => {
                                    entry.ins += 1;
                                }
                                rust_htslib::bam::pileup::Indel::Del(_) => {
                                    entry.del += 1;
                                    entry.depth += 1;
                                }
                                _ => {}
                            }
                            entry.depth += 1;
                            // Count effective_depth only if quality check passed (similar to shapemapper2)
                            if has_good_quality {
                                entry.effective_depth += 1;
                            }
                            // Only count mutation if quality check passed
                            if is_mutation && should_count_mutation {
                                entry.rfc += 1;
                            }
                        }
                    }
                }
            }

            // Batch update main data structure (merge to region level), no window-level filtering to avoid boundary sites being ignored
            let pos_map = pos_items_map.get_mut(&w.region).unwrap();
            let neg_map = neg_items_map.get_mut(&w.region).unwrap();
            for (k, v) in temp_pos_items {
                pos_map.insert(k, v);
            }
            for (k, v) in temp_neg_items {
                neg_map.insert(k, v);
            }

            // Sort and update positions (unconditionally merge sites observed within window)
            let mut filtered_positions = temp_positions;
            filtered_positions.sort();
            let pos_vec = region_positions.get_mut(&w.region).unwrap();

            // Debug information removed

            // Merge and deduplicate (simple method, global sorting later)
            for p in filtered_positions {
                if !pos_vec.contains(&p) {
                    pos_vec.push(p);
                }
            }
        }
    }

    Ok((pos_items_map, neg_items_map, region_positions))
}
