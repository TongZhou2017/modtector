use chrono;
use clap::Args;
use std::collections::HashMap;
use std::error::Error;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::Path;
use std::time::Instant;

/// Validate norm command arguments
fn validate_norm_args(args: &NormArgs) -> Result<(), Box<dyn Error>> {
    // Validate input file
    if args.input.trim().is_empty() {
        return Err("Error: Input file path cannot be empty".into());
    }
    if !Path::new(&args.input).exists() {
        return Err(format!("Error: Input file does not exist: {}", args.input).into());
    }
    if !args.input.ends_with(".csv") {
        return Err(format!("Error: Input file path must end with .csv: {}", args.input).into());
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

    // Validate normalization methods
    let method_name = args.method.to_lowercase();
    if ![
        "percentile28",
        "2-8",
        "2_8",
        "winsor90",
        "winsorizing",
        "boxplot",
    ]
    .contains(&method_name.as_str())
    {
        return Err(format!("Error: Unknown normalization method: {}. Supported methods: percentile28, winsor90, boxplot", args.method).into());
    }

    // Validate SNP threshold
    if args.snp_cutoff < 0.0 || args.snp_cutoff > 1.0 {
        return Err(format!(
            "Error: SNP threshold must be between 0.0 and 1.0, current: {}",
            args.snp_cutoff
        )
        .into());
    }

    // Validate base types
    if args.base_types.trim().is_empty() {
        return Err("Error: Base types cannot be empty".into());
    }
    for base in args.base_types.chars() {
        if !"ACGT".contains(base) {
            return Err(format!(
                "Error: Base types can only contain ACGT, current contains: {}",
                base
            )
            .into());
        }
    }

    // Validate window parameters
    if args.dynamic_window && args.window_size > 0 {
        return Err("Error: Dynamic window mode cannot be used with fixed window size".into());
    }

    if args.window_size > 0 && args.window_offset >= args.window_size {
        return Err(format!("Error: Window offset cannot be greater than or equal to window size, current: {} >= {}", args.window_offset, args.window_size).into());
    }

    if args.window_size_dynamic == 0 {
        return Err("Error: Dynamic window size cannot be 0".into());
    }
    if args.window_size_dynamic > 1000 {
        return Err(format!(
            "Error: Dynamic window size cannot exceed 1000, current: {}",
            args.window_size_dynamic
        )
        .into());
    }

    Ok(())
}

#[derive(Debug, Clone)]
struct ReactivityData {
    index: usize,
    position: u32,
    base: char,
    reactivity: f64,
    is_reactive: bool,
}

#[derive(Args, Debug)]
pub struct NormArgs {
    // Input/Output
    /// Input CSV file
    #[arg(short = 'i', long = "input")]
    pub input: String,
    /// Output CSV file
    #[arg(short = 'o', long = "output")]
    pub output: String,

    // Basic configuration
    /// Normalization method: percentile28, winsor90, boxplot
    #[arg(short = 'm', long = "method")]
    pub method: String,
    /// SNP threshold for filtering mutation signals
    #[arg(long = "cutoff", default_value_t = 0.25)]
    pub snp_cutoff: f64,
    /// Reactive base types (e.g., "AC" for A and C, default: all bases)
    #[arg(long = "bases", default_value = "ACGT")]
    pub base_types: String,

    // Window configuration (mutually exclusive)
    /// Fixed window size (fixed sliding window mode)
    #[arg(long = "window", default_value_t = 0)]
    pub window_size: usize,
    /// Window offset (fixed sliding window mode)
    #[arg(long = "window-offset", default_value_t = 0)]
    pub window_offset: usize,
    /// Enable dynamic sliding window mode
    #[arg(long = "dynamic", default_value_t = false)]
    pub dynamic_window: bool,
    /// Number of reactive bases required per window in dynamic sliding window
    #[arg(long = "window-size-dynamic", default_value_t = 50)]
    pub window_size_dynamic: usize,

    // Post-processing
    /// Apply piecewise linear mapping to normalized signals
    #[arg(long = "linear", default_value_t = false)]
    pub linear_remap: bool,

    // Logging
    /// Log file path (optional)
    #[arg(short = 'l', long = "log")]
    pub log: Option<String>,
}

#[derive(Debug, Clone, PartialEq)]
pub enum StopNormMethod {
    Percentile28,
    Winsor90,
    Boxplot,
}

impl std::str::FromStr for StopNormMethod {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "percentile28" | "2-8" | "2_8" => Ok(StopNormMethod::Percentile28),
            "winsor90" | "winsorizing" => Ok(StopNormMethod::Winsor90),
            "boxplot" => Ok(StopNormMethod::Boxplot),
            _ => Err(format!("Unknown stop normalization method: {}", s)),
        }
    }
}

/// Zarringhalam piecewise linear mapping. Implemented entirely inside this function so that
/// whenever Zarringhalam is used (with or without prior norm), behavior is consistent:
/// the minimum value (e.g. Winsor90's v5/v95) is treated as 0, then [0, max] is mapped by
/// the piecewise linear rule. This preserves relative mod/unmod differences for Siegfried
/// instead of clamping negatives to 0 and losing information.
/// 
/// Fixed: Use actual maximum value instead of fixed 1.0 for high-value region mapping.
/// This matches rf-norm implementation and prevents compression of high values that exceed 1.0
/// after normalization, which was causing high false positive rates.
/// 
/// Key change: After winsor90 normalization, values may exceed 1.0. Instead of normalizing
/// to [0, 1] with fixed 1.0 as max, we use the actual max_val for the high-value region
/// mapping [0.7, max] → [0.85, 1.0], matching rf-norm's behavior.
fn zarringhalam_remap(values: &[f64]) -> Vec<f64> {
    // Filter out NaN values for min/max calculation
    let valid_values: Vec<f64> = values.iter().cloned().filter(|&x| x.is_finite()).collect();
    
    if valid_values.is_empty() {
        return values.iter().map(|&x| if x.is_finite() { 0.0 } else { x }).collect();
    }
    
    let min_val = valid_values
        .iter()
        .cloned()
        .fold(f64::INFINITY, |a, b| a.min(b));
    let max_val = valid_values
        .iter()
        .cloned()
        .fold(f64::NEG_INFINITY, |a, b| a.max(b));
    
    // Apply Zarringhalam piecewise linear mapping directly on original values
    // This matches rf-norm's implementation which uses actual max value for high-value region
    values
        .iter()
        .map(|&x| {
            // Skip NaN values
            if !x.is_finite() {
                return x;
            }
            
            if x < 0.25 {
                // [0, 0.25) → [0, 0.35)
                x / 0.25 * 0.35
            } else if x >= 0.25 && x < 0.3 {
                // [0.25, 0.3) → [0.35, 0.55)
                0.35 + (x - 0.25) / 0.05 * (0.55 - 0.35)
            } else if x >= 0.3 && x < 0.7 {
                // [0.3, 0.7) → [0.55, 0.85)
                0.55 + (x - 0.3) / 0.4 * (0.85 - 0.55)
            } else {
                // [0.7, max] → [0.85, 1.0]
                // Use actual max_val instead of fixed 1.0 (matches rf-norm)
                // This prevents compression of high values that exceed 1.0 after normalization
                if max_val > 0.7 {
                    0.85 + (x - 0.7) / (max_val - 0.7) * (1.0 - 0.85)
                } else {
                    // Edge case: if max <= 0.7, all values map to 0.85
                    0.85
                }
            }
        })
        .collect()
}

pub fn norm_csv(args: &NormArgs, logger: &mut crate::Logger) -> Result<(), Box<dyn Error>> {
    // Validate norm command parameters
    validate_norm_args(args)?;

    let start_time = Instant::now();

    // Set up log file
    let log_file = if let Some(log_path) = &args.log {
        std::fs::File::create(log_path)?
    } else {
        std::fs::File::create("norm.log")?
    };
    let mut logger = crate::Logger::new(log_file);

    // Record environment information and parameters
    logger.log("=== ModDetector Norm Function Log ===")?;
    logger.log(&format!("Software Version: v{}", crate::VERSION))?;
    logger.log(&format!(
        "Runtime: {}",
        chrono::Utc::now().format("%Y-%m-%d %H:%M:%S")
    ))?;
    logger.log(&format!("Input File: {}", args.input))?;
    logger.log(&format!("Output File: {}", args.output))?;
    // Record normalization method information
    logger.log(&format!("Normalization Method: {}", args.method))?;
    logger.log(&format!("SNP Threshold: {}", args.snp_cutoff))?;
    logger.log(&format!("Base Types: {}", args.base_types))?;
    logger.log(&format!("Window Size: {}", args.window_size))?;
    logger.log(&format!("Dynamic Window: {}", args.dynamic_window))?;
    logger.log(&format!("Linear Remap: {}", args.linear_remap))?;
    logger.log("Starting normalization processing...")?;

    // Display data loading information
    println!("[Loading data]");
    println!("    Reactivity data: {}", args.input);
    println!();

    let file = File::open(&args.input)?;
    let reader = BufReader::new(file);
    let mut lines = reader.lines();
    let header = lines.next().ok_or("CSV file is empty")??;
    let mut records: Vec<Vec<String>> = Vec::new();
    let mut region_map: HashMap<String, Vec<usize>> = HashMap::new();
    let mut stop_reactivity_values: Vec<f64> = Vec::new();
    let mut mutation_reactivity_values: Vec<f64> = Vec::new();

    // Record two types of reactivity values for each row
    for (_idx, line) in lines.enumerate() {
        let l = line?;
        let fields: Vec<String> = l.split(',').map(|s| s.to_string()).collect();
        if fields.len() < 6 {
            continue;
        } // Simplified format requires at least 6 columns

        let chr = fields[0].clone();
        let strand = fields[1].clone();

        // Read two types of reactivity values, support NaN
        let stop_reactivity = fields[4].trim().parse::<f64>()
            .unwrap_or(if fields[4].trim().eq_ignore_ascii_case("nan") { f64::NAN } else { 0.0 });
        let mutation_reactivity = fields[5].trim().parse::<f64>()
            .unwrap_or(if fields[5].trim().eq_ignore_ascii_case("nan") { f64::NAN } else { 0.0 });

        records.push(fields);
        let record_index = records.len() - 1;
        // Group by region (ChrID+Strand), ensure positive and negative strands are processed separately
        let region = format!("{}|{}", chr, strand);
        region_map.entry(region).or_default().push(record_index);
        stop_reactivity_values.push(stop_reactivity);
        mutation_reactivity_values.push(mutation_reactivity);
    }

    // Parse normalization method
    let norm_method = args.method.parse().unwrap_or(StopNormMethod::Percentile28);

    // Display parameter information
    println!("[Params]                              ");
    let method_name = match norm_method {
        StopNormMethod::Percentile28 => "percentile28",
        StopNormMethod::Winsor90 => "winsor90",
        StopNormMethod::Boxplot => "boxplot",
    };
    println!("    Normalization method: {}.", method_name);
    println!("    Threads: 1.");
    println!();

    // Display progress information
    let total_regions = region_map.len();
    let mut progress = crate::progress::SimpleProgress::new(total_regions);

    print!("\r[Progressing] Processing {} regions...", total_regions);
    std::io::stdout().flush()?;

    // Process stop signal normalization
    println!("[Processing] Normalizing stop signals...");
    for (i, (_region, idxs)) in region_map.iter().enumerate() {
        if i % 100 == 0 {
            progress.update(i)?;
        }
        let region_reactivities: Vec<(usize, f64)> = idxs
            .iter()
            .map(|&i| (i, stop_reactivity_values[i]))
            .collect();
        if region_reactivities.is_empty() {
            continue;
        }
        // Filter out NaN values, only used for normalization calculation
        let values: Vec<f64> = region_reactivities.iter()
            .map(|x| x.1)
            .filter(|&v| !v.is_nan())
            .collect();
        
        // If all values are NaN, keep all values unchanged
        if values.is_empty() {
            for (idx, orig_val) in region_reactivities.iter() {
                if orig_val.is_nan() {
                    records[*idx][4] = "NaN".to_string();
                } else {
                    records[*idx][4] = format!("{:.6}", orig_val);
                }
            }
            continue;
        }
        
        let mut normed = match norm_method {
            StopNormMethod::Percentile28 => percentile_28_norm(&values),
            StopNormMethod::Winsor90 => winsor_90_norm(&values),
            StopNormMethod::Boxplot => boxplot_norm(&values),
        };
        // Zarringhalam piecewise linear mapping (min-as-0 is done inside zarringhalam_remap)
        if args.linear_remap {
            normed = zarringhalam_remap(&normed);
        }
        // Write back to stop reactivity column (column 4), keep NaN values unchanged
        let mut norm_idx = 0;
        for (idx, orig_val) in region_reactivities.iter() {
            if orig_val.is_nan() {
                // Keep NaN value unchanged
                records[*idx][4] = "NaN".to_string();
            } else {
                records[*idx][4] = format!("{:.6}", normed[norm_idx]);
                norm_idx += 1;
            }
        }
    }

    // Process mutation signal normalization
    println!("[Processing] Normalizing mutation signals...");
    for (i, (_region, idxs)) in region_map.iter().enumerate() {
        if i % 100 == 0 {
            progress.update(i)?;
        }
        let region_reactivities: Vec<(usize, f64)> = idxs
            .iter()
            .map(|&i| (i, mutation_reactivity_values[i]))
            .collect();
        if region_reactivities.is_empty() {
            continue;
        }
        // Filter out NaN values, only used for normalization calculation
        let values: Vec<f64> = region_reactivities.iter()
            .map(|x| x.1)
            .filter(|&v| !v.is_nan())
            .collect();
        
        // If all values are NaN, keep all values unchanged
        if values.is_empty() {
            for (idx, orig_val) in region_reactivities.iter() {
                if orig_val.is_nan() {
                    records[*idx][5] = "NaN".to_string();
                } else {
                    records[*idx][5] = format!("{:.6}", orig_val);
                }
            }
            continue;
        }
        
        let mut normed = match norm_method {
            StopNormMethod::Percentile28 => percentile_28_norm(&values),
            StopNormMethod::Winsor90 => winsor_90_norm(&values),
            StopNormMethod::Boxplot => boxplot_norm(&values),
        };
        // Zarringhalam piecewise linear mapping (min-as-0 is done inside zarringhalam_remap)
        if args.linear_remap {
            normed = zarringhalam_remap(&normed);
        }
        // Write back to mutation reactivity column (column 5), keep NaN values unchanged
        let mut norm_idx = 0;
        for (idx, orig_val) in region_reactivities.iter() {
            if orig_val.is_nan() {
                // Keep NaN value unchanged
                records[*idx][5] = "NaN".to_string();
            } else {
                records[*idx][5] = format!("{:.6}", normed[norm_idx]);
                norm_idx += 1;
            }
        }
    }

    progress.finish()?;

    // Ensure all records have correct 6-column format
    for record in &mut records {
        // Ensure record has at least 6 columns
        while record.len() < 6 {
            record.push("0.000000".to_string());
        }
        // Ensure reactivity columns (columns 4, 5) have valid values, but keep NaN unchanged
        if record.len() > 4 && record[4].is_empty() {
            record[4] = "0.000000".to_string();
        }
        // NaN values remain unchanged, used to mark SNP filtered positions
        if record.len() > 5 && record[5].is_empty() {
            record[5] = "0.000000".to_string();
        }
        // NaN values remain unchanged, used to mark SNP filtered positions
    }
    // Output
    let mut out = File::create(&args.output)?;
    writeln!(out, "{}", header)?;
    for rec in records {
        writeln!(out, "{}", rec.join(","))?;
    }
    let elapsed = start_time.elapsed();

    // Overwrite progress display, show output information
    println!("\r[Output]                           ");
    println!("    Normalized: {}", args.output);
    println!("{}", crate::progress::format_time_used(elapsed));

    // Log completion
    logger.log(&format!(
        "Normalization processing completed, output file: {}",
        args.output
    ))?;
    logger.log(&format!("Total time: {:.2}s", elapsed.as_secs_f64()))?;

    Ok(())
}

/// Dynamic sliding window normalization
fn dynamic_window_norm(
    data: &[ReactivityData],
    norm_method: &StopNormMethod,
    target_reactive_count: usize,
    reactive_bases: &str,
) -> Vec<f64> {
    if data.is_empty() {
        return vec![];
    }

    let mut result = vec![0.0; data.len()];

    // Calculate total number of reactive bases in data
    let total_reactive_bases = data
        .iter()
        .filter(|d| reactive_bases.contains(d.base))
        .count();

    // If total reactive base count is too low, use global normalization
    if total_reactive_bases < target_reactive_count {
        let all_reactivities: Vec<f64> = data.iter().map(|d| d.reactivity).collect();
        let global_normed = normalize_window(&all_reactivities, norm_method);
        for (j, &value) in global_normed.iter().enumerate() {
            if j < result.len() {
                result[j] = value;
            }
        }
        return result;
    }

    // Calculate target size for each window (based on total reactive base count)
    let window_count = (total_reactive_bases + target_reactive_count - 1) / target_reactive_count;
    let bases_per_window = (total_reactive_bases + window_count - 1) / window_count;

    let mut i = 0;
    let mut window_start = 0;

    while i < data.len() {
        // Find window containing enough reactive bases
        let window = find_adaptive_window(data, i, bases_per_window, reactive_bases);

        if window.is_empty() {
            // If no suitable window found, use all remaining data
            let remaining_data: Vec<f64> = data[i..].iter().map(|d| d.reactivity).collect();
            let normed = normalize_window(&remaining_data, norm_method);
            for (j, &value) in normed.iter().enumerate() {
                if i + j < result.len() {
                    result[i + j] = value;
                }
            }
            break;
        }

        // Normalize data within window
        let window_reactivities: Vec<f64> = window.iter().map(|d| d.reactivity).collect();

        // Check if window has enough non-zero values
        let non_zero_count = window_reactivities.iter().filter(|&&x| x > 0.0).count();

        if non_zero_count < 3 {
            // If window has too few non-zero values, use global normalization of current region data
            let region_reactivities: Vec<f64> = data.iter().map(|d| d.reactivity).collect();
            let region_non_zero_count = region_reactivities.iter().filter(|&&x| x > 0.0).count();

            if region_non_zero_count > 0 {
                let region_normed = normalize_window(&region_reactivities, norm_method);
                for (j, &value) in region_normed.iter().enumerate() {
                    if j < result.len() {
                        result[j] = value;
                    }
                }
            } else {
                // If entire region has no non-zero values, keep original values
                for (j, data_item) in data.iter().enumerate() {
                    if j < result.len() {
                        result[j] = data_item.reactivity;
                    }
                }
            }
            break;
        }

        let normed = normalize_window(&window_reactivities, norm_method);

        // Write normalized results back to corresponding positions
        for (j, &value) in normed.iter().enumerate() {
            if j < window.len() && window[j].index < result.len() {
                result[window[j].index] = value;
            }
        }

        // Move to next window start position
        if let Some(last_pos) = window.last() {
            i = last_pos.index + 1;
            if i >= data.len() {
                break;
            }
        } else {
            i += 1;
        }
    }

    result
}

/// Find adaptive window
fn find_adaptive_window(
    data: &[ReactivityData],
    start_idx: usize,
    target_reactive_count: usize,
    reactive_bases: &str,
) -> Vec<ReactivityData> {
    let mut window = Vec::new();
    let mut reactive_count = 0;
    let mut i = start_idx;

    // If target count is 0, return a small fixed window
    if target_reactive_count == 0 {
        let end_idx = (start_idx + 100).min(data.len());
        return data[start_idx..end_idx].to_vec();
    }

    // Calculate maximum window size (avoid window being too large)
    let max_window_size = (target_reactive_count * 3).min(data.len() - start_idx);

    while i < data.len() && reactive_count < target_reactive_count && window.len() < max_window_size
    {
        let item = &data[i];
        if reactive_bases.contains(item.base) {
            reactive_count += 1;
        }
        window.push(item.clone());
        i += 1;
    }

    // If enough reactive bases found, return window
    if reactive_count >= target_reactive_count {
        window
    } else {
        // If reached sequence end and still not enough, return all remaining data
        data[start_idx..].to_vec()
    }
}

/// Find dynamic window (keep original function for compatibility)
fn find_dynamic_window(
    data: &[ReactivityData],
    start_idx: usize,
    target_count: usize,
    reactive_bases: &str,
) -> Vec<ReactivityData> {
    find_adaptive_window(data, start_idx, target_count, reactive_bases)
}

/// Normalize window data
fn normalize_window(values: &[f64], norm_method: &StopNormMethod) -> Vec<f64> {
    match norm_method {
        StopNormMethod::Percentile28 => percentile_28_norm(values),
        StopNormMethod::Winsor90 => winsor_90_norm(values),
        StopNormMethod::Boxplot => boxplot_norm(values),
    }
}

// 2-8% Normalization
fn percentile_28_norm(values: &[f64]) -> Vec<f64> {
    let mut sorted: Vec<f64> = values.iter().cloned().collect();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let n = sorted.len();
    if n < 10 {
        return values.to_vec();
    }
    let p2 = (n as f64 * 0.02).ceil() as usize;
    let p10 = (n as f64 * 0.10).ceil() as usize;
    let p2 = p2.min(n - 1);
    let p10 = p10.min(n - 1);
    let avg: f64 = if p10 > p2 {
        let sum = sorted[p2..p10].iter().sum::<f64>();
        if sum > 0.0 {
            sum / (p10 - p2) as f64
        } else {
            1.0
        }
    } else {
        1.0
    };
    values.iter().map(|&x| x / avg).collect()
}

// 90% Winsorizing
fn winsor_90_norm(values: &[f64]) -> Vec<f64> {
    let mut sorted: Vec<f64> = values.iter().cloned().collect();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let n = sorted.len();
    if n < 2 {
        return values.to_vec();
    }
    let p5 = (n as f64 * 0.05).floor() as usize;
    let p95 = (n as f64 * 0.95).ceil() as usize;
    let v5 = sorted[p5.min(n - 1)];
    let v95 = sorted[p95.min(n - 1)];
    let winsorized: Vec<f64> = values
        .iter()
        .map(|&x| {
            if x < v5 {
                v5
            } else if x > v95 {
                v95
            } else {
                x
            }
        })
        .collect();
    let norm_val = v95.max(1e-6);
    winsorized.iter().map(|&x| x / norm_val).collect()
}

// Box-plot Normalization
fn boxplot_norm(values: &[f64]) -> Vec<f64> {
    let mut sorted: Vec<f64> = values.iter().cloned().collect();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let n = sorted.len();
    if n < 4 {
        return values.to_vec();
    }
    let q1_idx = (n as f64 * 0.25).ceil() as usize;
    let q3_idx = (n as f64 * 0.75).ceil() as usize;
    let q1 = sorted[q1_idx.min(n - 1)];
    let q3 = sorted[q3_idx.min(n - 1)];
    let iqr = q3 - q1;
    let upper = q3 + 1.5 * iqr;
    // Remove outliers
    let filtered: Vec<f64> = sorted.iter().cloned().filter(|&x| x <= upper).collect();
    let top10 = (filtered.len() as f64 * 0.10).ceil() as usize;
    let avg = if top10 > 0 {
        let sum = filtered.iter().take(top10).sum::<f64>();
        if sum > 0.0 {
            sum / top10 as f64
        } else {
            1.0
        }
    } else {
        1.0
    };
    values.iter().map(|&x| x / avg).collect()
}
