use std::error::Error;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use clap::Args;
use std::collections::HashMap;
use std::time::Instant;
use std::path::Path;
use chrono;

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
        return Err(format!("Error: Output file path must end with .csv: {}", args.output).into());
    }
    
    // Validate normalization methods
    let method_name = args.method.to_lowercase();
    if !["percentile28", "2-8", "2_8", "winsor90", "winsorizing", "boxplot"].contains(&method_name.as_str()) {
        return Err(format!("Error: Unknown normalization method: {}. Supported methods: percentile28, winsor90, boxplot", args.method).into());
    }
    
    // Validate SNP threshold
    if args.snp_cutoff < 0.0 || args.snp_cutoff > 1.0 {
        return Err(format!("Error: SNP threshold must be between 0.0 and 1.0, current: {}", args.snp_cutoff).into());
    }
    
    // Validate base types
    if args.base_types.trim().is_empty() {
        return Err("Error: Base types cannot be empty".into());
    }
    for base in args.base_types.chars() {
        if !"ACGT".contains(base) {
            return Err(format!("Error: Base types can only contain ACGT, current contains: {}", base).into());
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
        return Err(format!("Error: Dynamic window size cannot exceed 1000, current: {}", args.window_size_dynamic).into());
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

#[derive(Args)]
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

fn zarringhalam_remap(values: &[f64]) -> Vec<f64> {
    values.iter().map(|&x| {
        if x < 0.25 {
            x / 0.25 * 0.35
        } else if x < 0.3 {
            0.35 + (x - 0.25) / 0.05 * (0.55 - 0.35)
        } else if x < 0.7 {
            0.55 + (x - 0.3) / 0.4 * (0.85 - 0.55)
        } else {
            0.85 + (x - 0.7) / 0.3 * (1.0 - 0.85)
        }
    }).collect()
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
    logger.log(&format!("Runtime: {}", chrono::Utc::now().format("%Y-%m-%d %H:%M:%S")))?;
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
    for (idx, line) in lines.enumerate() {
        let l = line?;
        let fields: Vec<String> = l.split(',').map(|s| s.to_string()).collect();
        if fields.len() < 6 { continue; } // Simplified format requires at least 6 columns
        
        let chr = fields[0].clone();
        let strand = fields[1].clone();
        let pos = fields[2].parse::<u32>().unwrap_or(0);
        let base = if fields.len() > 3 { fields[3].chars().next().unwrap_or('N') } else { 'N' };
        
        // 读取两种reactivity值
        let stop_reactivity = fields[4].parse().unwrap_or(0.0);
        let mutation_reactivity = fields[5].parse().unwrap_or(0.0);
        
        records.push(fields);
        // 区域分组（ChrID+Strand），确保正负链分开处理
        let region = format!("{}|{}", chr, strand);
        region_map.entry(region).or_default().push(idx);
        stop_reactivity_values.push(stop_reactivity);
        mutation_reactivity_values.push(mutation_reactivity);
    }
    
    // 解析归一化方法
    let norm_method = args.method.parse().unwrap_or(StopNormMethod::Percentile28);
    
    // 显示参数信息
    println!("[Params]                              ");
    let method_name = match norm_method {
        StopNormMethod::Percentile28 => "percentile28",
        StopNormMethod::Winsor90 => "winsor90",
        StopNormMethod::Boxplot => "boxplot",
    };
    println!("    Normalization method: {}.", method_name);
    println!("    Threads: 1.");
    println!();
    
    // 显示进度信息
    let total_regions = region_map.len();
    let mut progress = crate::progress::SimpleProgress::new(total_regions);
    
    print!("\r[Progressing] Processing {} regions...", total_regions);
    std::io::stdout().flush()?;
    
    // 处理stop信号归一化
    println!("[Processing] Normalizing stop signals...");
    for (i, (_region, idxs)) in region_map.iter().enumerate() {
        if i % 100 == 0 {
            progress.update(i)?;
        }
        let region_reactivities: Vec<(usize, f64)> = idxs.iter().map(|&i| (i, stop_reactivity_values[i])).collect();
        if region_reactivities.is_empty() { continue; }
        let values: Vec<f64> = region_reactivities.iter().map(|x| x.1).collect();
        let mut normed = match norm_method {
            StopNormMethod::Percentile28 => {
                percentile_28_norm(&values)
            },
            StopNormMethod::Winsor90 => {
                winsor_90_norm(&values)
            },
            StopNormMethod::Boxplot => {
                boxplot_norm(&values)
            },
        };
        // Zarringhalam分段线性映射
        if args.linear_remap {
            normed = zarringhalam_remap(&normed);
        }
        // 写回stop reactivity列（第4列）
        for (j, (idx, _)) in region_reactivities.iter().enumerate() {
            records[*idx][4] = format!("{:.6}", normed[j]);
        }
    }
    
    // 处理mutation信号归一化
    println!("[Processing] Normalizing mutation signals...");
    for (i, (_region, idxs)) in region_map.iter().enumerate() {
        if i % 100 == 0 {
            progress.update(i)?;
        }
        let region_reactivities: Vec<(usize, f64)> = idxs.iter().map(|&i| (i, mutation_reactivity_values[i])).collect();
        if region_reactivities.is_empty() { continue; }
        let values: Vec<f64> = region_reactivities.iter().map(|x| x.1).collect();
        let mut normed = match norm_method {
            StopNormMethod::Percentile28 => {
                percentile_28_norm(&values)
            },
            StopNormMethod::Winsor90 => {
                winsor_90_norm(&values)
            },
            StopNormMethod::Boxplot => {
                boxplot_norm(&values)
            },
        };
        // Zarringhalam分段线性映射
        if args.linear_remap {
            normed = zarringhalam_remap(&normed);
        }
        // 写回mutation reactivity列（第5列）
        for (j, (idx, _)) in region_reactivities.iter().enumerate() {
            records[*idx][5] = format!("{:.6}", normed[j]);
        }
    }
    
    progress.finish()?;
    
    // 确保所有记录都有正确的6列格式
    for record in &mut records {
        // 确保记录至少有6列
        while record.len() < 6 {
            record.push("0.000000".to_string());
        }
        // 确保reactivity列（第4、5列）有有效值
        if record.len() > 4 && (record[4].is_empty() || record[4] == "NaN") {
            record[4] = "0.000000".to_string();
        }
        if record.len() > 5 && (record[5].is_empty() || record[5] == "NaN") {
            record[5] = "0.000000".to_string();
        }
    }
    // 输出
    let mut out = File::create(&args.output)?;
    writeln!(out, "{}", header)?;
    for rec in records {
        writeln!(out, "{}", rec.join(","))?;
    }
    let elapsed = start_time.elapsed();
    
    // 覆盖进度显示，显示输出信息
    println!("\r[Output]                           ");
    println!("    Normalized: {}", args.output);
    println!("{}", crate::progress::format_time_used(elapsed));
    
    // 记录完成日志
    logger.log(&format!("Normalization processing completed, output file: {}", args.output))?;
    logger.log(&format!("Total time: {:.2}s", elapsed.as_secs_f64()))?;
    
    Ok(())
}

/// 动态滑窗归一化
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
    
    // 计算总数据中反应性碱基的数量
    let total_reactive_bases = data.iter()
        .filter(|d| reactive_bases.contains(d.base))
        .count();
    
    // 如果总反应性碱基数太少，使用全局归一化
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
    
    // 计算每个窗口的目标大小（基于总反应性碱基数）
    let window_count = (total_reactive_bases + target_reactive_count - 1) / target_reactive_count;
    let bases_per_window = (total_reactive_bases + window_count - 1) / window_count;
    
    let mut i = 0;
    let mut window_start = 0;
    
    while i < data.len() {
        // 寻找包含足够反应性碱基的窗口
        let window = find_adaptive_window(data, i, bases_per_window, reactive_bases);
        
        if window.is_empty() {
            // 如果找不到合适的窗口，使用剩余所有数据
            let remaining_data: Vec<f64> = data[i..].iter().map(|d| d.reactivity).collect();
            let normed = normalize_window(&remaining_data, norm_method);
            for (j, &value) in normed.iter().enumerate() {
                if i + j < result.len() {
                    result[i + j] = value;
                }
            }
            break;
        }
        
        // 对窗口内的数据进行归一化
        let window_reactivities: Vec<f64> = window.iter().map(|d| d.reactivity).collect();
        
        // 检查窗口内是否有足够的非零值
        let non_zero_count = window_reactivities.iter().filter(|&&x| x > 0.0).count();
        
        if non_zero_count < 3 {
            // 如果窗口内非零值太少，使用当前区域数据的全局归一化
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
                // 如果整个区域都没有非零值，保持原值
                for (j, data_item) in data.iter().enumerate() {
                    if j < result.len() {
                        result[j] = data_item.reactivity;
                    }
                }
            }
            break;
        }
        
        let normed = normalize_window(&window_reactivities, norm_method);
        
        // 将归一化结果写回对应位置
        for (j, &value) in normed.iter().enumerate() {
            if j < window.len() && window[j].index < result.len() {
                result[window[j].index] = value;
            }
        }
        
        // 移动到下一个窗口的起始位置
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

/// 寻找自适应窗口
fn find_adaptive_window(
    data: &[ReactivityData],
    start_idx: usize,
    target_reactive_count: usize,
    reactive_bases: &str,
) -> Vec<ReactivityData> {
    let mut window = Vec::new();
    let mut reactive_count = 0;
    let mut i = start_idx;
    
    // 如果目标数量为0，返回一个小的固定窗口
    if target_reactive_count == 0 {
        let end_idx = (start_idx + 100).min(data.len());
        return data[start_idx..end_idx].to_vec();
    }
    
    // 计算最大窗口大小（避免窗口过大）
    let max_window_size = (target_reactive_count * 3).min(data.len() - start_idx);
    
    while i < data.len() && reactive_count < target_reactive_count && window.len() < max_window_size {
        let item = &data[i];
        if reactive_bases.contains(item.base) {
            reactive_count += 1;
        }
        window.push(item.clone());
        i += 1;
    }
    
    // 如果找到了足够的反应性碱基，返回窗口
    if reactive_count >= target_reactive_count {
        window
    } else {
        // 如果到达序列末尾仍不够，返回所有剩余数据
        data[start_idx..].to_vec()
    }
}

/// 寻找动态窗口（保留原函数以兼容）
fn find_dynamic_window(
    data: &[ReactivityData],
    start_idx: usize,
    target_count: usize,
    reactive_bases: &str,
) -> Vec<ReactivityData> {
    find_adaptive_window(data, start_idx, target_count, reactive_bases)
}

/// 对窗口数据进行归一化
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
    if n < 10 { return values.to_vec(); }
    let p2 = (n as f64 * 0.02).ceil() as usize;
    let p10 = (n as f64 * 0.10).ceil() as usize;
    let p2 = p2.min(n-1);
    let p10 = p10.min(n-1);
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
    if n < 2 { return values.to_vec(); }
    let p5 = (n as f64 * 0.05).floor() as usize;
    let p95 = (n as f64 * 0.95).ceil() as usize;
    let v5 = sorted[p5.min(n-1)];
    let v95 = sorted[p95.min(n-1)];
    let winsorized: Vec<f64> = values.iter().map(|&x| {
        if x < v5 { v5 } else if x > v95 { v95 } else { x }
    }).collect();
    let norm_val = v95.max(1e-6);
    winsorized.iter().map(|&x| x / norm_val).collect()
}

// Box-plot Normalization
fn boxplot_norm(values: &[f64]) -> Vec<f64> {
    let mut sorted: Vec<f64> = values.iter().cloned().collect();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let n = sorted.len();
    if n < 4 { return values.to_vec(); }
    let q1_idx = (n as f64 * 0.25).ceil() as usize;
    let q3_idx = (n as f64 * 0.75).ceil() as usize;
    let q1 = sorted[q1_idx.min(n-1)];
    let q3 = sorted[q3_idx.min(n-1)];
    let iqr = q3 - q1;
    let upper = q3 + 1.5 * iqr;
    // 去除极端值
    let filtered: Vec<f64> = sorted.iter().cloned().filter(|&x| x <= upper).collect();
    let top10 = (filtered.len() as f64 * 0.10).ceil() as usize;
    let avg = if top10 > 0 {
        let sum = filtered.iter().take(top10).sum::<f64>();
        if sum > 0.0 {
            sum / top10 as f64
        } else {
            1.0
        }
    } else { 1.0 };
    values.iter().map(|&x| x / avg).collect()
}
