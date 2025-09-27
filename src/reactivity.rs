use std::error::Error;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use clap::Args;
use std::collections::HashMap;
use rayon::prelude::*;
use std::time::Instant;
use std::path::Path;
use chrono;

/// Validate reactivity command parameters
fn validate_reactivity_args(args: &ReactivityArgs) -> Result<(), Box<dyn Error>> {
    // Validate modified sample file
    if args.mod_csv.trim().is_empty() {
        return Err("Error: Modified sample file path cannot be empty".into());
    }
    if !Path::new(&args.mod_csv).exists() {
        return Err(format!("Error: Modified sample file does not exist: {}", args.mod_csv).into());
    }
    if !args.mod_csv.ends_with(".csv") {
        return Err(format!("Error: Modified sample file path must end with .csv: {}", args.mod_csv).into());
    }
    
    // Validate unmodified sample file
    if args.unmod_csv.trim().is_empty() {
        return Err("Error: Unmodified sample file path cannot be empty".into());
    }
    if !Path::new(&args.unmod_csv).exists() {
        return Err(format!("Error: Unmodified sample file does not exist: {}", args.unmod_csv).into());
    }
    if !args.unmod_csv.ends_with(".csv") {
        return Err(format!("Error: Unmodified sample file path must end with .csv: {}", args.unmod_csv).into());
    }
    
    // 验证输出文件
    if args.output.trim().is_empty() {
        return Err("错误: reactivity输出文件路径不能为空".into());
    }
    if !args.output.ends_with(".csv") {
        return Err(format!("错误: reactivity输出文件路径必须以.csv结尾: {}", args.output).into());
    }
    
    // 验证stop方法
    match args.stop_method.to_lowercase().as_str() {
        "current" | "ding" | "rouskin" => {},
        _ => return Err(format!("错误: 未知的stop信号方法: {}。支持的方法: current, ding, rouskin", args.stop_method).into()),
    }
    
    // 验证mutation方法
    match args.mutation_method.to_lowercase().as_str() {
        "current" | "siegfried" | "zubradt" => {},
        _ => return Err(format!("错误: 未知的mutation信号方法: {}。支持的方法: current, siegfried, zubradt", args.mutation_method).into()),
    }
    
    // 验证Ding方法相关参数
    if args.stop_method.to_lowercase() == "ding" {
        if args.pseudocount <= 0.0 {
            return Err("错误: Ding方法的伪计数参数必须大于0".into());
        }
        if args.maxscore <= 0.0 {
            return Err("错误: Ding方法的最大值限制参数必须大于0".into());
        }
        if args.maxscore > 100.0 {
            return Err(format!("错误: Ding方法的最大值限制参数不能超过100，当前: {}", args.maxscore).into());
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
    
    Ok(())
}

#[derive(Args)]
pub struct ReactivityArgs {
    // Input files
    /// Modified sample CSV
    #[arg(short = 'M', long = "mod")]
    pub mod_csv: String,
    /// Unmodified sample CSV
    #[arg(short = 'U', long = "unmod")]
    pub unmod_csv: String,
    
    // Output files
    /// reactivity output file (containing stop and mutation signals)
    #[arg(short = 'O', long = "output")]
    pub output: String,
    
    // Method selection
    /// stop signal reactivity method: current, ding, rouskin
    #[arg(short = 's', long = "stop-method", default_value = "current")]
    pub stop_method: String,
    /// mutation signal reactivity method: current, siegfried, zubradt
    #[arg(short = 'm', long = "mutation-method", default_value = "current")]
    pub mutation_method: String,
    
    // Method-specific parameters (only used with specific methods)
    /// Pseudocount parameter to avoid zero-value issues in logarithmic calculations (Ding method)
    #[arg(long = "pseudocount", default_value_t = 1.0)]
    pub pseudocount: f64,
    /// Maximum score limit parameter to limit the upper bound of reactivity values (Ding method)
    #[arg(long = "maxscore", default_value_t = 10.0)]
    pub maxscore: f64,
    
    // Performance configuration
    /// Number of parallel threads
    #[arg(short = 't', long = "threads")]
    pub threads: Option<usize>,
    
    // Logging configuration
    /// Log file path (optional)
    #[arg(short = 'l', long = "log")]
    pub log: Option<String>,
}



#[derive(Debug, Clone, PartialEq)]
pub enum StopReactivityMethod {
    Current,  // 现有方法：差分
    Ding,     // Ding方法
    Rouskin,  // Rouskin方法：只用mod样本stop信号
}

impl std::str::FromStr for StopReactivityMethod {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "current" => Ok(StopReactivityMethod::Current),
            "ding" => Ok(StopReactivityMethod::Ding),
            "rouskin" => Ok(StopReactivityMethod::Rouskin),
            _ => Err(format!("未知的stop信号reactivity方法: {}", s)),
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub enum MutationReactivityMethod {
    Current,   // 现有方法
    Siegfried, // Siegfried方法
    Zubradt,   // Zubradt方法：只用mod样本mutation信号
}

impl std::str::FromStr for MutationReactivityMethod {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "current" => Ok(MutationReactivityMethod::Current),
            "siegfried" => Ok(MutationReactivityMethod::Siegfried),
            "zubradt" => Ok(MutationReactivityMethod::Zubradt),
            _ => Err(format!("未知的mutation信号reactivity方法: {}", s)),
        }
    }
}



/// 计算带k因子校正的reactivity
/// 公式：Rᵢ = max(0, Tᵢ - k * Uᵢ)
/// 其中 k 是基于背景区信号强度比值估算的校正因子
fn calculate_k_corrected_reactivity(
    mod_signal: f64,
    unmod_signal: f64,
    k_factor: f64,
) -> f64 {
    (mod_signal - k_factor * unmod_signal).max(0.0)
}

/// Ding et al., 2014 方法的reactivity计算
/// 正确的Ding方法：对数归一化后直接计算差异
/// 公式：S_i = max(0, T_i - U_i)
/// 其中 T_i = ln(n_Ti + p) / mean(ln(n_Tj + p)), U_i = ln(n_Ui + p) / mean(ln(n_Uj + p))
fn calculate_reactivity_ding_style(
    mod_signal: f64,
    unmod_signal: f64,
    _k_factor: f64, // 不使用k因子校正
    pseudocount: f64,
    _maxscore: f64, // 不使用maxscore限制
) -> f64 {
    // 对原始信号进行对数处理（加上伪计数避免ln(0)）
    let log_mod = (mod_signal + pseudocount).ln();
    let log_unmod = (unmod_signal + pseudocount).ln();
    
    // 直接计算差异（不使用k因子校正）
    let reactivity = log_mod - log_unmod;
    
    // 返回最大值0
    reactivity.max(0.0)
}

fn parse_signal(fields: &Vec<String>, rfc_idx: usize, pipc_idx: usize, depth_idx: usize) -> (f64, f64) {
    let depth: f64 = fields[depth_idx].parse().unwrap_or(0.0);
    let rfc: f64 = fields[rfc_idx].parse().unwrap_or(0.0);
    let pipc: f64 = fields[pipc_idx].parse().unwrap_or(0.0);
    let mutation = if depth > 0.0 { rfc / depth } else { 0.0 };
    let stop = if depth > 0.0 { pipc / depth } else { 0.0 };
    (mutation, stop)
}

pub fn calc_reactivity(args: &ReactivityArgs, logger: &mut crate::Logger) -> Result<(), Box<dyn Error>> {
    // 验证reactivity命令参数
    validate_reactivity_args(args)?;
    let start_time = Instant::now();
    

    
    // 记录环境信息和参数
    logger.log("=== ModDetector Reactivity Function Log ===")?;
    logger.log(&format!("Software Version: v{}", crate::VERSION))?;
    logger.log(&format!("Runtime: {}", chrono::Utc::now().format("%Y-%m-%d %H:%M:%S")))?;
    logger.log(&format!("Modified Sample File: {}", args.mod_csv))?;
    logger.log(&format!("Unmodified Sample File: {}", args.unmod_csv))?;
    logger.log(&format!("Reactivity Output: {}", args.output))?;
    logger.log(&format!("Stop Method: {}", args.stop_method))?;
    logger.log(&format!("Mutation Method: {}", args.mutation_method))?;
    logger.log(&format!("Pseudocount Parameter: {}", args.pseudocount))?;
    logger.log(&format!("Max Score Limit: {}", args.maxscore))?;
    if let Some(threads) = args.threads {
        logger.log(&format!("Threads: {}", threads))?;
    } else {
        logger.log("Threads: Auto-detect")?;
    }
    logger.log("Starting reactivity calculation...")?;
    
    // 解析stop信号方法
    let stop_method: StopReactivityMethod = args.stop_method.parse().unwrap_or(StopReactivityMethod::Current);
    
    // 解析mutation信号方法
    let mutation_method: MutationReactivityMethod = args.mutation_method.parse().unwrap_or(MutationReactivityMethod::Current);
    
    // 设置线程数
    let thread_count = args.threads.unwrap_or(8);
    rayon::ThreadPoolBuilder::new().num_threads(thread_count).build_global().ok();
    
    println!("[Loading data] Threads: {}", thread_count);
    println!("    Mod={}", args.mod_csv);
    println!("    Unmod={}", args.unmod_csv);
    println!();
    
    // 使用分块并行数据加载
    let (mod_records, unmod_records, _mod_key_map, unmod_key_map, _header) = 
        load_data_in_chunks(&args.mod_csv, &args.unmod_csv, thread_count)?;
    
    // 1. 并行收集背景区信号，计算k因子
    let mut stop_mod_bg = Vec::new();
    let mut stop_unmod_bg = Vec::new();
    let mut mut_mod_bg = Vec::new();
    let mut mut_unmod_bg = Vec::new();
    let mut all_pairs: Vec<(usize, f64, f64, f64, f64)> = Vec::new(); // (mod_idx, m_mut, m_stop, u_mut, u_stop)
    
    let total_mod = mod_records.len();
    let mut progress = crate::progress::SimpleProgress::new(total_mod);
    
    // 显示处理进度
    print!("\r[Processing] ");
    std::io::stdout().flush()?;
    
    // 并行处理数据收集
    let chunk_size = (total_mod + thread_count - 1) / thread_count;
    let chunks: Vec<Vec<usize>> = (0..total_mod)
        .collect::<Vec<usize>>()
        .chunks(chunk_size)
        .map(|chunk| chunk.to_vec())
        .collect();
    
    let results: Vec<(Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>, Vec<(usize, f64, f64, f64, f64)>)> = chunks
        .into_par_iter()
        .map(|chunk| {
            let mut local_stop_mod_bg = Vec::new();
            let mut local_stop_unmod_bg = Vec::new();
            let mut local_mut_mod_bg = Vec::new();
            let mut local_mut_unmod_bg = Vec::new();
            let mut local_pairs = Vec::new();
            
            for &i in &chunk {
                let m = &mod_records[i];
                let key = format!("{}|{}|{}", m[0], m[1], m[2]);
                if let Some(&u_idx) = unmod_key_map.get(&key) {
                    let u = &unmod_records[u_idx];
                    let (m_mut, m_stop) = parse_signal(m, 4, 5, 6);
                    let (u_mut, u_stop) = parse_signal(u, 4, 5, 6);
                    let is_bg = m[0].contains("rRNA") || m[0].contains("16S") || m[0].contains("18S") || args.mod_csv.contains("16S") || args.mod_csv.contains("18S");
                    if is_bg {
                        local_stop_mod_bg.push(m_stop);
                        local_stop_unmod_bg.push(u_stop);
                        local_mut_mod_bg.push(m_mut);
                        local_mut_unmod_bg.push(u_mut);
                    }
                    local_pairs.push((i, m_mut, m_stop, u_mut, u_stop));
                }
            }
            
            (local_stop_mod_bg, local_stop_unmod_bg, local_mut_mod_bg, local_mut_unmod_bg, local_pairs)
        })
        .collect();
    
    // 合并结果
    for (local_stop_mod, local_stop_unmod, local_mut_mod, local_mut_unmod, local_pairs) in results {
        stop_mod_bg.extend(local_stop_mod);
        stop_unmod_bg.extend(local_stop_unmod);
        mut_mod_bg.extend(local_mut_mod);
        mut_unmod_bg.extend(local_mut_unmod);
        all_pairs.extend(local_pairs);
    }
    
    progress.finish()?;
    println!("\r[Data info]                      ");
    println!("    mod records: {}, unmod records: {}", mod_records.len(), unmod_records.len());
    println!("    matched positions: {}, background regions: {}", all_pairs.len(), stop_mod_bg.len());
    
    // 计算k因子
    print!("\r[Calculating] Computing k coefficients...");
    std::io::stdout().flush()?;
    let mut k_stop = 1.0;
    let mut k_mut = 1.0;
    
    if !stop_mod_bg.is_empty() && !stop_unmod_bg.is_empty() {
        let mod_mean = stop_mod_bg.iter().sum::<f64>() / stop_mod_bg.len() as f64;
        let unmod_mean = stop_unmod_bg.iter().sum::<f64>() / stop_unmod_bg.len() as f64;
        
        if unmod_mean > 0.0 {
            k_stop = mod_mean / unmod_mean;
        }
    }
    
    if !mut_mod_bg.is_empty() && !mut_unmod_bg.is_empty() {
        let mod_mean = mut_mod_bg.iter().sum::<f64>() / mut_mod_bg.len() as f64;
        let unmod_mean = mut_unmod_bg.iter().sum::<f64>() / mut_unmod_bg.len() as f64;
        
        if unmod_mean > 0.0 {
            k_mut = mod_mean / unmod_mean;
        }
    }
    
    println!("\r[k value]                                ");
    println!("    stop signal mean: mod={:.4}, unmod={:.4}. k_stop = {:.3}", 
             stop_mod_bg.iter().sum::<f64>() / stop_mod_bg.len().max(1) as f64,
             stop_unmod_bg.iter().sum::<f64>() / stop_unmod_bg.len().max(1) as f64,
             k_stop);
    println!("    mutation signal mean: mod={:.4}, unmod={:.4}. k_mut = {:.3}", 
             mut_mod_bg.iter().sum::<f64>() / mut_mod_bg.len().max(1) as f64,
             mut_unmod_bg.iter().sum::<f64>() / mut_unmod_bg.len().max(1) as f64,
             k_mut);
    
    // 2. 并行计算reactivity
    let total_pairs = all_pairs.len();
    let mut calc_progress = crate::progress::SimpleProgress::new(total_pairs);
    
    print!("\r[Reactivity calculating] ");
    std::io::stdout().flush()?;
    
    // 并行计算reactivity
    let reactivity_chunk_size = (total_pairs + thread_count - 1) / thread_count;
    let reactivity_chunks: Vec<Vec<(usize, f64, f64, f64, f64)>> = all_pairs
        .chunks(reactivity_chunk_size)
        .map(|chunk| chunk.to_vec())
        .collect();
    
    let reactivity_results: Vec<Vec<(f64, f64)>> = reactivity_chunks
        .into_par_iter()
        .map(|chunk| {
            let mut local_results = Vec::with_capacity(chunk.len());
            
            for (_mod_idx, m_mut, m_stop, u_mut, u_stop) in chunk {
                // stop信号方法选择
                let stop_reactivity = match stop_method {
                    StopReactivityMethod::Current => calculate_k_corrected_reactivity(m_stop, u_stop, k_stop),
                    StopReactivityMethod::Ding => calculate_reactivity_ding_style(m_stop, u_stop, k_stop, args.pseudocount, args.maxscore),
                    StopReactivityMethod::Rouskin => m_stop, // 只用mod样本stop信号
                };
                
                // mutation信号方法选择
                let mut_reactivity = match mutation_method {
                    MutationReactivityMethod::Current => calculate_k_corrected_reactivity(m_mut, u_mut, k_mut),
                    MutationReactivityMethod::Siegfried => {
                        // Siegfried方法：S_i = (n_Ti / c_Ti) - (n_Ui / c_Ui)
                        // 即直接计算mutation rate的差值，不使用k因子校正
                        (m_mut - u_mut).max(0.0)
                    },
                    MutationReactivityMethod::Zubradt => m_mut, // 只用mod样本mutation信号
                };
                
                local_results.push((stop_reactivity, mut_reactivity));
            }
            
            local_results
        })
        .collect();
    
    calc_progress.finish()?;
    
    // 3. 并行生成输出文件
    let mut output_file = File::create(&args.output)?;
    
    // 写入表头
    let header_line = "ChrID,Strand,Position,Base,reactivity_stop,reactivity_mutation".to_string();
    writeln!(output_file, "{}", header_line)?;
    
    // 并行准备输出数据
    let write_chunk_size = (total_pairs + thread_count - 1) / thread_count;
    let write_chunks: Vec<Vec<usize>> = (0..total_pairs)
        .collect::<Vec<usize>>()
        .chunks(write_chunk_size)
        .map(|chunk| chunk.to_vec())
        .collect();
    
    let prepared_outputs: Vec<Vec<String>> = write_chunks
        .into_par_iter()
        .map(|chunk| {
            let mut lines = Vec::with_capacity(chunk.len());
            
            for &i in &chunk {
                let (mod_idx, _m_mut, _m_stop, _u_mut, _u_stop) = all_pairs[i];
                let rec = &mod_records[mod_idx];
                
                // 获取并行计算的reactivity结果
                let chunk_index = i / reactivity_chunk_size;
                let local_index = i % reactivity_chunk_size;
                let (stop_reactivity, mut_reactivity) = reactivity_results[chunk_index][local_index];
                
                // 简化输出格式：只保留核心列
                let chr_id = &rec[0];
                let strand = &rec[1];
                let position = &rec[2];
                let base = &rec[3];
                
                let line = format!("{},{},{},{},{:.6},{:.6}", 
                    chr_id, strand, position, base, stop_reactivity, mut_reactivity);
                lines.push(line);
            }
            
            lines
        })
        .collect();
    
    // 并行写入文件
    let mut write_progress = crate::progress::SimpleProgress::new(total_pairs);
    print!("\r[Writing files] ");
    std::io::stdout().flush()?;
    
    // 并行写入策略：将输出分成多个临时文件，然后合并
    let temp_dir = std::env::temp_dir().join("moddetector_reactivity");
    std::fs::create_dir_all(&temp_dir)?;
    
    let write_results: Vec<Vec<String>> = prepared_outputs
        .into_par_iter()
        .enumerate()
        .map(|(chunk_id, lines)| {
            let mut chunk_output = Vec::new();
            
            // 创建临时文件
            let temp_file = temp_dir.join(format!("reactivity_chunk_{}.tmp", chunk_id));
            
            // 写入临时文件
            let mut temp_out = File::create(&temp_file).unwrap();
            
            for line in lines.into_iter() {
                writeln!(temp_out, "{}", line).unwrap();
                chunk_output.push(line);
            }
            
            // 关闭临时文件
            drop(temp_out);
            
            chunk_output
        })
        .collect();
    
    // 合并所有临时文件到最终输出
    let mut processed_count = 0;
    for (chunk_id, chunk_output) in write_results.into_iter().enumerate() {
        // 读取临时文件并写入最终输出
        let temp_file = temp_dir.join(format!("reactivity_chunk_{}.tmp", chunk_id));
        
        if temp_file.exists() {
            let content = std::fs::read_to_string(&temp_file)?;
            
            output_file.write_all(content.as_bytes())?;
            
            // 删除临时文件
            std::fs::remove_file(temp_file).ok();
            
            processed_count += chunk_output.len();
            if processed_count % 10000 == 0 {
                write_progress.update(processed_count)?;
            }
        }
    }
    
    // 清理临时目录
    std::fs::remove_dir_all(temp_dir).ok();
    
    write_progress.finish()?;
    
    let elapsed = start_time.elapsed();
    
    // 根据实际使用的方法显示输出信息
    let stop_method_name = match stop_method {
        StopReactivityMethod::Current => "Current",
        StopReactivityMethod::Ding => "Ding",
        StopReactivityMethod::Rouskin => "Rouskin",
    };
    
    let mutation_method_name = match mutation_method {
        MutationReactivityMethod::Current => "Current",
        MutationReactivityMethod::Siegfried => "Siegfried",
        MutationReactivityMethod::Zubradt => "Zubradt",
    };
    
    println!("\r[Reactivity outputs] stop={}, mutation={}", stop_method_name, mutation_method_name);
    println!("    output={}", args.output);
    println!("{}", crate::progress::format_time_used(elapsed));
    
    // 记录完成日志
    logger.log(&format!("Reactivity calculation completed"))?;
    logger.log(&format!("Reactivity Output: {}", args.output))?;
    logger.log(&format!("Stop Method: {}", stop_method_name))?;
    logger.log(&format!("Mutation Method: {}", mutation_method_name))?;
    logger.log(&format!("Total time: {:.2}s", elapsed.as_secs_f64()))?;
    
    Ok(())
} 

/// 分块并行加载CSV数据（原有实现，作为备用）
fn load_data_in_chunks(mod_csv: &str, unmod_csv: &str, thread_count: usize) -> Result<(Vec<Vec<String>>, Vec<Vec<String>>, HashMap<String, usize>, HashMap<String, usize>, String), Box<dyn Error>> {
    use std::sync::atomic::{AtomicUsize, Ordering};
    use std::sync::Arc;
    
    println!("[Loading] Starting parallel data loading...");
    
    // 读取mod文件
    let mod_file = File::open(mod_csv)?;
    let mut mod_lines: Vec<String> = BufReader::new(mod_file).lines().collect::<Result<Vec<String>, _>>()?;
    let header = mod_lines.remove(0);
    
    // 读取unmod文件
    let unmod_file = File::open(unmod_csv)?;
    let mut unmod_lines: Vec<String> = BufReader::new(unmod_file).lines().collect::<Result<Vec<String>, _>>()?;
    let _ = unmod_lines.remove(0); // 跳过unmod表头
    
    println!("[Loading] Mod records: {}, Unmod records: {}", mod_lines.len(), unmod_lines.len());
    
    // 分块处理mod数据
    let mod_chunk_size = (mod_lines.len() + thread_count - 1) / thread_count;
    let mod_chunks: Vec<Vec<String>> = mod_lines.chunks(mod_chunk_size).map(|chunk| chunk.to_vec()).collect();
    
    // 分块处理unmod数据
    let unmod_chunk_size = (unmod_lines.len() + thread_count - 1) / thread_count;
    let unmod_chunks: Vec<Vec<String>> = unmod_lines.chunks(unmod_chunk_size).map(|chunk| chunk.to_vec()).collect();
    
    let num_chunks = mod_chunks.len() + unmod_chunks.len();
    let completed = Arc::new(AtomicUsize::new(0));
    
    print!("\r[Loading] Processing {} chunks with {} threads", num_chunks, thread_count);
    std::io::stdout().flush().ok();
    
    // 并行处理mod数据
    let mod_results: Vec<(Vec<Vec<String>>, HashMap<String, usize>)> = mod_chunks
        .into_par_iter()
        .map(|chunk| {
            let mut records = Vec::with_capacity(chunk.len());
            let mut key_map = HashMap::with_capacity(chunk.len());
            
            for (i, line) in chunk.iter().enumerate() {
                let fields: Vec<String> = line.split(',').map(|s| s.to_string()).collect();
                if fields.len() >= 3 {
                    let key = format!("{}|{}|{}", fields[0], fields[1], fields[2]);
                    key_map.insert(key, i);
                    records.push(fields);
                }
            }
            
            // 更新进度
            let completed_count = completed.fetch_add(1, Ordering::Relaxed) + 1;
            let percentage = (completed_count as f64 * 100.0) / num_chunks as f64;
            print!("\r[Loading] Processing {}/{} chunks ({:.1}%)", completed_count, num_chunks, percentage);
            std::io::stdout().flush().ok();
            
            (records, key_map)
        })
        .collect();
    
    // 并行处理unmod数据
    let unmod_results: Vec<(Vec<Vec<String>>, HashMap<String, usize>)> = unmod_chunks
        .into_par_iter()
        .map(|chunk| {
            let mut records = Vec::with_capacity(chunk.len());
            let mut key_map = HashMap::with_capacity(chunk.len());
            
            for (i, line) in chunk.iter().enumerate() {
                let fields: Vec<String> = line.split(',').map(|s| s.to_string()).collect();
                if fields.len() >= 3 {
                    let key = format!("{}|{}|{}", fields[0], fields[1], fields[2]);
                    key_map.insert(key, i);
                    records.push(fields);
                }
            }
            
            // 更新进度
            let completed_count = completed.fetch_add(1, Ordering::Relaxed) + 1;
            let percentage = (completed_count as f64 * 100.0) / num_chunks as f64;
            print!("\r[Loading] Processing {}/{} chunks ({:.1}%)", completed_count, num_chunks, percentage);
            std::io::stdout().flush().ok();
            
            (records, key_map)
        })
        .collect();
    
    println!("\r[Loading] Data loading completed");
    
    // 合并结果
    let mut mod_records = Vec::with_capacity(mod_lines.len());
    let mut mod_key_map = HashMap::with_capacity(mod_lines.len());
    let mut offset = 0;
    
    for (records, key_map) in mod_results {
        for (key, local_idx) in key_map {
            mod_key_map.insert(key, offset + local_idx);
        }
        let records_len = records.len();
        mod_records.extend(records);
        offset += records_len;
    }
    
    let mut unmod_records = Vec::with_capacity(unmod_lines.len());
    let mut unmod_key_map = HashMap::with_capacity(unmod_lines.len());
    let mut offset = 0;
    
    for (records, key_map) in unmod_results {
        for (key, local_idx) in key_map {
            unmod_key_map.insert(key, offset + local_idx);
        }
        let records_len = records.len();
        unmod_records.extend(records);
        offset += records_len;
    }
    
    println!("[Loading] Loaded {} mod records, {} unmod records", mod_records.len(), unmod_records.len());
    
    Ok((mod_records, unmod_records, mod_key_map, unmod_key_map, header))
}