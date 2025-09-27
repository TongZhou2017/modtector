use std::collections::HashMap;
use std::error::Error;
use std::fs::File;
use std::io::{BufRead, BufReader, Write, BufWriter};
use plotters::prelude::*;
use rayon::prelude::*;
use chrono;
use regex::Regex;

#[derive(Debug, Clone)]
struct SignalData {
    position: u32,
    stop_signal: f64,    // PIPC signal
    mutation_signal: f64, // RFC signal
    depth: usize,
}

#[derive(Debug, Clone)]
pub struct ReactivityData {
    pub position: u32,
    pub reactivity_stop: f64,
    pub reactivity_mutation: f64,
}

#[derive(Debug, Clone)]
struct GenePlotTask {
    gene_id: String,
    mod_signals: Vec<SignalData>,
    unmod_signals: Vec<SignalData>,
    reactivity_data: Option<Vec<ReactivityData>>,
    output_dir: String,
    is_high: bool,
}

#[derive(Debug, Clone)]
struct GeneStat {
    gene_id: String,
    coverage: f64,
    avg_depth: f64,
    cov_x_depth: f64,
}

fn load_signal_data(csv_file: &str) -> Result<HashMap<String, Vec<SignalData>>, Box<dyn Error>> {
    let mut gene_data: HashMap<String, Vec<SignalData>> = HashMap::new();
    let file = File::open(csv_file)?;
    let reader = BufReader::new(file);
    
    for (line_num, line) in reader.lines().enumerate() {
        if line_num == 0 {
            continue; // Skip header row
        }
        
        let line = line?;
        let fields: Vec<&str> = line.split(',').collect();
        
        if fields.len() >= 13 {
            let chr_id = fields[0].trim_matches('"').to_string();
            let strand = fields[1].trim_matches('"').to_string();
            let position = fields[2].parse::<u32>()?;
            let rfc = fields[4].parse::<usize>()?;
            // fields[5] is pipe_truncation_count (float), needs to be converted to integer
            let pipc = fields[5].parse::<f64>()?.round() as usize;
            let depth = fields[6].parse::<usize>()?;
            
            // Calculate signal strength (normalized)
            let stop_signal = if depth > 0 {
                pipc as f64 / depth as f64
            } else {
                0.0
            };
            
            let mutation_signal = if depth > 0 {
                rfc as f64 / depth as f64
            } else {
                0.0
            };
            
            let signal_data = SignalData {
                position,
                stop_signal,
                mutation_signal,
                depth,
            };
            
            // Use chr_id+strand as gene_id to ensure separate processing of positive and negative strands
            let gene_id = format!("{}_{}", chr_id, strand);
            gene_data.entry(gene_id).or_insert_with(Vec::new).push(signal_data);
        }
    }
    
    // Sort data by position for each gene
    for data in gene_data.values_mut() {
        data.sort_by_key(|d| d.position);
    }
    
    Ok(gene_data)
}

pub fn load_reactivity_data(csv_file: &str) -> Result<HashMap<String, Vec<ReactivityData>>, Box<dyn Error>> {
    let mut gene_data: HashMap<String, Vec<ReactivityData>> = HashMap::new();
    let file = File::open(csv_file)?;
    let reader = BufReader::new(file);
    
    // Determine file type (stop or mutation)
    let is_stop_file = csv_file.contains("stop");
    let is_mutation_file = csv_file.contains("mutation");
    
    for (line_num, line) in reader.lines().enumerate() {
        if line_num == 0 {
            continue; // Skip header row
        }
        
        let line = line?;
        let fields: Vec<&str> = line.split(',').collect();
        
        if fields.len() >= 6 { // 包含reactivity列
            let chr_id = fields[0].trim_matches('"').to_string();
            let strand = fields[1].trim_matches('"').to_string();
            let position = fields[2].parse::<u32>()?;
            
            // 简化格式：第4列是reactivity_stop，第5列是reactivity_mutation
            let stop_val = fields[4].parse::<f64>()?;
            let mut_val = fields[5].parse::<f64>()?;
            let (reactivity_stop, reactivity_mutation) = (stop_val, mut_val);
            
            let reactivity_data = ReactivityData {
                position,
                reactivity_stop,
                reactivity_mutation,
            };
            
            // Use chr_id+strand as gene_id to ensure separate processing of positive and negative strands
            let gene_id = format!("{}_{}", chr_id, strand);
            gene_data.entry(gene_id).or_insert_with(Vec::new).push(reactivity_data);
        }
    }
    
    // Sort data by position for each gene
    for data in gene_data.values_mut() {
        data.sort_by_key(|d| d.position);
    }
    
    Ok(gene_data)
}

fn gene_stats(gene_data: &HashMap<String, Vec<SignalData>>) -> HashMap<String, GeneStat> {
    let mut stats = HashMap::new();
    for (gene_id, signals) in gene_data {
        let n = signals.len();
        if n == 0 { continue; }
        let covered = signals.iter().filter(|d| d.depth > 0).count();
        let avg_depth = signals.iter().map(|d| d.depth as f64).sum::<f64>() / n as f64;
        let coverage = covered as f64 / n as f64;
        let cov_x_depth = coverage * avg_depth;
        stats.insert(gene_id.clone(), GeneStat {
            gene_id: gene_id.clone(),
            coverage,
            avg_depth,
            cov_x_depth,
        });
    }
    stats
}

/// 解析GFF/GTF文件，返回gene_id->strand的映射
/// 支持GFF3格式：locus_tag=XXX 或 ID=XXX
/// 支持GTF格式：gene_id "XXX" 或 transcript_id "XXX"
pub fn parse_gff_strand_map(gff_file: &str) -> Result<HashMap<String, String>, Box<dyn Error>> {
    let mut map = HashMap::new();
    let file = File::open(gff_file)?;
    let reader = BufReader::new(file);
    
    // 检测文件格式
    let mut is_gtf = false;
    let mut first_data_line = true;
    
    for line in reader.lines() {
        let line = line?;
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 9 { continue; }
        
        // 检测格式（基于第一行数据）
        if first_data_line {
            let attr = fields[8];
            // GTF格式：gene_id "ENSG00000243485"
            // GFF3格式：ID=BOGIMF_00001
            is_gtf = attr.contains("gene_id \"") || attr.contains("transcript_id \"");
            first_data_line = false;
        }
        
        let strand = fields[6].trim();
        let attr = fields[8];
        
        // 根据格式提取gene_id
        let mut gene_id = None;
        
        if is_gtf {
            // GTF格式解析 - 优先使用transcript_id，因为norm文件使用转录本ID
            for kv in attr.split(';') {
                let kv = kv.trim();
                if kv.starts_with("transcript_id \"") {
                    // 优先使用transcript_id
                    if let Some(start) = kv.find('"') {
                        if let Some(end) = kv[start+1..].find('"') {
                            gene_id = Some(kv[start+1..start+1+end].to_string());
                            break;
                        }
                    }
                } else if kv.starts_with("gene_id \"") {
                    // 如果没有transcript_id，使用gene_id
                    if gene_id.is_none() {
                        if let Some(start) = kv.find('"') {
                            if let Some(end) = kv[start+1..].find('"') {
                                gene_id = Some(kv[start+1..start+1+end].to_string());
                            }
                        }
                    }
                }
            }
        } else {
            // GFF3格式解析
            for kv in attr.split(';') {
                let kv = kv.trim();
                if kv.starts_with("locus_tag=") {
                    gene_id = Some(kv.trim_start_matches("locus_tag=").to_string());
                    break;
                } else if kv.starts_with("ID=") {
                    gene_id = Some(kv.trim_start_matches("ID=").to_string());
                }
            }
        }
        
        if let Some(id) = gene_id {
            map.insert(id, strand.to_string());
        }
    }
    
    Ok(map)
}

// 修改plot_signal_distributions函数签名，增加gff_file参数
pub fn plot_signal_distributions(
    mod_file: &str,
    unmod_file: &str,
    output_dir: &str,
    num_threads: Option<usize>,
    coverage_threshold: Option<f64>,
    depth_threshold: Option<f64>,
    reactivity_file: Option<&str>,
    gff_file: Option<&str>, // 新增参数
    logger: &mut crate::Logger,
) -> Result<(), Box<dyn Error>> {
    let start_time = std::time::Instant::now();
    

    
    // 记录环境信息和参数
    logger.log("=== ModDetector Plot Function Log ===")?;
    logger.log(&format!("Software Version: v{}", crate::VERSION))?;
    logger.log(&format!("Runtime: {}", chrono::Utc::now().format("%Y-%m-%d %H:%M:%S")))?;
    logger.log(&format!("Modified Sample File: {}", mod_file))?;
    logger.log(&format!("Unmodified Sample File: {}", unmod_file))?;
    logger.log(&format!("Output Directory: {}", output_dir))?;
    if let Some(coverage) = coverage_threshold {
        logger.log(&format!("Coverage Threshold: {}", coverage))?;
    } else {
        logger.log("Coverage Threshold: Default 0.2")?;
    }
    if let Some(depth) = depth_threshold {
        logger.log(&format!("Depth Threshold: {}", depth))?;
    } else {
        logger.log("Depth Threshold: Default 50")?;
    }
    if let Some(reactivity_file) = reactivity_file {
        logger.log(&format!("Reactivity File: {}", reactivity_file))?;
    }
    if let Some(gff) = gff_file {
        logger.log(&format!("GFF Annotation File: {}", gff))?;
    }
    if let Some(threads) = num_threads {
        logger.log(&format!("Threads: {}", threads))?;
    } else {
        logger.log("Threads: Auto-detect")?;
    }
    logger.log("Starting signal distribution plotting...")?;
    
    // 显示加载数据信息
    println!("[Loading data]");
    println!("    mod: {}", mod_file);
    println!("    unmod: {}", unmod_file);
    if let Some(reactivity_file) = reactivity_file {
        println!("    reactivity: {}", reactivity_file);
    }
    if let Some(gff_file) = gff_file {
        println!("    gene annotation: {}", gff_file);
    }
    println!();
    
    // 加载数据
    let mut mod_data = load_signal_data(mod_file)?;
    let mut unmod_data = load_signal_data(unmod_file)?;
    
    // 保存未过滤的基因数量
    let total_genes_before_filter = mod_data.len();
    let mut plus_strand_before_filter = 0;
    let mut minus_strand_before_filter = 0;
    for gene_id in mod_data.keys() {
        if gene_id.ends_with('+') {
            plus_strand_before_filter += 1;
        } else if gene_id.ends_with('-') {
            minus_strand_before_filter += 1;
        }
    }
    
    // ====== 新增：GFF注释链筛选 ======
    if let Some(gff_path) = gff_file {
        let gff_map = parse_gff_strand_map(gff_path)?;
        // 只保留注释链
        let filter_keys: Vec<String> = mod_data.keys().filter(|k| {
            // k格式: geneid[_EupX_EdownY]_strand
            if let Some(pos) = k.rfind('_') {
                let mut geneid = &k[..pos];
                let strand = &k[pos+1..];
                // 去除_Eup/Edown等上下游信息
                if let Some(eup_pos) = geneid.find("_Eup") {
                    geneid = &geneid[..eup_pos];
                }
                if let Some(anno_strand) = gff_map.get(geneid) {
                    return strand == anno_strand;
                }
            }
            false
        }).cloned().collect();
        mod_data = mod_data.into_iter().filter(|(k,_)| filter_keys.contains(k)).collect();
        unmod_data = unmod_data.into_iter().filter(|(k,_)| filter_keys.contains(k)).collect();
    }
    
    // 加载reactivity数据
    let mut reactivity_combined_data: Option<HashMap<String, Vec<ReactivityData>>> = None;
    
    if let Some(reactivity_file_path) = reactivity_file {
        if std::path::Path::new(reactivity_file_path).exists() {
            let reactivity_data = load_reactivity_data(reactivity_file_path)?;
            reactivity_combined_data = Some(reactivity_data);
        }
    }
    
    // 确保输出目录存在
    std::fs::create_dir_all(output_dir)?;
    
    // 设置线程数
    let threads = num_threads.unwrap_or(rayon::current_num_threads());
    if let Some(threads) = num_threads {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build_global()
            .unwrap_or_else(|_| println!("Warning: Unable to set thread count, using default settings"));
    }
    
    // 统计基因覆盖度和深度
    let stats = gene_stats(&mod_data);
    
    // 设置coverage和depth的阈值
    let coverage_threshold = coverage_threshold.unwrap_or(0.2); // 默认0.2
    let depth_threshold = depth_threshold.unwrap_or(50.0); // 默认50 reads
    
    // 显示参数信息
    println!("[Params]");
    println!("    Threads: {}.", threads);
    println!("    Coverage threshold: {:.1} %.", coverage_threshold * 100.0);
    println!("    Depth threshold: {:.0} reads.", depth_threshold);
    println!();
    
    // 统计高低组数量
    let mut high_count = 0;
    let mut low_count = 0;
    for stat in stats.values() {
        if stat.coverage >= coverage_threshold && stat.avg_depth >= depth_threshold {
            high_count += 1;
        } else {
            low_count += 1;
        }
    }
    
    // 统计正负链数量（过滤后）
    let mut plus_strand = 0;
    let mut minus_strand = 0;
    for gene_id in mod_data.keys() {
        if gene_id.ends_with('+') {
            plus_strand += 1;
        } else if gene_id.ends_with('-') {
            minus_strand += 1;
        }
    }
    
    // 显示数据信息
    println!("[Data info]");
    println!("    Total gene regions: {}, {} (+ strand), {} (- strand).", total_genes_before_filter, plus_strand_before_filter, minus_strand_before_filter);
    if gff_file.is_some() {
        println!("    Gene annotation matched: {}, {} (+ strand), {} (- strand).", mod_data.len(), plus_strand, minus_strand);
    }
    println!("    High coverage: {}. Low coverage: {}.", high_count, low_count);
    println!();
    
    // 显示进度信息
    print!("\r[Progressing] Plotting {}/{} ({:.1}%)", 0, mod_data.len(), 0.0);
    std::io::stdout().flush()?;
    
    // 创建高/低目录
    let high_dir = format!("{}/high", output_dir);
    let low_dir = format!("{}/low", output_dir);
    std::fs::create_dir_all(&high_dir)?;
    std::fs::create_dir_all(&low_dir)?;
    
    // 准备绘图任务
    let mut plot_tasks: Vec<GenePlotTask> = Vec::new();
    for (gene_id, mod_signals) in &mod_data {
        if let Some(unmod_signals) = unmod_data.get(gene_id) {
            let stat = stats.get(gene_id).unwrap();
            let is_high = stat.coverage >= coverage_threshold && stat.avg_depth >= depth_threshold;
            
            // 获取reactivity数据
            let reactivity_data = if let Some(ref combined_data) = reactivity_combined_data {
                if let Some(reactivity) = combined_data.get(gene_id) {
                    Some(reactivity.clone())
                } else {
                    None
                }
            } else {
                None
            };
            
            plot_tasks.push(GenePlotTask {
                gene_id: gene_id.clone(),
                mod_signals: mod_signals.clone(),
                unmod_signals: unmod_signals.clone(),
                reactivity_data,
                output_dir: if is_high { high_dir.clone() } else { low_dir.clone() },
                is_high,
            });
        }
    }
    
    // 使用原子计数器跟踪进度
    use std::sync::atomic::{AtomicUsize, Ordering};
    use std::sync::Arc;
    
    let completed = Arc::new(AtomicUsize::new(0));
    let total_tasks = plot_tasks.len();
    
    // 并行处理绘图任务
    let results: Vec<Result<(), Box<dyn Error + Send + Sync>>> = plot_tasks
        .into_par_iter()
        .map(|task| {
            let gene_id = &task.gene_id;
            
            // 绘制stop信号图
            plot_stop_signals(gene_id, &task.mod_signals, &task.unmod_signals, &task.output_dir)?;
            
            // 绘制mutation信号图
            plot_mutation_signals(gene_id, &task.mod_signals, &task.unmod_signals, &task.output_dir)?;
            
            // 绘制reactivity图
            if let Some(ref reactivity_data) = task.reactivity_data {
                plot_reactivity_stop_signals(gene_id, reactivity_data, &task.output_dir)?;
                plot_reactivity_mutation_signals(gene_id, reactivity_data, &task.output_dir)?;
            }
            
            // 更新进度
            let completed_count = completed.fetch_add(1, Ordering::Relaxed) + 1;
            let percentage = (completed_count as f64 * 100.0) / total_tasks as f64;
            print!("\r[Progressing] Plotting {}/{} ({:.1}%)", completed_count, total_tasks, percentage);
            std::io::stdout().flush()?;
            
            Ok(())
        })
        .collect();
    
    // 检查结果
    let mut success_count = 0;
    let mut error_count = 0;
    for result in results {
        match result {
            Ok(_) => success_count += 1,
            Err(e) => {
                error_count += 1;
                eprintln!("绘图错误: {}", e);
            }
        }
    }
    
    let elapsed = start_time.elapsed();
    
    // 覆盖进度显示，显示输出信息
    println!("\r[Outputs]  Success: {}. Failed: {}.", success_count, error_count);
    println!("    summary plot: {}/overall_scatter.png", output_dir);
    println!("    each gene plots: {}/", output_dir);
    println!("{}", crate::progress::format_time_used(elapsed));
    
    // 绘制总体scatter
    let scatter_path = format!("{}/overall_scatter.png", output_dir);
    plot_overall_scatter(&stats, coverage_threshold, depth_threshold, &scatter_path)?;
    
    // 记录完成日志
    logger.log(&format!("Signal distribution plotting completed"))?;
    logger.log(&format!("Successfully plotted genes: {}", success_count))?;
    if error_count > 0 {
        logger.log(&format!("Failed to plot genes: {}", error_count))?;
    }
    logger.log(&format!("Output Directory: {}", output_dir))?;
    logger.log(&format!("Total time: {:.2}s", elapsed.as_secs_f64()))?;
    
    if error_count > 0 {
        return Err(format!("有 {} 个基因区域绘图失败", error_count).into());
    }
    
    Ok(())
}

fn plot_stop_signals(
    gene_id: &str,
    mod_signals: &[SignalData],
    unmod_signals: &[SignalData],
    output_dir: &str,
) -> Result<(), Box<dyn Error + Send + Sync>> {
    let filename = format!("{}/{}_stop_signals.png", output_dir, gene_id.replace('/', "_"));
    let root = BitMapBackend::new(&filename, (1200, 800)).into_drawing_area();
    root.fill(&WHITE)?;
    
    let root = root.margin(10, 10, 10, 10);
    
    // 确定数据范围
    let all_positions: Vec<u32> = mod_signals.iter()
        .chain(unmod_signals.iter())
        .map(|d| d.position)
        .collect();
    
    let min_pos = all_positions.iter().min().unwrap_or(&0);
    let max_pos = all_positions.iter().max().unwrap_or(&0);
    
    let max_stop_signal = mod_signals.iter()
        .chain(unmod_signals.iter())
        .map(|d| d.stop_signal)
        .fold(0.0, f64::max);
    
    let mut chart = ChartBuilder::on(&root)
        .caption(format!("{} - Stop Signal Distribution", gene_id), ("sans-serif", 30))
        .x_label_area_size(40)
        .y_label_area_size(60)
        .build_cartesian_2d(*min_pos..*max_pos, 0.0..max_stop_signal * 1.1)?;
    
    chart
        .configure_mesh()
        .x_desc("Position")
        .y_desc("Stop Signal (PIPC/Depth)")
        .draw()?;
    
    // 绘制修饰样本数据（绿色）
    if !mod_signals.is_empty() {
        let mod_points: Vec<(u32, f64)> = mod_signals.iter()
            .map(|d| (d.position, d.stop_signal))
            .collect();
        
        chart.draw_series(LineSeries::new(
            mod_points.iter().map(|&(x, y)| (x, y)),
            GREEN.mix(0.8).stroke_width(2),
        ))?.label("Modified").legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], GREEN));
        
        chart.draw_series(mod_points.iter().map(|&(x, y)| {
            Circle::new((x, y), 2, GREEN.mix(0.8).filled())
        }))?;
    }
    
    // 绘制未修饰样本数据（蓝色）
    if !unmod_signals.is_empty() {
        let unmod_points: Vec<(u32, f64)> = unmod_signals.iter()
            .map(|d| (d.position, d.stop_signal))
            .collect();
        
        chart.draw_series(LineSeries::new(
            unmod_points.iter().map(|&(x, y)| (x, y)),
            BLUE.mix(0.8).stroke_width(2),
        ))?.label("Unmodified").legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], BLUE));
        
        chart.draw_series(unmod_points.iter().map(|&(x, y)| {
            Circle::new((x, y), 2, BLUE.mix(0.8).filled())
        }))?;
    }
    
    chart.configure_series_labels()
        .background_style(WHITE.mix(0.8))
        .border_style(BLACK)
        .draw()?;
    
    root.present()?;
    
    Ok(())
}

fn plot_mutation_signals(
    gene_id: &str,
    mod_signals: &[SignalData],
    unmod_signals: &[SignalData],
    output_dir: &str,
) -> Result<(), Box<dyn Error + Send + Sync>> {
    let filename = format!("{}/{}_mutation_signals.png", output_dir, gene_id.replace('/', "_"));
    let root = BitMapBackend::new(&filename, (1200, 800)).into_drawing_area();
    root.fill(&WHITE)?;
    
    let root = root.margin(10, 10, 10, 10);
    
    // 确定数据范围
    let all_positions: Vec<u32> = mod_signals.iter()
        .chain(unmod_signals.iter())
        .map(|d| d.position)
        .collect();
    
    let min_pos = all_positions.iter().min().unwrap_or(&0);
    let max_pos = all_positions.iter().max().unwrap_or(&0);
    
    let max_mutation_signal = mod_signals.iter()
        .chain(unmod_signals.iter())
        .map(|d| d.mutation_signal)
        .fold(0.0, f64::max);
    
    let mut chart = ChartBuilder::on(&root)
        .caption(format!("{} - Mutation Signal Distribution", gene_id), ("sans-serif", 30))
        .x_label_area_size(40)
        .y_label_area_size(60)
        .build_cartesian_2d(*min_pos..*max_pos, 0.0..max_mutation_signal * 1.1)?;
    
    chart
        .configure_mesh()
        .x_desc("Position")
        .y_desc("Mutation Signal (RFC/Depth)")
        .draw()?;
    
    // 绘制修饰样本数据（绿色）
    if !mod_signals.is_empty() {
        let mod_points: Vec<(u32, f64)> = mod_signals.iter()
            .map(|d| (d.position, d.mutation_signal))
            .collect();
        
        chart.draw_series(LineSeries::new(
            mod_points.iter().map(|&(x, y)| (x, y)),
            GREEN.mix(0.8).stroke_width(2),
        ))?.label("Modified").legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], GREEN));
        
        chart.draw_series(mod_points.iter().map(|&(x, y)| {
            Circle::new((x, y), 2, GREEN.mix(0.8).filled())
        }))?;
    }
    
    // 绘制未修饰样本数据（蓝色）
    if !unmod_signals.is_empty() {
        let unmod_points: Vec<(u32, f64)> = unmod_signals.iter()
            .map(|d| (d.position, d.mutation_signal))
            .collect();
        
        chart.draw_series(LineSeries::new(
            unmod_points.iter().map(|&(x, y)| (x, y)),
            BLUE.mix(0.8).stroke_width(2),
        ))?.label("Unmodified").legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], BLUE));
        
        chart.draw_series(unmod_points.iter().map(|&(x, y)| {
            Circle::new((x, y), 2, BLUE.mix(0.8).filled())
        }))?;
    }
    
    chart.configure_series_labels()
        .background_style(WHITE.mix(0.8))
        .border_style(BLACK)
        .draw()?;
    
    root.present()?;
    
    Ok(())
}

fn plot_reactivity_stop_signals(
    gene_id: &str,
    reactivity_data: &[ReactivityData],
    output_dir: &str,
) -> Result<(), Box<dyn Error + Send + Sync>> {
    let filename = format!("{}/{}_reactivity_stop_signals.png", output_dir, gene_id.replace('/', "_"));
    let root = BitMapBackend::new(&filename, (1200, 800)).into_drawing_area();
    root.fill(&WHITE)?;
    
    let root = root.margin(10, 10, 10, 10);
    
    // 确定数据范围
    let all_positions: Vec<u32> = reactivity_data.iter().map(|d| d.position).collect();
    let min_pos = all_positions.iter().min().unwrap_or(&0);
    let max_pos = all_positions.iter().max().unwrap_or(&0);
    
    let min_pos_f = *min_pos as f64;
    let max_pos_f = *max_pos as f64;
    
    let max_reactivity = reactivity_data.iter()
        .map(|d| d.reactivity_stop)
        .fold(0.0, f64::max);
    
    let mut chart = ChartBuilder::on(&root)
        .caption(format!("{} - Stop Reactivity Signal Distribution", gene_id), ("sans-serif", 30))
        .x_label_area_size(40)
        .y_label_area_size(60)
        .build_cartesian_2d(min_pos_f..max_pos_f, 0.0..max_reactivity * 1.1)?;
    
    chart
        .configure_mesh()
        .x_desc("Position")
        .y_desc("Stop Reactivity Signal")
        .draw()?;
    
    // 绘制stop reactivity柱形图
    for data in reactivity_data {
        let x = data.position as f64;
        let y_stop = data.reactivity_stop;
        
        // 根据reactivity值选择颜色
        let stop_color = if y_stop <= 0.2 {
            RGBColor(128, 128, 128) // 灰色
        } else if y_stop <= 0.4 {
            RGBColor(255, 255, 0) // 黄色
        } else {
            RGBColor(255, 0, 0) // 红色
        };
        
        // 绘制stop reactivity柱形图
        if y_stop > 0.0 {
            chart.draw_series(std::iter::once(Rectangle::new(
                [(x - 0.2, 0.0), (x + 0.2, y_stop)],
                stop_color.filled(),
            )))?;
        }
    }
    
    // 添加图例
    chart.draw_series(std::iter::once(Rectangle::new(
        [(0.0, 0.0), (10.0, 0.0)],
        RGBColor(128, 128, 128).filled(),
    )))?.label("0-0.2 (Gray)");
    
    chart.draw_series(std::iter::once(Rectangle::new(
        [(0.0, 0.0), (10.0, 0.0)],
        RGBColor(255, 255, 0).filled(),
    )))?.label("0.2-0.4 (Yellow)");
    
    chart.draw_series(std::iter::once(Rectangle::new(
        [(0.0, 0.0), (10.0, 0.0)],
        RGBColor(255, 0, 0).filled(),
    )))?.label("0.4-1.0 (Red)");
    
    chart.configure_series_labels()
        .background_style(WHITE.mix(0.8))
        .border_style(BLACK)
        .draw()?;
    
    root.present()?;
    
    Ok(())
}

fn plot_reactivity_mutation_signals(
    gene_id: &str,
    reactivity_data: &[ReactivityData],
    output_dir: &str,
) -> Result<(), Box<dyn Error + Send + Sync>> {
    let filename = format!("{}/{}_reactivity_mutation_signals.png", output_dir, gene_id.replace('/', "_"));
    let root = BitMapBackend::new(&filename, (1200, 800)).into_drawing_area();
    root.fill(&WHITE)?;
    
    let root = root.margin(10, 10, 10, 10);
    
    // 确定数据范围
    let all_positions: Vec<u32> = reactivity_data.iter().map(|d| d.position).collect();
    let min_pos = all_positions.iter().min().unwrap_or(&0);
    let max_pos = all_positions.iter().max().unwrap_or(&0);
    
    let min_pos_f = *min_pos as f64;
    let max_pos_f = *max_pos as f64;
    
    let max_reactivity = reactivity_data.iter()
        .map(|d| d.reactivity_mutation)
        .fold(0.0, f64::max);
    
    let mut chart = ChartBuilder::on(&root)
        .caption(format!("{} - Mutation Reactivity Signal Distribution", gene_id), ("sans-serif", 30))
        .x_label_area_size(40)
        .y_label_area_size(60)
        .build_cartesian_2d(min_pos_f..max_pos_f, 0.0..max_reactivity * 1.1)?;
    
    chart
        .configure_mesh()
        .x_desc("Position")
        .y_desc("Mutation Reactivity Signal")
        .draw()?;
    
    // 绘制mutation reactivity柱形图
    for data in reactivity_data {
        let x = data.position as f64;
        let y_mut = data.reactivity_mutation;
        
        // 根据reactivity值选择颜色
        let mut_color = if y_mut <= 0.2 {
            RGBColor(128, 128, 128) // 灰色
        } else if y_mut <= 0.4 {
            RGBColor(255, 255, 0) // 黄色
        } else {
            RGBColor(255, 0, 0) // 红色
        };
        
        // 绘制mutation reactivity柱形图
        if y_mut > 0.0 {
            chart.draw_series(std::iter::once(Rectangle::new(
                [(x - 0.2, 0.0), (x + 0.2, y_mut)],
                mut_color.filled(),
            )))?;
        }
    }
    
    // 添加图例
    chart.draw_series(std::iter::once(Rectangle::new(
        [(0.0, 0.0), (10.0, 0.0)],
        RGBColor(128, 128, 128).filled(),
    )))?.label("0-0.2 (Gray)");
    
    chart.draw_series(std::iter::once(Rectangle::new(
        [(0.0, 0.0), (10.0, 0.0)],
        RGBColor(255, 255, 0).filled(),
    )))?.label("0.2-0.4 (Yellow)");
    
    chart.draw_series(std::iter::once(Rectangle::new(
        [(0.0, 0.0), (10.0, 0.0)],
        RGBColor(255, 0, 0).filled(),
    )))?.label("0.4-1.0 (Red)");
    
    chart.configure_series_labels()
        .background_style(WHITE.mix(0.8))
        .border_style(BLACK)
        .draw()?;
    
    root.present()?;
    
    Ok(())
}

fn plot_overall_scatter(
    stats: &HashMap<String, GeneStat>,
    coverage_threshold: f64,
    depth_threshold: f64,
    output_path: &str,
) -> Result<(), Box<dyn Error>> {
    let root = BitMapBackend::new(output_path, (1200, 900)).into_drawing_area();
    root.fill(&WHITE)?;
    
    let mut all_coverage: Vec<f64> = Vec::new();
    let mut all_depth: Vec<f64> = Vec::new();
    for s in stats.values() {
        all_coverage.push(s.coverage);
        all_depth.push(s.avg_depth);
    }
    
    // coverage不使用log10转换，depth仍然使用log10转换
    let log_depth: Vec<f64> = all_depth.iter().map(|&x| (x + 1.0).log10()).collect();
    
    // 强制x轴和y轴从0开始
    let min_cov = 0.0;
    let max_cov = *all_coverage.iter().max_by(|a, b| a.partial_cmp(b).unwrap()).unwrap_or(&1.0);
    let min_depth = 0.0;
    let max_depth = *log_depth.iter().max_by(|a, b| a.partial_cmp(b).unwrap()).unwrap_or(&2.0);
    
    let mut chart = ChartBuilder::on(&root)
        .caption("Gene Coverage-Depth Scatter", ("sans-serif", 36))
        .x_label_area_size(60)
        .y_label_area_size(80)
        .margin(20)
        .build_cartesian_2d(min_cov..max_cov, min_depth..max_depth)?;
    
    chart.configure_mesh()
        .x_desc("Coverage")
        .y_desc("log10(Average Depth)")
        .draw()?;
    
    // 绘制阈值线 (coverage不使用log10，depth使用log10)
    let log_depth_threshold = (depth_threshold + 1.0).log10();
    
    // 垂直coverage阈值线
    chart.draw_series(LineSeries::new(
        vec![(coverage_threshold, min_depth), (coverage_threshold, max_depth)],
        RGBColor(255, 0, 0).mix(0.7).stroke_width(2),
    ))?;
    
    // 水平depth阈值线
    chart.draw_series(LineSeries::new(
        vec![(min_cov, log_depth_threshold), (max_cov, log_depth_threshold)],
        RGBColor(0, 0, 255).mix(0.7).stroke_width(2),
    ))?;
    
    // 统计高低组数量
    let mut high_count = 0;
    let mut low_count = 0;
    for s in stats.values() {
        if s.coverage >= coverage_threshold && s.avg_depth >= depth_threshold {
            high_count += 1;
        } else {
            low_count += 1;
        }
    }
    
    // 画点
    for s in stats.values() {
        let cov = s.coverage;
        let log_depth = (s.avg_depth + 1.0).log10();
        let color = if s.coverage >= coverage_threshold && s.avg_depth >= depth_threshold { 
            RED 
        } else { 
            RGBColor(180, 180, 180) 
        };
        chart.draw_series(std::iter::once(Circle::new((cov, log_depth), 7, color.filled())))?;
    }
    
    // 添加图例和数量标注
    chart.draw_series(std::iter::once(Circle::new((0.0, 0.0), 0, RED.filled())))?
        .label(&format!("High (N={})", high_count));
    chart.draw_series(std::iter::once(Circle::new((0.0, 0.0), 0, RGBColor(180, 180, 180).filled())))?
        .label(&format!("Low (N={})", low_count));
    
    chart.configure_series_labels()
        .background_style(WHITE.mix(0.8))
        .border_style(BLACK)
        .draw()?;
    
    // 添加阈值标注
    chart.draw_series(std::iter::once(Text::new(
        format!("Coverage ≥ {:.3}", coverage_threshold),
        (max_cov - 0.1, max_depth - 0.3),
        ("sans-serif", 20).into_font().color(&RGBColor(255, 0, 0)),
    )))?;
    
    chart.draw_series(std::iter::once(Text::new(
        format!("Depth ≥ {:.3}", depth_threshold),
        (min_cov + 0.1, max_depth - 0.3),
        ("sans-serif", 20).into_font().color(&RGBColor(0, 0, 255)),
    )))?;
    
    root.present()?;
    Ok(())
}

// ==================== SVG Plotting Functions ====================

#[derive(Debug, Clone)]
pub struct SvgReactivityData {
    pub position: u32,
    pub base: char,
    pub reactivity: f64,
}

#[derive(Debug, Clone)]
pub struct ColorRange {
    pub min: f64,
    pub max: f64,
    pub color: String,
}

#[derive(Debug, Clone)]
pub struct BaseColorRanges {
    pub a: Vec<ColorRange>,
    pub t: Vec<ColorRange>,
    pub c: Vec<ColorRange>,
    pub g: Vec<ColorRange>,
}

impl Default for BaseColorRanges {
    fn default() -> Self {
        Self {
            a: vec![
                ColorRange { min: 0.0, max: 0.3, color: "rgba(245, 245, 245, 0.8)".to_string() },
                ColorRange { min: 0.3, max: 0.65, color: "rgba(255, 165, 0, 0.8)".to_string() },
                ColorRange { min: 0.65, max: 1.0, color: "rgba(255, 0, 0, 0.8)".to_string() },
            ],
            t: vec![
                ColorRange { min: 0.0, max: 0.2, color: "rgba(245, 245, 245, 0.8)".to_string() },
                ColorRange { min: 0.2, max: 0.55, color: "rgba(255, 165, 0, 0.8)".to_string() },
                ColorRange { min: 0.55, max: 1.0, color: "rgba(255, 0, 0, 0.8)".to_string() },
            ],
            c: vec![
                ColorRange { min: 0.0, max: 0.3, color: "rgba(245, 245, 245, 0.8)".to_string() },
                ColorRange { min: 0.3, max: 0.6, color: "rgba(255, 165, 0, 0.8)".to_string() },
                ColorRange { min: 0.6, max: 1.0, color: "rgba(255, 0, 0, 0.8)".to_string() },
            ],
            g: vec![
                ColorRange { min: 0.0, max: 0.2, color: "rgba(245, 245, 245, 0.8)".to_string() },
                ColorRange { min: 0.2, max: 0.5, color: "rgba(255, 165, 0, 0.8)".to_string() },
                ColorRange { min: 0.5, max: 1.0, color: "rgba(255, 0, 0, 0.8)".to_string() },
            ],
        }
    }
}

/// Load reactivity data from CSV file for SVG plotting
pub fn load_svg_reactivity_data(
    csv_file: &str, 
    signal_type: &str, 
    strand_filter: &str
) -> Result<Vec<SvgReactivityData>, Box<dyn Error>> {
    let mut data = Vec::new();
    let file = File::open(csv_file)?;
    let reader = BufReader::new(file);
    let lines: Vec<String> = reader.lines().collect::<Result<Vec<_>, _>>()?;
    
    if lines.is_empty() {
        return Ok(data);
    }
    
    // Auto-detect CSV format from first line
    let first_line = &lines[0];
    let fields: Vec<&str> = first_line.split(',').collect();
    
    let has_header = if fields.len() >= 3 {
        fields[2].parse::<u32>().is_err() // If position field is not numeric, it's a header
    } else {
        false
    };
    
    let has_strand_column = fields.len() >= 5;
    
    // Auto-detect signal columns based on header
    let mut signal_columns: Vec<(String, usize)> = Vec::new();
    if has_header {
        // Parse header to find signal columns
        for (i, field) in fields.iter().enumerate() {
            if field.contains("reactivity_") {
                let signal_name = field.replace("reactivity_", "");
                signal_columns.push((signal_name, i));
            } else if i >= 4 && *field != "Strand" && *field != "Position" && *field != "Base" {
                // Assume this is a signal column if it's not a standard column
                signal_columns.push((field.to_string(), i));
            }
        }
    } else {
        // For files without header, assume last column is the signal
        if fields.len() >= 5 {
            signal_columns.push(("score".to_string(), fields.len() - 1));
        }
    }
    
    // If no signal columns detected, use default
    if signal_columns.is_empty() {
        signal_columns.push(("score".to_string(), if fields.len() == 6 { 4 } else if fields.len() == 5 { 4 } else { 3 }));
    }
    
    // Process data lines
    let start_line = if has_header { 1 } else { 0 };
    
    for line in &lines[start_line..] {
        let fields: Vec<&str> = line.split(',').collect();
        
        // Skip empty lines
        if fields.is_empty() || fields[0].trim().is_empty() {
            continue;
        }
        
        // Check strand filter
        let should_process = if has_strand_column && fields.len() >= 5 {
            match strand_filter {
                "+" => fields[1].trim() == "+",
                "-" => fields[1].trim() == "-",
                "both" => true,
                _ => fields[1].trim() == "+", // Default to positive strand
            }
        } else {
            true // No strand column, process all data
        };
        
        if should_process {
            // Process each signal type
            for (signal_name, column_index) in &signal_columns {
                // Check if this signal type should be processed
                let should_process_signal = match signal_type {
                    "all" => true,
                    signal => signal_name == signal,
                };
                
                if should_process_signal && column_index < &fields.len() {
                    if let Ok(position) = fields[2].parse::<u32>() {
                        if let Ok(reactivity) = fields[*column_index].parse::<f64>() {
                            let base = normalize_base(fields[3].chars().next().unwrap_or('N'));
                            data.push(SvgReactivityData {
                                position,
                                base,
                                reactivity,
                            });
                        }
                    }
                }
            }
        }
    }
    
    Ok(data)
}

/// Normalize base character (U -> T, uppercase)
fn normalize_base(base: char) -> char {
    let base = base.to_uppercase().next().unwrap_or('N');
    if base == 'U' { 'T' } else { base }
}

/// Load reference sequence from FASTA or plain text file
pub fn load_reference_sequence(ref_path: &str) -> Result<HashMap<u32, char>, Box<dyn Error>> {
    let mut seq = String::new();
    let lower_path = ref_path.to_lowercase();
    
    if lower_path.ends_with(".fa") || lower_path.ends_with(".fasta") || 
       lower_path.ends_with(".fa.gz") || lower_path.ends_with(".fasta.gz") {
        // Parse FASTA file
        let file = File::open(ref_path)?;
        let reader = BufReader::new(file);
        let mut in_sequence = false;
        
        for line in reader.lines() {
            let line = line?;
            if line.starts_with('>') {
                in_sequence = true;
                continue;
            }
            if in_sequence && !line.trim().is_empty() {
                seq.push_str(&line.trim());
            }
        }
    } else {
        // Treat as plain text sequence
        let file = File::open(ref_path)?;
        let reader = BufReader::new(file);
        
        for line in reader.lines() {
            let line = line?;
            if !line.trim().is_empty() && !line.starts_with('>') {
                seq.push_str(&line.trim());
            }
        }
    }
    
    if seq.is_empty() {
        return Err("Reference sequence is empty or not found".into());
    }
    
    // Convert to position->base map (1-based positions)
    let mut pos_to_base = HashMap::new();
    for (i, base) in seq.chars().enumerate() {
        pos_to_base.insert((i + 1) as u32, normalize_base(base));
    }
    
    Ok(pos_to_base)
}

/// Find best shift for a given base type by maximizing matches with reference sequence
pub fn find_best_shift_for_base(
    base_data: &[&SvgReactivityData],
    ref_pos2base: &HashMap<u32, char>,
    max_shift: u32,
) -> i32 {
    if base_data.is_empty() {
        return 0;
    }
    
    // Build observed map
    let obs: HashMap<u32, char> = base_data
        .iter()
        .map(|data| (data.position, data.base))
        .collect();
    
    let mut best_shift = 0;
    let mut best_match = -1;
    
    for s in -(max_shift as i32)..=(max_shift as i32) {
        let mut match_count = 0;
        
        // Compare positions that exist in both after shifting
        for (&pos, &base) in &obs {
            let ref_pos = (pos as i32 + s) as u32;
            if let Some(&ref_base) = ref_pos2base.get(&ref_pos) {
                if base == ref_base {
                    match_count += 1;
                }
            }
        }
        
        if match_count > best_match {
            best_match = match_count;
            best_shift = s;
        }
    }
    
    best_shift
}

/// Determine color based on reactivity value and color ranges
fn determine_color(shape_value: f64, color_ranges: &[ColorRange]) -> String {
    if shape_value == -1.0 {
        return "rgba(169, 169, 169, 0.8)".to_string();
    }
    
    for range in color_ranges {
        if shape_value >= range.min && shape_value <= range.max {
            return range.color.clone();
        }
    }
    
    "rgba(169, 169, 169, 0.8)".to_string()
}

/// Draw color bar legend for a base type
fn draw_color_bar(
    writer: &mut BufWriter<File>,
    color_ranges: &[ColorRange],
    base: char,
    position: usize,
) -> Result<(), Box<dyn Error>> {
    let bar_width = 10.0;
    let bar_height = 100.0;
    let bar_x = 20.0;
    let bar_y = 20.0 + (position as f64) * (bar_height + 50.0);
    let segment_height = bar_height / (color_ranges.len() + 1) as f64; // +1 for -1 value
    
    writeln!(writer, r#"<g id="color-bar-{}">"#, base)?;
    
    // Draw -1 value bar
    writeln!(writer, r#"<rect x="{}" y="{}" width="{}" height="{}" fill="rgba(169, 169, 169, 0.8)" stroke="black" stroke-width="1" />"#,
        bar_x, bar_y, bar_width, segment_height)?;
    writeln!(writer, r#"<text x="{}" y="{}" font-family="Arial" font-size="6" fill="black" transform="rotate(90, {}, {})">-1</text>"#,
        bar_x + bar_width + 5.0, bar_y + segment_height / 2.0, 
        bar_x + bar_width + 5.0, bar_y + segment_height / 2.0)?;
    
    // Draw base label
    writeln!(writer, r#"<text x="{}" y="{}" font-family="Arial" font-size="14" fill="black">{}</text>"#,
        bar_x - 5.0, bar_y - 10.0, base)?;
    
    // Draw color range bars
    for (i, range) in color_ranges.iter().enumerate() {
        let segment_y = bar_y + segment_height + (i as f64) * segment_height;
        writeln!(writer, r#"<rect x="{}" y="{}" width="{}" height="{}" fill="{}" stroke="black" stroke-width="1" />"#,
            bar_x, segment_y, bar_width, segment_height, range.color)?;
        writeln!(writer, r#"<text x="{}" y="{}" font-family="Arial" font-size="6" fill="black" transform="rotate(90, {}, {})">{:.2} - {:.2}</text>"#,
            bar_x + bar_width + 5.0, segment_y + segment_height / 2.0,
            bar_x + bar_width + 5.0, segment_y + segment_height / 2.0,
            range.min, range.max)?;
    }
    
    writeln!(writer, "</g>")?;
    Ok(())
}

/// Plot reactivity data to SVG without alignment
pub fn plot_reactivity_to_svg(
    reactivity_file: &str,
    svg_template: &str,
    output_file: &str,
    bases_filter: &str,
    signal_type: &str,
    strand_filter: &str,
    color_ranges: Option<BaseColorRanges>,
) -> Result<(), Box<dyn Error>> {
    let reactivity_data = load_svg_reactivity_data(reactivity_file, signal_type, strand_filter)?;
    let svg_content = std::fs::read_to_string(svg_template)?;
    let color_ranges = color_ranges.unwrap_or_default();
    
    // Filter data by bases and group by base type
    let bases: Vec<char> = bases_filter.chars().collect();
    let mut base_groups: HashMap<char, Vec<&SvgReactivityData>> = HashMap::new();
    
    for data in &reactivity_data {
        if bases.contains(&data.base) {
            base_groups.entry(data.base).or_insert_with(Vec::new).push(data);
        }
    }
    
    // Create output file
    let output = File::create(output_file)?;
    let mut writer = BufWriter::new(output);
    
    // Regex patterns for parsing SVG
    let position_regex = Regex::new(r"title>(\d+)")?;
    let x_regex = Regex::new(r#"x="([^"]+)""#)?;
    let y_regex = Regex::new(r#"y="([^"]+)""#)?;
    
    // First, write all original SVG content except the closing </svg> tag (like Python script)
    let svg_lines: Vec<&str> = svg_content.lines().collect();
    for line in &svg_lines {
        if line.contains("</svg>") {
            break; // Stop before the closing tag
        }
        writeln!(writer, "{}", line)?;
    }
    
    // Process each base type in the order specified by bases_filter (like Python script)
    for (base_idx, base) in bases_filter.chars().enumerate() {
        if let Some(data_list) = base_groups.get(&base) {
            let color_ranges_for_base = match base {
                'A' => &color_ranges.a,
                'T' => &color_ranges.t,
                'C' => &color_ranges.c,
                'G' => &color_ranges.g,
                _ => continue,
            };
        
            // Create position to reactivity mapping
            let pos_to_reactivity: HashMap<u32, f64> = data_list
                .iter()
                .map(|data| (data.position, data.reactivity))
                .collect();
            
            // Process entire SVG content for this base (like Python add_shape_to_svg function)
            for template_line in svg_content.lines() {
                if template_line.contains("title>") && !template_line.contains("numbering-label") {
                    if let Some(caps) = position_regex.captures(template_line) {
                        if let Ok(svg_position) = caps[1].parse::<u32>() {
                            // CSV position and SVG position are both 0-based, direct mapping
                            let csv_position = svg_position;
                            if csv_position > 0 {
                                if let Some(&reactivity) = pos_to_reactivity.get(&csv_position) {
                                    let color = determine_color(reactivity, color_ranges_for_base);
                                    
                                    // Extract x and y coordinates
                                    if let (Some(x_caps), Some(y_caps)) = (x_regex.captures(template_line), y_regex.captures(template_line)) {
                                        let x = &x_caps[1];
                                        let y = &y_caps[1];
                                        
                                        // Add circle element (wrapped in <g> like Python)
                                        writeln!(writer, r#"<g><circle cx="{}" cy="{}" r="2.8" stroke="none" stroke-width="0.2" fill="{}"/></g>"#,
                                            x, y, color)?;
                                    }
                                }
                            }
                        }
                    }
                } else {
                    // Write all other lines (like Python script)
                    writeln!(writer, "{}", template_line)?;
                }
            }
            
            // Draw color bar for this base
            draw_color_bar(&mut writer, color_ranges_for_base, base, base_idx)?;
        }
    }
    
    // Write final closing SVG tag (like Python script)
    writeln!(writer, "</svg>")?;
    
    writer.flush()?;
    Ok(())
}

/// Plot multiple signal types to SVG (auto-detect signals from CSV)
pub fn plot_multiple_signals_to_svg(
    reactivity_file: &str,
    svg_template: &str,
    output_dir: &str,
    bases_filter: &str,
    strand_filter: &str,
    ref_sequence_file: Option<&str>,
    max_shift: u32,
    color_ranges: Option<BaseColorRanges>,
) -> Result<(), Box<dyn Error>> {
    // Auto-detect signal types from CSV file
    let signal_types = detect_signal_types(reactivity_file)?;
    
    // Create output directory
    std::fs::create_dir_all(output_dir)?;
    
    // Plot each signal type
    for signal_type in &signal_types {
        let output_file = format!("{}/rna_structure_colored_{}.svg", output_dir, signal_type);
        
        if let Some(ref_file) = ref_sequence_file {
            // Use reference sequence alignment
            plot_reactivity_to_svg_with_alignment(
                reactivity_file,
                svg_template,
                &output_file,
                bases_filter,
                signal_type,
                strand_filter,
                ref_file,
                max_shift,
                color_ranges.clone(),
            )?;
        } else {
            // No alignment
            plot_reactivity_to_svg(
                reactivity_file,
                svg_template,
                &output_file,
                bases_filter,
                signal_type,
                strand_filter,
                color_ranges.clone(),
            )?;
        }
    }
    
    Ok(())
}

/// Auto-detect signal types from CSV file header
fn detect_signal_types(csv_file: &str) -> Result<Vec<String>, Box<dyn Error>> {
    let file = File::open(csv_file)?;
    let reader = BufReader::new(file);
    let first_line = reader.lines().next().ok_or("Empty CSV file")??;
    let fields: Vec<&str> = first_line.split(',').collect();
    
    let mut signal_types = Vec::new();
    for field in &fields {
        if field.contains("reactivity_") {
            let signal_name = field.replace("reactivity_", "");
            signal_types.push(signal_name);
        }
    }
    
    // If no reactivity columns found, assume single signal
    if signal_types.is_empty() {
        signal_types.push("score".to_string());
    }
    
    Ok(signal_types)
}

/// Plot reactivity data to SVG with alignment
pub fn plot_reactivity_to_svg_with_alignment(
    reactivity_file: &str,
    svg_template: &str,
    output_file: &str,
    bases_filter: &str,
    signal_type: &str,
    strand_filter: &str,
    ref_sequence_file: &str,
    max_shift: u32,
    color_ranges: Option<BaseColorRanges>,
) -> Result<(), Box<dyn Error>> {
    let reactivity_data = load_svg_reactivity_data(reactivity_file, signal_type, strand_filter)?;
    let svg_content = std::fs::read_to_string(svg_template)?;
    let color_ranges = color_ranges.unwrap_or_default();
    
    // Load reference sequence
    let ref_pos2base = load_reference_sequence(ref_sequence_file)?;
    
    // Filter data by bases and group by base type
    let bases: Vec<char> = bases_filter.chars().collect();
    let mut base_groups: HashMap<char, Vec<&SvgReactivityData>> = HashMap::new();
    
    for data in &reactivity_data {
        if bases.contains(&data.base) {
            base_groups.entry(data.base).or_insert_with(Vec::new).push(data);
        }
    }
    
    // Calculate best shift for each base type
    let mut base_to_shift: HashMap<char, i32> = HashMap::new();
    for base in &bases {
        let base_data: Vec<&SvgReactivityData> = reactivity_data
            .iter()
            .filter(|data| data.base == *base)
            .collect();
        let shift = find_best_shift_for_base(&base_data, &ref_pos2base, max_shift);
        base_to_shift.insert(*base, shift);
    }
    
    // Create output file
    let output = File::create(output_file)?;
    let mut writer = BufWriter::new(output);
    
    // Regex patterns for parsing SVG
    let position_regex = Regex::new(r"title>(\d+)")?;
    let x_regex = Regex::new(r#"x="([^"]+)""#)?;
    let y_regex = Regex::new(r#"y="([^"]+)""#)?;
    
    // First, write all original SVG content except the closing </svg> tag (like Python script)
    let svg_lines: Vec<&str> = svg_content.lines().collect();
    for line in &svg_lines {
        if line.contains("</svg>") {
            break; // Stop before the closing tag
        }
        writeln!(writer, "{}", line)?;
    }
    
    // Process each base type (like Python script) - each base gets its own SVG
    for (base_idx, (base, data_list)) in base_groups.iter().enumerate() {
        let color_ranges_for_base = match base {
            'A' => &color_ranges.a,
            'T' => &color_ranges.t,
            'C' => &color_ranges.c,
            'G' => &color_ranges.g,
            _ => continue,
        };
        
        // Get shift for this base
        let shift = *base_to_shift.get(base).unwrap_or(&0);
        
        // Create position to reactivity mapping
        let pos_to_reactivity: HashMap<u32, f64> = data_list
            .iter()
            .map(|data| (data.position, data.reactivity))
            .collect();
        
        // Process entire SVG content for this base (like Python add_shape_to_svg function)
        for template_line in svg_content.lines() {
            if template_line.contains("title>") && !template_line.contains("numbering-label") {
                if let Some(caps) = position_regex.captures(template_line) {
                    if let Ok(svg_position) = caps[1].parse::<u32>() {
                        // Apply shift: pos_data = pos_svg - shift (like Python version)
                        let csv_position = (svg_position as i32 - shift) as u32;
                        if csv_position > 0 {
                            if let Some(&reactivity) = pos_to_reactivity.get(&csv_position) {
                                let color = determine_color(reactivity, color_ranges_for_base);
                                
                                // Extract x and y coordinates
                                if let (Some(x_caps), Some(y_caps)) = (x_regex.captures(template_line), y_regex.captures(template_line)) {
                                    let x = &x_caps[1];
                                    let y = &y_caps[1];
                                    
                                    // Add circle element (wrapped in <g> like Python)
                                    writeln!(writer, r#"<g><circle cx="{}" cy="{}" r="2.8" stroke="none" stroke-width="0.2" fill="{}"/></g>"#,
                                        x, y, color)?;
                                }
                            }
                        }
                    }
                }
            } else {
                // Write all other lines (like Python script)
                writeln!(writer, "{}", template_line)?;
            }
        }
        
        // Draw color bar for this base
        draw_color_bar(&mut writer, color_ranges_for_base, *base, base_idx)?;
    }
    
    // Write final closing SVG tag (like Python script)
    writeln!(writer, "</svg>")?;
    
    writer.flush()?;
    Ok(())
} 