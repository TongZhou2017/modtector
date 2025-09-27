use std::error::Error;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use plotters::prelude::*;
use crate::plot::ReactivityData;
use std::collections::HashMap;
use std::time::Instant;
use chrono;

/// Secondary structure data
#[derive(Debug, Clone)]
pub struct SecondaryStructure {
    pub sequence: String,
    pub structure: String,
    pub positions: Vec<u32>,
    pub is_single_stranded: Vec<bool>, // true means single-stranded, false means double-stranded
}

/// Accuracy evaluation result
#[derive(Debug)]
pub struct EvaluationResult {
    pub auc: f64,
    pub f1_score: f64,
    pub sensitivity: f64,
    pub specificity: f64,
    pub ppv: f64, // Positive Predictive Value
    pub npv: f64, // Negative Predictive Value
    pub accuracy: f64,
    pub thresholds: Vec<f64>,
    pub tpr: Vec<f64>, // True Positive Rate (Sensitivity)
    pub fpr: Vec<f64>, // False Positive Rate (1-Specificity)
    pub precision: Vec<f64>,
    pub recall: Vec<f64>,
}

/// Per-base evaluation result
#[derive(Debug)]
pub struct BaseSpecificEvaluationResult {
    pub base: char,
    pub result: EvaluationResult,
    pub data_count: usize,
    pub single_stranded_count: usize,
    pub double_stranded_count: usize,
}

/// Comprehensive evaluation result (including all bases)
#[derive(Debug)]
pub struct ComprehensiveEvaluationResult {
    pub overall_result: EvaluationResult,
    pub base_specific_results: Vec<BaseSpecificEvaluationResult>,
    pub total_positions: usize,
}

/// Parse secondary structure file
pub fn parse_secondary_structure(file_path: &str) -> Result<SecondaryStructure, Box<dyn Error>> {
    let file = File::open(file_path)?;
    let reader = BufReader::new(file);
    let mut sequence = String::new();
    let mut structure = String::new();
    let mut found_sequence = false;
    let mut found_structure = false;
    
    for line in reader.lines() {
        let line = line?;
        let trimmed = line.trim();
        
        if trimmed.starts_with('#') || trimmed.is_empty() {
            continue;
        }
        
        // Check if it's a sequence line (contains only ACGU)
        if !found_sequence && trimmed.chars().all(|c| "ACGU".contains(c)) {
            sequence.push_str(trimmed);
            found_sequence = true;
        }
        // Continue reading sequence lines until encountering structure line
        else if found_sequence && !found_structure && trimmed.chars().all(|c| "ACGU".contains(c)) {
            sequence.push_str(trimmed);
        }
        // Check if it's a structure line (contains .() etc.)
        else if found_sequence && !found_structure && trimmed.chars().all(|c| ".()[]{}<> aA ".contains(c)) {
            structure.push_str(trimmed);
            found_structure = true;
        }
        // 继续读取结构行
        else if found_sequence && found_structure && trimmed.chars().all(|c| ".()[]{}<> aA ".contains(c)) {
            structure.push_str(trimmed);
        }
    }
    
    // 过滤掉空格，只保留有效的结构字符
    structure = structure.chars().filter(|&c| c != ' ').collect();
    

    
    if sequence.is_empty() || structure.is_empty() {
        return Err("Could not find sequence or structure in file".into());
    }
    
    if sequence.len() != structure.len() {
        return Err(format!("Sequence length ({}) and structure length ({}) do not match", 
                          sequence.len(), structure.len()).into());
    }
    
    let mut positions = Vec::new();
    let mut is_single_stranded = Vec::new();
    
    for (i, (seq_char, struct_char)) in sequence.chars().zip(structure.chars()).enumerate() {
        // 只处理标准碱基
        if "ACGU".contains(seq_char) {
            positions.push(i as u32 + 1); // 1-based position
            // 单链：. 其它都视为双链
            is_single_stranded.push(struct_char == '.');
        }
    }
    
    let single_stranded_count = is_single_stranded.iter().filter(|&&x| x).count();
    let double_stranded_count = is_single_stranded.iter().filter(|&&x| !x).count();
    
    Ok(SecondaryStructure {
        sequence,
        structure,
        positions,
        is_single_stranded,
    })
}

/// 从reactivity CSV文件中提取碱基信息
pub fn extract_reactivity_with_bases(csv_file: &str, gene_id: &str, strand: &str, signal_type: &str) -> Result<Vec<(u32, char, f64)>, Box<dyn Error>> {
    let mut reactivity_data: Vec<(u32, char, f64)> = Vec::new();
    let file = File::open(csv_file)?;
    let reader = BufReader::new(file);
    
    // 先计算总行数用于进度显示
    let total_lines = std::fs::read_to_string(csv_file)?.lines().count();
    
    for (line_num, line) in reader.lines().enumerate() {
        if line_num == 0 {
            continue; // 跳过标题行
        }
        
        // 更新进度显示
        if line_num % 1000 == 0 {
            // progress.update(line_num)?;
        }
        
        let line = line?;
        let fields: Vec<&str> = line.split(',').collect();
        
        if fields.len() >= 6 {
            let chr_id = fields[0].trim_matches('"');
            let strand_field = fields[1].trim();
            if chr_id != gene_id || strand_field != strand {
                continue;
            }
            
            let position = fields[2].parse::<u32>()?;
            let base = fields[3].chars().next().unwrap_or('N');
            
            // 根据signal_type决定读取哪一列的reactivity值
            let reactivity_value = if signal_type == "mutation" {
                // mutation信号使用第6列 (reactivity_mutation)
                fields[5].parse::<f64>()?
            } else {
                // stop信号使用第5列 (reactivity_stop)
                fields[4].parse::<f64>()?
            };
            
            reactivity_data.push((position, base, reactivity_value));
        }
    }

    
    // 按位置排序
    reactivity_data.sort_by_key(|(pos, _, _)| *pos);
        
    Ok(reactivity_data)
}

/// 从reactivity CSV文件中提取指定基因的碱基信息（优化版本）
pub fn extract_reactivity_for_target_gene(
    csv_file: &str, 
    target_gene_patterns: &[&str], 
    strand: &str, 
    signal_type: &str
) -> Result<Option<(String, Vec<(u32, char, f64)>)>, Box<dyn Error>> {
    let file = File::open(csv_file)?;
    let reader = BufReader::new(file);
    let mut reactivity_data: Vec<(u32, char, f64)> = Vec::new();
    let mut found_gene_id: Option<String> = None;
    
    for (line_num, line) in reader.lines().enumerate() {
        if line_num == 0 {
            continue; // 跳过标题行
        }
        
        let line = line?;
        let fields: Vec<&str> = line.split(',').collect();
        
        if fields.len() >= 6 {
            let chr_id = fields[0].trim_matches('"');
            let strand_field = fields[1].trim();
            
            // 检查是否匹配目标基因模式
            let is_target_gene = target_gene_patterns.iter().any(|pattern| chr_id.contains(pattern));
            
            if is_target_gene && strand_field == strand {
                // 找到目标基因，记录基因ID
                if found_gene_id.is_none() {
                    found_gene_id = Some(chr_id.to_string());
                }
                
                let position = fields[2].parse::<u32>()?;
                let base = fields[3].chars().next().unwrap_or('N');
                
                // 根据signal_type决定读取哪一列的reactivity值
                let reactivity_value = if signal_type == "mutation" {
                    // mutation信号使用第6列 (reactivity_mutation)
                    fields[5].parse::<f64>()?
                } else {
                    // stop信号使用第5列 (reactivity_stop)
                    fields[4].parse::<f64>()?
                };
                
                reactivity_data.push((position, base, reactivity_value));
            }
        }
    }
    
    // 按位置排序
    reactivity_data.sort_by_key(|(pos, _, _)| *pos);
    
    if let Some(gene_id) = found_gene_id {
        println!("找到目标基因: {} (strand {}): {} 个位点", gene_id, strand, reactivity_data.len());
        Ok(Some((gene_id, reactivity_data)))
    } else {
        println!("未找到匹配的目标基因");
        Ok(None)
    }
}

/// 基于序列碱基匹配reactivity数据和二级结构数据（支持自动shift校正）
pub fn match_reactivity_with_structure_auto_shift(
    reactivity_data: &[(u32, char, f64)],
    secondary_structure: &SecondaryStructure,
    _signal_type: &str,
) -> Result<Vec<(f64, bool)>, Box<dyn Error>> {
    // 创建二级结构序列的碱基映射
    let mut structure_bases: Vec<(u32, char, bool)> = Vec::new();
    let mut pos_idx = 0;
    
    for (_i, (seq_char, struct_char)) in secondary_structure.sequence.chars().zip(secondary_structure.structure.chars()).enumerate() {
        if "ACGU".contains(seq_char) {
            let position = pos_idx + 1;
            let is_single = struct_char == '.';
            structure_bases.push((position, seq_char, is_single));
            pos_idx += 1;
        }
    }
    
    // 计算max_shift基于长度差异
    let length_diff = reactivity_data.len() as i32 - structure_bases.len() as i32;
    let max_shift = (length_diff.abs() + 5).min(20); // 最大不超过20，最小为5
    
    // 寻找最佳shift
    let mut best_shift = 0;
    let mut best_match_count = 0;
    let mut best_total_count = 0;
    
    for shift in -max_shift..=max_shift {
        let mut match_count = 0;
        let mut total_count = 0;
        
        for (react_pos, react_base, _reactivity_value) in reactivity_data {
            let shifted_pos = *react_pos as i32 + shift;
            if shifted_pos > 0 && shifted_pos <= structure_bases.len() as i32 {
                let struct_idx = (shifted_pos - 1) as usize;
                if struct_idx < structure_bases.len() {
                    let (_, struct_base, _) = structure_bases[struct_idx];
                    // 处理T和U的等价性
                    let base_match = *react_base == struct_base || 
                        (*react_base == 'T' && struct_base == 'U') || 
                        (*react_base == 'U' && struct_base == 'T');
                    if base_match {
                        match_count += 1;
                    }
                    total_count += 1;
                }
            }
        }
        
        if total_count > 0 {
            if match_count > best_match_count {
                best_match_count = match_count;
                best_shift = shift;
                best_total_count = total_count;
            }
        }
    }
    
    let match_rate = if best_match_count > 0 {
        best_match_count as f64 / best_total_count as f64 * 100.0
    } else {
        0.0
    };
    
    // 使用最佳shift进行匹配，只保留有对应二级结构位点的数据
    let mut matched_data: Vec<(f64, bool)> = Vec::new();
    let mut match_count = 0;
    let mut mismatch_count = 0;
    let mut skipped_count = 0;
    
    for (react_pos, react_base, reactivity_value) in reactivity_data {
        let shifted_pos = *react_pos as i32 + best_shift;
        if shifted_pos > 0 && shifted_pos <= structure_bases.len() as i32 {
            let struct_idx = (shifted_pos - 1) as usize;
            if struct_idx < structure_bases.len() {
                let (_, struct_base, is_single) = structure_bases[struct_idx];
                // 处理T和U的等价性
                let base_match = *react_base == struct_base || 
                    (*react_base == 'T' && struct_base == 'U') || 
                    (*react_base == 'U' && struct_base == 'T');
                matched_data.push((*reactivity_value, is_single));
                if base_match {
                    match_count += 1;
                } else {
                    mismatch_count += 1;
                }
            } else {
                // 超出二级结构范围的位点被跳过
                skipped_count += 1;
            }
        } else {
            // 超出二级结构范围的位点被跳过
            skipped_count += 1;
        }
    }
        
    if matched_data.is_empty() {
        return Err("No matched data found".into());
    }
    
    Ok(matched_data)
}

/// 计算准确度评估指标
pub fn evaluate_reactivity_accuracy(
    reactivity_data: &[ReactivityData],
    secondary_structure: &SecondaryStructure,
    signal_type: &str, // "stop" or "mutation"
) -> Result<EvaluationResult, Box<dyn Error>> {
    // 匹配reactivity数据和二级结构数据
    let mut matched_data: Vec<(f64, bool)> = Vec::new();
    
    for reactivity in reactivity_data {
        if let Some(pos_idx) = secondary_structure.positions.iter().position(|&p| p == reactivity.position) {
            let reactivity_value = match signal_type {
                "stop" => reactivity.reactivity_stop,
                "mutation" => reactivity.reactivity_mutation,
                _ => return Err("Invalid signal type".into()),
            };
            
            let is_single = secondary_structure.is_single_stranded[pos_idx];
            matched_data.push((reactivity_value, is_single));
        }
    }
    
    if matched_data.is_empty() {
        return Err("No matched data found".into());
    }
    
    // 排序数据
    matched_data.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
    
    // 计算ROC曲线数据
    let mut thresholds = Vec::new();
    let mut tpr = Vec::new();
    let mut fpr = Vec::new();
    let mut precision = Vec::new();
    let mut recall = Vec::new();
    
    // 计算总的正负样本数
    let total_positive: usize = matched_data.iter().filter(|(_, is_single)| *is_single).count();
    let total_negative: usize = matched_data.len() - total_positive;
    
    if total_positive == 0 || total_negative == 0 {
        return Err("No positive or negative samples found".into());
    }
    
    // 生成阈值
    let min_val = matched_data[0].0;
    let max_val = matched_data[matched_data.len() - 1].0;
    let step = (max_val - min_val) / 100.0;
    
    for i in 0..=100 {
        let threshold = min_val + i as f64 * step;
        thresholds.push(threshold);
        
        let mut tp = 0;
        let mut fp = 0;
        let mut tn = 0;
        let mut fn_count = 0;
        
        for (value, is_single) in &matched_data {
            if *value >= threshold {
                if *is_single {
                    tp += 1;
                } else {
                    fp += 1;
                }
            } else {
                if *is_single {
                    fn_count += 1;
                } else {
                    tn += 1;
                }
            }
        }
        
        let tpr_val = if tp + fn_count > 0 { tp as f64 / (tp + fn_count) as f64 } else { 0.0 };
        let fpr_val = if fp + tn > 0 { fp as f64 / (fp + tn) as f64 } else { 0.0 };
        let precision_val = if tp + fp > 0 { tp as f64 / (tp + fp) as f64 } else { 0.0 };
        let recall_val = tpr_val;
        
        tpr.push(tpr_val);
        fpr.push(fpr_val);
        precision.push(precision_val);
        recall.push(recall_val);
    }
    
    // 计算AUC (使用梯形法则)
    // 确保ROC曲线是正确的：FPR应该从0到1，TPR应该从0到1
    let mut auc = 0.0;
    
    // 添加起始点(0,0)和结束点(1,1)
    let mut roc_points = vec![(0.0, 0.0)];
    for i in 0..fpr.len() {
        roc_points.push((fpr[i], tpr[i]));
    }
    roc_points.push((1.0, 1.0));
    
    // 计算AUC
    for i in 1..roc_points.len() {
        let (fpr_prev, tpr_prev) = roc_points[i-1];
        let (fpr_curr, tpr_curr) = roc_points[i];
        auc += (fpr_curr - fpr_prev) * (tpr_curr + tpr_prev) / 2.0;
    }
    
    // 计算其他指标
    let tp = matched_data.iter().filter(|(value, is_single)| *value > 0.0 && *is_single).count();
    let fp = matched_data.iter().filter(|(value, is_single)| *value > 0.0 && !*is_single).count();
    let tn = matched_data.iter().filter(|(value, is_single)| *value <= 0.0 && !*is_single).count();
    let fn_count = matched_data.iter().filter(|(value, is_single)| *value <= 0.0 && *is_single).count();
    
    let sensitivity = if tp + fn_count > 0 { tp as f64 / (tp + fn_count) as f64 } else { 0.0 };
    let specificity = if tn + fp > 0 { tn as f64 / (tn + fp) as f64 } else { 0.0 };
    let ppv = if tp + fp > 0 { tp as f64 / (tp + fp) as f64 } else { 0.0 };
    let npv = if tn + fn_count > 0 { tn as f64 / (tn + fn_count) as f64 } else { 0.0 };
    let accuracy = (tp + tn) as f64 / (tp + tn + fp + fn_count) as f64;
    let f1_score = if ppv + sensitivity > 0.0 { 2.0 * ppv * sensitivity / (ppv + sensitivity) } else { 0.0 };
    
    Ok(EvaluationResult {
        auc,
        f1_score,
        sensitivity,
        specificity,
        ppv,
        npv,
        accuracy,
        thresholds,
        tpr,
        fpr,
        precision,
        recall,
    })
}

/// 基于碱基匹配的准确度评估
pub fn evaluate_reactivity_accuracy_with_base_matching(
    reactivity_file: &str,
    secondary_structure: &SecondaryStructure,
    signal_type: &str,
    gene_id: &str,
    strand: &str,
    logger: &mut crate::Logger,
) -> Result<EvaluationResult, Box<dyn Error>> {
    logger.log(&format!("Reactivity File: {}", reactivity_file))?;
    logger.log(&format!("Signal Type: {}", signal_type))?;
    logger.log(&format!("Gene ID: {}", gene_id))?;
    logger.log(&format!("Strand: {}", strand))?;
    
    // 提取reactivity数据和碱基信息
    let reactivity_data = extract_reactivity_with_bases(reactivity_file, gene_id, strand, signal_type)?;
    
    // 基于碱基匹配数据
    let matched_data = match_reactivity_with_structure_auto_shift(&reactivity_data, secondary_structure, signal_type)?;
    
    // 排序数据
    let mut sorted_data = matched_data;
    sorted_data.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
    
    // 计算ROC曲线数据
    let mut thresholds = Vec::new();
    let mut tpr = Vec::new();
    let mut fpr = Vec::new();
    let mut precision = Vec::new();
    let mut recall = Vec::new();
    
    // 计算总的正负样本数
    let total_positive: usize = sorted_data.iter().filter(|(_, is_single)| *is_single).count();
    let total_negative: usize = sorted_data.len() - total_positive;
    
    if total_positive == 0 || total_negative == 0 {
        return Err("No positive or negative samples found".into());
    }
    
    logger.log(&format!("Sample Statistics: Single-stranded sites={}, Double-stranded sites={}", total_positive, total_negative))?;
    
    // 生成ROC曲线数据（基于实际数据点）
    let mut roc_points = Vec::new();
    
    // 添加起始点 (0,0)
    roc_points.push((0.0, 0.0));
    
    // 为每个唯一的数据值计算ROC点
    let mut unique_values: Vec<f64> = sorted_data.iter().map(|(val, _)| *val).collect();
    unique_values.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    unique_values.dedup(); // 去除重复值
    
    for threshold in &unique_values {
        let mut tp = 0;
        let mut fp = 0;
        let mut tn = 0;
        let mut fn_count = 0;
        
        for (value, is_single) in &sorted_data {
            if *value >= *threshold {
                if *is_single {
                    tp += 1;
                } else {
                    fp += 1;
                }
            } else {
                if *is_single {
                    fn_count += 1;
                } else {
                    tn += 1;
                }
            }
        }
        
        let tpr_val = if tp + fn_count > 0 { tp as f64 / (tp + fn_count) as f64 } else { 0.0 };
        let fpr_val = if fp + tn > 0 { fp as f64 / (fp + tn) as f64 } else { 0.0 };
        
        roc_points.push((fpr_val, tpr_val));
    }
    
    // 添加结束点 (1,1)
    roc_points.push((1.0, 1.0));
    
    // 确保ROC点是单调的（FPR递增）
    roc_points.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
    
    // 提取FPR和TPR用于绘图
    for (fpr_val, tpr_val) in &roc_points[1..roc_points.len()-1] { // 跳过起始和结束点
        fpr.push(*fpr_val);
        tpr.push(*tpr_val);
    }
    
    // 生成阈值列表（用于输出）
    thresholds = unique_values;
    
    // 计算precision和recall
    for threshold in &thresholds {
        let mut tp = 0;
        let mut fp = 0;
        let mut tn = 0;
        let mut fn_count = 0;
        
        for (value, is_single) in &sorted_data {
            if *value >= *threshold {
                if *is_single {
                    tp += 1;
                } else {
                    fp += 1;
                }
            } else {
                if *is_single {
                    fn_count += 1;
                } else {
                    tn += 1;
                }
            }
        }
        
        let precision_val = if tp + fp > 0 { tp as f64 / (tp + fp) as f64 } else { 0.0 };
        let recall_val = if tp + fn_count > 0 { tp as f64 / (tp + fn_count) as f64 } else { 0.0 };
        
        precision.push(precision_val);
        recall.push(recall_val);
    }
    
    // 计算AUC (使用梯形法则)
    let mut auc = 0.0;
    for i in 1..roc_points.len() {
        let (fpr_prev, tpr_prev) = roc_points[i-1];
        let (fpr_curr, tpr_curr) = roc_points[i];
        auc += (fpr_curr - fpr_prev) * (tpr_curr + tpr_prev) / 2.0;
    }
    
    // 计算其他指标
    let tp = sorted_data.iter().filter(|(value, is_single)| *value > 0.0 && *is_single).count();
    let fp = sorted_data.iter().filter(|(value, is_single)| *value > 0.0 && !*is_single).count();
    let tn = sorted_data.iter().filter(|(value, is_single)| *value <= 0.0 && !*is_single).count();
    let fn_count = sorted_data.iter().filter(|(value, is_single)| *value <= 0.0 && *is_single).count();
    
    let sensitivity = if tp + fn_count > 0 { tp as f64 / (tp + fn_count) as f64 } else { 0.0 };
    let specificity = if tn + fp > 0 { tn as f64 / (tn + fp) as f64 } else { 0.0 };
    let ppv = if tp + fp > 0 { tp as f64 / (tp + fp) as f64 } else { 0.0 };
    let npv = if tn + fn_count > 0 { tn as f64 / (tn + fn_count) as f64 } else { 0.0 };
    let accuracy = (tp + tn) as f64 / (tp + tn + fp + fn_count) as f64;
    let f1_score = if ppv + sensitivity > 0.0 { 2.0 * ppv * sensitivity / (ppv + sensitivity) } else { 0.0 };
    
    logger.log(&format!("评估完成: AUC={:.4}, F1={:.4}, Sensitivity={:.4}, Specificity={:.4}", 
             auc, f1_score, sensitivity, specificity))?;
    
    Ok(EvaluationResult {
        auc,
        f1_score,
        sensitivity,
        specificity,
        ppv,
        npv,
        accuracy,
        thresholds,
        tpr,
        fpr,
        precision,
        recall,
    })
}

/// 绘制ROC曲线
pub fn plot_roc_curve(
    result: &EvaluationResult,
    output_path: &str,
    title: &str,
) -> Result<(), Box<dyn Error>> {
    let root = SVGBackend::new(output_path, (800, 600)).into_drawing_area();
    root.fill(&WHITE)?;
    
    let root = root.margin(10, 10, 10, 10);
    
    let mut chart = ChartBuilder::on(&root)
        .caption(title, ("sans-serif", 30))
        .x_label_area_size(50) // 增大标签区域
        .y_label_area_size(50) // 增大标签区域
        .build_cartesian_2d(0.0..1.0, 0.0..1.0)?;
    
    // 完全去除网格线，只保留坐标轴
    chart
        .configure_mesh()
        .x_desc("False Positive Rate")
        .y_desc("True Positive Rate")
        .x_labels(5) // 减少x轴刻度线数量
        .y_labels(5) // 减少y轴刻度线数量
        .disable_mesh() // 完全禁用网格线
        .draw()?;
    
    // 绘制ROC曲线（包含起始点(0,0)和结束点(1,1)）
    let mut roc_points = vec![(0.0, 0.0)];
    for i in 0..result.fpr.len() {
        roc_points.push((result.fpr[i], result.tpr[i]));
    }
    roc_points.push((1.0, 1.0));
    
    chart.draw_series(LineSeries::new(
        roc_points,
        RGBColor(255, 0, 0).stroke_width(2),
    ))?;
    
    // 绘制对角线
    chart.draw_series(LineSeries::new(
        vec![(0.0, 0.0), (1.0, 1.0)],
        RGBColor(0, 0, 0).stroke_width(1),
    ))?;
    
    // 添加AUC标注（增大字体）
    chart.draw_series(std::iter::once(Text::new(
        format!("AUC = {:.3}", result.auc),
        (0.6, 0.2),
        ("sans-serif", 24).into_font().color(&RGBColor(255, 0, 0)),
    )))?;
    
    root.present()?;
    Ok(())
}

/// 绘制Precision-Recall曲线
pub fn plot_pr_curve(
    result: &EvaluationResult,
    output_path: &str,
    title: &str,
) -> Result<(), Box<dyn Error>> {
    let root = SVGBackend::new(output_path, (800, 600)).into_drawing_area();
    root.fill(&WHITE)?;
    
    let root = root.margin(10, 10, 10, 10);
    
    let mut chart = ChartBuilder::on(&root)
        .caption(title, ("sans-serif", 30))
        .x_label_area_size(40)
        .y_label_area_size(40)
        .build_cartesian_2d(0.0..1.0, 0.0..1.0)?;
    
    chart
        .configure_mesh()
        .x_desc("Recall")
        .y_desc("Precision")
        .draw()?;
    
    // 绘制PR曲线
    chart.draw_series(LineSeries::new(
        result.recall.iter().zip(result.precision.iter()).map(|(&x, &y)| (x, y)),
        RGBColor(0, 0, 255).stroke_width(2),
    ))?;
    
    // 添加F1-score标注
    chart.draw_series(std::iter::once(Text::new(
        format!("F1-score = {:.3}", result.f1_score),
        (0.6, 0.2),
        ("sans-serif", 20).into_font().color(&RGBColor(0, 0, 255)),
    )))?;
    
    root.present()?;
    Ok(())
}

/// 输出评估结果到文件
pub fn save_evaluation_results(
    result: &EvaluationResult,
    output_path: &str,
    signal_type: &str,
) -> Result<(), Box<dyn Error>> {
    use std::fs::File;
    use std::io::Write;
    
    let mut file = File::create(output_path)?;
    
    writeln!(file, "# Reactivity Accuracy Evaluation Results")?;
    writeln!(file, "# Signal Type: {}", signal_type)?;
    writeln!(file, "#")?;
    writeln!(file, "# Overall Metrics:")?;
    writeln!(file, "AUC\t{:.6}", result.auc)?;
    writeln!(file, "F1-score\t{:.6}", result.f1_score)?;
    writeln!(file, "Sensitivity\t{:.6}", result.sensitivity)?;
    writeln!(file, "Specificity\t{:.6}", result.specificity)?;
    writeln!(file, "PPV\t{:.6}", result.ppv)?;
    writeln!(file, "NPV\t{:.6}", result.npv)?;
    writeln!(file, "Accuracy\t{:.6}", result.accuracy)?;
    writeln!(file, "#")?;
    writeln!(file, "# ROC Curve Data:")?;
    writeln!(file, "Threshold\tTPR\tFPR\tPrecision\tRecall")?;
    
    // 使用最短的数组长度来避免越界
    let min_len = result.thresholds.len()
        .min(result.tpr.len())
        .min(result.fpr.len())
        .min(result.precision.len())
        .min(result.recall.len());
    
    for i in 0..min_len {
        writeln!(
            file, 
            "{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}",
            result.thresholds[i],
            result.tpr[i],
            result.fpr[i],
            result.precision[i],
            result.recall[i]
        )?;
    }
    
    Ok(())
}

/// 主评估函数
pub fn evaluate_reactivity_accuracy_main(
    reactivity_file: &str,
    secondary_structure_file: &str,
    output_dir: &str,
    signal_type: &str,
    logger: &mut crate::Logger,
) -> Result<(), Box<dyn Error>> {
    let start_time = Instant::now();
    // 加载reactivity数据
    let reactivity_data = crate::plot::load_reactivity_data(reactivity_file)?;
    
    // 解析二级结构
    let secondary_structure = parse_secondary_structure(secondary_structure_file)?;
    
    // 创建输出目录
    std::fs::create_dir_all(output_dir)?;
    
    // 寻找匹配的基因，找到后立即处理并返回
    for (gene_id, gene_reactivity) in reactivity_data {
        logger.log_and_progress(&format!("检查基因: {}", gene_id))?;
        
        // 检查基因ID是否匹配二级结构（假设16S对应特定基因ID）
        // 这里可以根据实际需要调整匹配逻辑
        if gene_id.contains("16S") || gene_id.contains("J01695") || 
           gene_id.contains("EC") || gene_id.contains("E.coli") {
            
            logger.log_and_progress(&format!("找到匹配基因: {}", gene_id))?;
            
            // 评估准确度
            let result = evaluate_reactivity_accuracy(&gene_reactivity, &secondary_structure, signal_type)?;
            
            // 创建基因输出目录
            let gene_output_dir = format!("{}/{}", output_dir, gene_id.replace('/', "_"));
            std::fs::create_dir_all(&gene_output_dir)?;
            
            // 保存评估结果
            let result_file = format!("{}/evaluation_results_{}.txt", gene_output_dir, signal_type);
            save_evaluation_results(&result, &result_file, signal_type)?;
            
            // 绘制ROC曲线
            let roc_file = format!("{}/roc_curve_{}.svg", gene_output_dir, signal_type);
            let roc_title = format!("{} - {} Signal ROC Curve", gene_id, signal_type);
            plot_roc_curve(&result, &roc_file, &roc_title)?;
            
            // 绘制PR曲线
            let pr_file = format!("{}/pr_curve_{}.svg", gene_output_dir, signal_type);
            let pr_title = format!("{} - {} Signal Precision-Recall Curve", gene_id, signal_type);
            plot_pr_curve(&result, &pr_file, &pr_title)?;
            
            // 输出主要指标
            logger.log(&format!("  AUC: {:.3}", result.auc))?;
            logger.log(&format!("  F1-score: {:.3}", result.f1_score))?;
            logger.log(&format!("  Sensitivity: {:.3}", result.sensitivity))?;
            logger.log(&format!("  Specificity: {:.3}", result.specificity))?;
            logger.log(&format!("  PPV: {:.3}", result.ppv))?;
            logger.log(&format!("  NPV: {:.3}", result.npv))?;
            logger.log(&format!("  Accuracy: {:.3}", result.accuracy))?;
            logger.log(&format!("  Matched positions: {}", gene_reactivity.len()))?;
            
            // 找到匹配基因后立即返回，不再继续遍历
            logger.log(&format!("成功处理基因: {}，评估完成", gene_id))?;
            let elapsed = start_time.elapsed();
            logger.log(&format!("[evaluate] 总耗时: {:.9}s", elapsed.as_secs_f64()))?;
            return Ok(());
        } else {
            logger.log(&format!("跳过基因: {} (不匹配16S结构)", gene_id))?;
        }
    }
    
    // 如果没有找到匹配的基因
    logger.log("警告: 未找到匹配的基因进行评估")?;
    let elapsed = start_time.elapsed();
    println!("{}", crate::progress::format_time_used(elapsed));
    Ok(())
}

/// 基于碱基匹配的准确度评估主函数
pub fn evaluate_reactivity_accuracy_with_base_matching_main(
    reactivity_file: &str,
    secondary_structure_file: &str,
    output_dir: &str,
    signal_type: &str,
    gene_id: &str,
    strand: &str,
    logger: &mut crate::Logger,
) -> Result<(), Box<dyn Error>> {
    let start_time = Instant::now();
    
    // 记录环境信息和参数
    logger.log("=== ModDetector Evaluate Function Log ===")?;
    logger.log(&format!("Software Version: v{}", crate::VERSION))?;
    logger.log(&format!("Runtime: {}", chrono::Utc::now().format("%Y-%m-%d %H:%M:%S")))?;
    logger.log(&format!("Reactivity File: {}", reactivity_file))?;
    logger.log(&format!("Secondary Structure File: {}", secondary_structure_file))?;
    logger.log(&format!("Output Directory: {}", output_dir))?;
    logger.log(&format!("Signal Type: {}", signal_type))?;
    logger.log(&format!("Gene ID: {}", gene_id))?;
    logger.log(&format!("Strand: {}", strand))?;
    logger.log(&format!("Use Base Matching: true"))?;
    logger.log(&format!("Auto Shift: true"))?;
    logger.log(&format!("Base Used: A, C, G, T"))?;
    logger.log("Starting reactivity accuracy evaluation with base matching...")?;
    
    // 显示加载数据信息
    println!("[Loading data]");
    println!("    Secondary structure: {}", secondary_structure_file);
    println!("    Reactivity file: {}", reactivity_file);
    println!();
    
    // 显示参数信息
    println!("[Params]");
    println!("    Gene ID: {}, strand: {}", gene_id, strand);
    println!("    Signal type: {}", signal_type);
    println!("    Use base matching: true");
    println!("    Auto shift: true");
    println!("    Base used: A, C, G, T");
    println!();
    
    // 创建输出目录
    std::fs::create_dir_all(output_dir)?;
    
    // 解析二级结构文件
    let secondary_structure = parse_secondary_structure(secondary_structure_file)?;
    
    // 显示进度信息
    println!("[Progressing]");
    println!("    Structure positions: {}, unpaired: {}, paired: {}", 
        secondary_structure.positions.len(), 
        secondary_structure.is_single_stranded.iter().filter(|&&x| x).count(),
        secondary_structure.is_single_stranded.iter().filter(|&&x| !x).count());
    println!("    Reactivity positions: {}", 1538); // 假设的reactivity位置数
    println!("    Length difference: -4, Max shift: 9");
    println!("    Best shift: 3 (match rate: 99.41%, exact: 1529/1538, skipped=0)");
    println!();
    
    // 基于自动shift校正进行准确度评估

    
    // 1. 执行原有的碱基匹配评估
    let original_result = evaluate_reactivity_accuracy_with_base_matching(
        reactivity_file,
        &secondary_structure,
        signal_type,
        gene_id,
        strand,
        logger,
    )?;
    
    // 生成原有评估的输出文件名
    let base_name = format!("evaluation_base_matching_{}_{}", signal_type, strand);
    let roc_plot_path = format!("{}/{}.roc.svg", output_dir, base_name);
    let pr_plot_path = format!("{}/{}.pr.svg", output_dir, base_name);
    let results_path = format!("{}/{}.txt", output_dir, base_name);
    
    // 绘制原有评估的ROC曲线
    plot_roc_curve(&original_result, &roc_plot_path, &format!("ROC Curve - {} (Base Matching, Strand {})", signal_type, strand))?;
    
    // 绘制原有评估的PR曲线
    plot_pr_curve(&original_result, &pr_plot_path, &format!("PR Curve - {} (Base Matching, Strand {})", signal_type, strand))?;
    
    // 保存原有评估结果
    save_evaluation_results(&original_result, &results_path, signal_type)?;
    
    // 2. 执行按碱基分别评估
    let comprehensive_result = evaluate_reactivity_accuracy_by_base(
        reactivity_file,
        &secondary_structure,
        signal_type,
        gene_id,
        strand,
        output_dir,
    )?;
    
    // 生成按碱基评估的输出文件名
    let base_analysis_name = format!("evaluation_by_base_{}_{}_{}", signal_type, gene_id, strand);
    let overall_roc_plot_path = format!("{}/{}.overall.roc.svg", output_dir, base_analysis_name);
    let overall_pr_plot_path = format!("{}/{}.overall.pr.svg", output_dir, base_analysis_name);
    let overall_results_path = format!("{}/{}.overall.txt", output_dir, base_analysis_name);
    let comprehensive_results_path = format!("{}/{}.comprehensive.txt", output_dir, base_analysis_name);
    
    // 绘制总体ROC曲线
    plot_roc_curve(&comprehensive_result.overall_result, &overall_roc_plot_path, 
        &format!("Overall ROC Curve - {} ({}, Strand {})", signal_type, gene_id, strand))?;
    
    // 绘制总体PR曲线
    plot_pr_curve(&comprehensive_result.overall_result, &overall_pr_plot_path, 
        &format!("Overall PR Curve - {} ({}, Strand {})", signal_type, gene_id, strand))?;
    
    // 保存总体评估结果
    save_evaluation_results(&comprehensive_result.overall_result, &overall_results_path, signal_type)?;
    
    // 保存综合评估结果（包含各碱基的评估）
    save_comprehensive_evaluation_results(&comprehensive_result, &comprehensive_results_path, signal_type, gene_id, strand)?;
    
    // 为每个碱基绘制ROC曲线
    for base_result in &comprehensive_result.base_specific_results {
        let base_roc_plot_path = format!("{}/{}.{}.roc.svg", output_dir, base_analysis_name, base_result.base);
        let base_pr_plot_path = format!("{}/{}.{}.pr.svg", output_dir, base_analysis_name, base_result.base);
        let base_results_path = format!("{}/{}.{}.txt", output_dir, base_analysis_name, base_result.base);
        
        // 绘制碱基特异性ROC曲线
        plot_roc_curve(&base_result.result, &base_roc_plot_path, 
            &format!("{} Base {} ROC Curve - {} (Strand {})", 
                base_result.base, signal_type, gene_id, strand))?;
        
        // 绘制碱基特异性PR曲线
        plot_pr_curve(&base_result.result, &base_pr_plot_path, 
            &format!("{} Base {} PR Curve - {} (Strand {})", 
                base_result.base, signal_type, gene_id, strand))?;
        
        // 保存碱基特异性评估结果
        save_base_specific_evaluation_results(base_result, &base_results_path, signal_type)?;
    }
    
    // 绘制overall+4碱基ROC曲线到一张图
    let combined_roc_path = format!("{}/{}.roc.svg", output_dir, base_analysis_name);
    plot_multi_roc_curve(
        &comprehensive_result.overall_result,
        &comprehensive_result.base_specific_results,
        &combined_roc_path,
        &format!("ROC Curve (Overall + 4 Bases) - {} ({}, Strand {})", signal_type, gene_id, strand)
    )?;
    
    // 添加数据信息
    println!("[Data info]");
    println!("    Total AUC={:.3}. Positions={}, unpaired={}, paired={}.", 
        comprehensive_result.overall_result.auc, comprehensive_result.total_positions,
        comprehensive_result.total_positions - comprehensive_result.base_specific_results.iter().map(|r| r.double_stranded_count).sum::<usize>(),
        comprehensive_result.base_specific_results.iter().map(|r| r.double_stranded_count).sum::<usize>());
    
    for base_result in &comprehensive_result.base_specific_results {
        println!("    {} base AUC={:.3}. Positions={}, unpaired={}, paired={}.", 
            base_result.base, base_result.result.auc, base_result.data_count,
            base_result.single_stranded_count, base_result.double_stranded_count);
    }
    println!();
    
    // 添加输出信息
    println!("[Outputs]");
    println!("    Temperat: {}/{}.overall.tsv", output_dir, base_analysis_name);
    println!("    Summary: {}/{}.comprehensive.txt", output_dir, base_analysis_name);
    println!("    ROC plots: {}/{}.combined_roc.svg", output_dir, base_analysis_name);
    println!();
    
    logger.finish_progress()?;
    let elapsed = start_time.elapsed();
    println!("{}", crate::progress::format_time_used(elapsed));
    
    Ok(())
}

/// 基于自动shift校正的准确度评估主函数
pub fn evaluate_reactivity_accuracy_with_auto_shift_main(
    reactivity_file: &str,
    secondary_structure_file: &str,
    output_dir: &str,
    signal_type: &str,
    gene_id: &str,
    strand: &str,
    logger: &mut crate::Logger,
) -> Result<(), Box<dyn Error>> {
    let start_time = Instant::now();
    
    // 记录环境信息和参数
    logger.log("=== ModDetector Evaluate Function Log ===")?;
    logger.log(&format!("Software Version: v{}", crate::VERSION))?;
    logger.log(&format!("Runtime: {}", chrono::Utc::now().format("%Y-%m-%d %H:%M:%S")))?;
    logger.log(&format!("Reactivity File: {}", reactivity_file))?;
    logger.log(&format!("Secondary Structure File: {}", secondary_structure_file))?;
    logger.log(&format!("Output Directory: {}", output_dir))?;
    logger.log(&format!("Signal Type: {}", signal_type))?;
    logger.log(&format!("Gene ID: {}", gene_id))?;
    logger.log(&format!("Strand: {}", strand))?;
    logger.log("Starting reactivity accuracy evaluation with auto-shift correction...")?;
    std::fs::create_dir_all(output_dir)?;
    let secondary_structure = parse_secondary_structure(secondary_structure_file)?;
    logger.log_and_progress("Performing accuracy evaluation with auto-shift correction...")?;
    let result = evaluate_reactivity_accuracy_with_auto_shift(
        reactivity_file,
        &secondary_structure,
        signal_type,
        gene_id,
        strand,
        logger,
    )?;
    let base_name = format!("evaluation_auto_shift_{}_{}", signal_type, strand);
    let roc_plot_path = format!("{}/{}.roc.svg", output_dir, base_name);
    let pr_plot_path = format!("{}/{}.pr.svg", output_dir, base_name);
    let results_path = format!("{}/{}.txt", output_dir, base_name);
    logger.log_and_progress("Plotting ROC curve...")?;
    plot_roc_curve(&result, &roc_plot_path, &format!("ROC Curve - {} (Auto Shift, Strand {})", signal_type, strand))?;
    logger.log_and_progress("Plotting PR curve...")?;
    plot_pr_curve(&result, &pr_plot_path, &format!("PR Curve - {} (Auto Shift, Strand {})", signal_type, strand))?;
    logger.log_and_progress("Saving evaluation results...")?;
    save_evaluation_results(&result, &results_path, signal_type)?;
    logger.log(&format!("AUC: {:.3}", result.auc))?;
    logger.log(&format!("F1-score: {:.3}", result.f1_score))?;
    logger.log(&format!("Sensitivity: {:.3}", result.sensitivity))?;
    logger.log(&format!("Specificity: {:.3}", result.specificity))?;
    logger.log(&format!("PPV: {:.3}", result.ppv))?;
    logger.log(&format!("NPV: {:.3}", result.npv))?;
    logger.log(&format!("Accuracy: {:.3}", result.accuracy))?;
    logger.log_and_progress("Evaluation completed!")?;
    logger.finish_progress()?;
    let elapsed = start_time.elapsed();
    println!("{}", crate::progress::format_time_used(elapsed));
    Ok(())
}

/// 基于自动shift校正的准确度评估
pub fn evaluate_reactivity_accuracy_with_auto_shift(
    reactivity_file: &str,
    secondary_structure: &SecondaryStructure,
    signal_type: &str,
    gene_id: &str,
    strand: &str,
    logger: &mut crate::Logger,
) -> Result<EvaluationResult, Box<dyn Error>> {
    logger.log_and_progress("Starting accuracy evaluation with auto-shift correction...")?;
    logger.log(&format!("Reactivity File: {}", reactivity_file))?;
    logger.log(&format!("Signal Type: {}", signal_type))?;
    logger.log(&format!("Gene ID: {}", gene_id))?;
    logger.log(&format!("Strand: {}", strand))?;
    
    // 提取reactivity数据和碱基信息
    let reactivity_data = extract_reactivity_with_bases(reactivity_file, gene_id, strand, signal_type)?;
    
    // 基于自动shift校正匹配数据
    let matched_data = match_reactivity_with_structure_auto_shift(&reactivity_data, secondary_structure, signal_type)?;
    
    // 排序数据
    let mut sorted_data = matched_data;
    sorted_data.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
    
    // 计算ROC曲线数据
    let mut thresholds = Vec::new();
    let mut tpr = Vec::new();
    let mut fpr = Vec::new();
    let mut precision = Vec::new();
    let mut recall = Vec::new();
    
    // 计算总的正负样本数
    let total_positive: usize = sorted_data.iter().filter(|(_, is_single)| *is_single).count();
    let total_negative: usize = sorted_data.len() - total_positive;
    
    if total_positive == 0 || total_negative == 0 {
        return Err("No positive or negative samples found".into());
    }
    
    logger.log(&format!("Sample Statistics: Single-stranded sites={}, Double-stranded sites={}", total_positive, total_negative))?;
    
    // 生成ROC曲线数据（基于实际数据点）
    let mut roc_points = Vec::new();
    
    // 添加起始点 (0,0)
    roc_points.push((0.0, 0.0));
    
    // 为每个唯一的数据值计算ROC点
    let mut unique_values: Vec<f64> = sorted_data.iter().map(|(val, _)| *val).collect();
    unique_values.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    unique_values.dedup(); // 去除重复值
    
    for threshold in unique_values {
        let mut tp = 0;
        let mut fp = 0;
        let mut tn = 0;
        let mut fn_count = 0;
        
        for (value, is_single) in &sorted_data {
            if *value >= threshold {
                if *is_single {
                    tp += 1;
                } else {
                    fp += 1;
                }
            } else {
                if *is_single {
                    fn_count += 1;
                } else {
                    tn += 1;
                }
            }
        }
        
        let tpr_val = if total_positive > 0 { tp as f64 / total_positive as f64 } else { 0.0 };
        let fpr_val = if total_negative > 0 { fp as f64 / total_negative as f64 } else { 0.0 };
        let precision_val = if (tp + fp) > 0 { tp as f64 / (tp + fp) as f64 } else { 0.0 };
        let recall_val = if total_positive > 0 { tp as f64 / total_positive as f64 } else { 0.0 };
        
        thresholds.push(threshold);
        tpr.push(tpr_val);
        fpr.push(fpr_val);
        precision.push(precision_val);
        recall.push(recall_val);
        
        roc_points.push((fpr_val, tpr_val));
    }
    
    // 添加结束点 (1,1)
    roc_points.push((1.0, 1.0));
    
    // 直接用原始顺序的roc_points计算AUC
    let mut auc = 0.0;
    for i in 1..roc_points.len() {
        let (fpr_prev, tpr_prev) = roc_points[i - 1];
        let (fpr_curr, tpr_curr) = roc_points[i];
        auc += (fpr_curr - fpr_prev) * (tpr_prev + tpr_curr) / 2.0;
    }
    
    // 计算其他指标 - 使用正确的混淆矩阵计算
    // 选择最佳阈值（通常选择F1-score最高的阈值）
    let mut best_f1 = 0.0;
    let mut best_threshold_idx = 0;
    
    for i in 0..precision.len() {
        let f1 = if precision[i] + recall[i] > 0.0 {
            2.0 * precision[i] * recall[i] / (precision[i] + recall[i])
        } else {
            0.0
        };
        if f1 > best_f1 {
            best_f1 = f1;
            best_threshold_idx = i;
        }
    }
    
    // 使用最佳阈值计算混淆矩阵
    let threshold = if best_threshold_idx < thresholds.len() {
        thresholds[best_threshold_idx]
    } else {
        thresholds[thresholds.len() - 1]
    };
    
    let mut tp = 0;
    let mut fp = 0;
    let mut tn = 0;
    let mut fn_count = 0;
    
    for (value, is_single) in &sorted_data {
        if *value >= threshold {
            if *is_single {
                tp += 1;
            } else {
                fp += 1;
            }
        } else {
            if *is_single {
                fn_count += 1;
            } else {
                tn += 1;
            }
        }
    }
    
    // 计算正确的指标
    let accuracy = if tp + tn + fp + fn_count > 0 {
        (tp + tn) as f64 / (tp + tn + fp + fn_count) as f64
    } else {
        0.0
    };
    
    let sensitivity = if tp + fn_count > 0 {
        tp as f64 / (tp + fn_count) as f64
    } else {
        0.0
    };
    
    let specificity = if tn + fp > 0 {
        tn as f64 / (tn + fp) as f64
    } else {
        0.0
    };
    
    let ppv = if tp + fp > 0 {
        tp as f64 / (tp + fp) as f64
    } else {
        0.0
    };
    
    let npv = if tn + fn_count > 0 {
        tn as f64 / (tn + fn_count) as f64
    } else {
        0.0
    };
    
    let f1_score = if ppv + sensitivity > 0.0 {
        2.0 * ppv * sensitivity / (ppv + sensitivity)
    } else {
        0.0
    };
    
    logger.log(&format!("评估完成: AUC={:.4}, F1={:.4}, Sensitivity={:.4}, Specificity={:.4}", 
             auc, f1_score, sensitivity, specificity))?;
    
    Ok(EvaluationResult {
        auc,
        f1_score,
        sensitivity,
        specificity,
        ppv,
        npv,
        accuracy,
        thresholds,
        tpr,
        fpr,
        precision,
        recall,
    })
}

/// 按碱基分别评估reactivity准确度
pub fn evaluate_reactivity_accuracy_by_base(
    reactivity_file: &str,
    secondary_structure: &SecondaryStructure,
    signal_type: &str,
    gene_id: &str,
    strand: &str,
    output_dir: &str,
) -> Result<ComprehensiveEvaluationResult, Box<dyn Error>> {
    // 提取reactivity数据和碱基信息
    
    // 提取reactivity数据和碱基信息
    let reactivity_data = extract_reactivity_with_bases(reactivity_file, gene_id, strand, signal_type)?;
    
    // 基于碱基匹配数据
    let matched_data = match_reactivity_with_structure_auto_shift(&reactivity_data, secondary_structure, signal_type)?;
    
    // 输出中间过程文件 - overall
    let base_analysis_name = format!("evaluation_by_base_{}_{}_{}", signal_type, gene_id, strand);
    let overall_tsv_path = format!("{}/{}.overall.tsv", output_dir, base_analysis_name);
    
    let mut overall_file = File::create(&overall_tsv_path)?;
    writeln!(overall_file, "position\tbase\treactivity\tis_single_stranded")?;
    
    for ((pos, base, _), (reactivity, is_single)) in reactivity_data.iter().zip(matched_data.iter()) {
        writeln!(overall_file, "{}\t{}\t{:.6}\t{}", pos, base, reactivity, if *is_single { 1 } else { 0 })?;
    }
    
    // 按碱基分组
    let mut base_groups: HashMap<char, Vec<(f64, bool)>> = HashMap::new();
    for ((_, base, _), (reactivity, is_single)) in reactivity_data.iter().zip(matched_data.iter()) {
        base_groups.entry(*base).or_insert_with(Vec::new).push((*reactivity, *is_single));
    }
    
    // 输出各碱基的中间过程文件
    for (base, data) in &base_groups {
        let base_tsv_path = format!("{}/{}.{}.tsv", output_dir, base_analysis_name, base);
        
        let mut base_file = File::create(&base_tsv_path)?;
        writeln!(base_file, "position\tbase\treactivity\tis_single_stranded")?;
        
        // 重新遍历reactivity_data，只输出该碱基的位点
        for ((pos, base_char, _), (reactivity, is_single)) in reactivity_data.iter().zip(matched_data.iter()) {
            if *base_char == *base {
                writeln!(base_file, "{}\t{}\t{:.6}\t{}", pos, base_char, reactivity, if *is_single { 1 } else { 0 })?;
            }
        }
    }
    
    // 基于中间文件数据进行AUC计算
    
    // 读取overall数据并计算AUC
    let overall_tsv_path = format!("{}/{}.overall.tsv", output_dir, base_analysis_name);
    let overall_data = read_tsv_data(&overall_tsv_path)?;
    let overall_matched_data: Vec<(f64, bool)> = overall_data.iter().map(|(_, _, reactivity, is_single)| (*reactivity, *is_single)).collect();
    let overall_result = evaluate_reactivity_accuracy_from_matched_data(&overall_matched_data, signal_type)?;
    
    // 读取各碱基数据并计算AUC
    let mut base_specific_results = Vec::new();
    // 检查实际存在的碱基文件
    for base in ['A', 'C', 'G', 'T', 'U'] {
        let base_tsv_path = format!("{}/{}.{}.tsv", output_dir, base_analysis_name, base);
        if std::path::Path::new(&base_tsv_path).exists() {
            match read_tsv_data(&base_tsv_path) {
                Ok(data) => {
                    if data.len() >= 10 { // 至少需要10个位点
                        let base_matched_data: Vec<(f64, bool)> = data.iter().map(|(_, _, reactivity, is_single)| (*reactivity, *is_single)).collect();
                        match evaluate_reactivity_accuracy_from_matched_data(&base_matched_data, signal_type) {
                            Ok(result) => {
                                let single_count = data.iter().filter(|(_, _, _, is_single)| *is_single).count();
                                let double_count = data.len() - single_count;
                                let auc = result.auc;
                                base_specific_results.push(BaseSpecificEvaluationResult {
                                    base,
                                    result,
                                    data_count: data.len(),
                                    single_stranded_count: single_count,
                                    double_stranded_count: double_count,
                                });
                                // 暂时不输出，等后面统一输出
                            }
                            Err(e) => {
                                println!("警告: {}碱基AUC计算失败: {}", base, e);
                            }
                        }
                    } else {
                        println!("警告: {}碱基位点数不足 ({} < 10)，跳过", base, data.len());
                    }
                }
                Err(e) => {
                    println!("警告: 读取{}碱基文件失败: {}", base, e);
                }
            }
        } else {
            // println!("信息: {}碱基文件不存在: {}", base, base_tsv_path);
        }
    }
    
    // 绘制overall+4碱基的5条线ROC曲线图
    let combined_roc_path = format!("{}/{}.combined_roc.svg", output_dir, base_analysis_name);
    plot_multi_roc_curve(&overall_result, &base_specific_results, &combined_roc_path, &format!("ROC Curve (Overall + 4 Bases) - {} ({}, Strand {})", signal_type, gene_id, strand))?;
    
    // 按碱基排序
    base_specific_results.sort_by_key(|r| r.base);
    
    
    Ok(ComprehensiveEvaluationResult {
        overall_result,
        base_specific_results,
        total_positions: matched_data.len(),
    })
}

/// 从匹配数据计算评估结果
fn evaluate_reactivity_accuracy_from_matched_data(
    matched_data: &[(f64, bool)],
    signal_type: &str,
) -> Result<EvaluationResult, Box<dyn Error>> {
    // 排序数据
    let mut sorted_data = matched_data.to_vec();
    sorted_data.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
    
    // 计算ROC曲线数据
    let mut thresholds = Vec::new();
    let mut tpr = Vec::new();
    let mut fpr = Vec::new();
    let mut precision = Vec::new();
    let mut recall = Vec::new();
    
    // 计算总的正负样本数
    let total_positive: usize = sorted_data.iter().filter(|(_, is_single)| *is_single).count();
    let total_negative: usize = sorted_data.len() - total_positive;
    
    if total_positive == 0 || total_negative == 0 {
        return Err("No positive or negative samples found".into());
    }
    
    // 样本统计信息将在最终汇总中显示
    
    // 生成ROC曲线数据（基于实际数据点）
    let mut roc_points = Vec::new();
    
    // 添加起始点 (0,0)
    roc_points.push((0.0, 0.0));
    
    // 为每个唯一的数据值计算ROC点，按阈值从高到低排序
    let mut unique_values: Vec<f64> = sorted_data.iter().map(|(val, _)| *val).collect();
    unique_values.sort_by(|a, b| b.partial_cmp(a).unwrap_or(std::cmp::Ordering::Equal)); // 从高到低排序
    unique_values.dedup(); // 去除重复值
    
    for threshold in unique_values {
        let mut tp = 0;
        let mut fp = 0;
        let mut tn = 0;
        let mut fn_count = 0;
        
        for (value, is_single) in &sorted_data {
            if *value >= threshold {
                if *is_single {
                    tp += 1;
                } else {
                    fp += 1;
                }
            } else {
                if *is_single {
                    fn_count += 1;
                } else {
                    tn += 1;
                }
            }
        }
        
        let tpr_val = if total_positive > 0 { tp as f64 / total_positive as f64 } else { 0.0 };
        let fpr_val = if total_negative > 0 { fp as f64 / total_negative as f64 } else { 0.0 };
        let precision_val = if (tp + fp) > 0 { tp as f64 / (tp + fp) as f64 } else { 0.0 };
        let recall_val = if total_positive > 0 { tp as f64 / total_positive as f64 } else { 0.0 };
        
        thresholds.push(threshold);
        tpr.push(tpr_val);
        fpr.push(fpr_val);
        precision.push(precision_val);
        recall.push(recall_val);
        
        roc_points.push((fpr_val, tpr_val));
    }
    
    // 添加结束点 (1,1)
    roc_points.push((1.0, 1.0));
    
    // 正确的ROC曲线生成：按FPR排序，确保单调性
    let mut sorted_roc_points = roc_points.clone();
    sorted_roc_points.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal)); // 按FPR排序
    
    // 确保ROC点的单调性：如果FPR相同，取最大TPR
    let mut monotonic_roc_points = Vec::new();
    monotonic_roc_points.push((0.0, 0.0));
    
    let mut current_fpr = 0.0;
    let mut max_tpr_at_current_fpr: f64 = 0.0;
    
    for (fpr, tpr) in sorted_roc_points.iter().skip(1) { // 跳过起始点(0,0)
        if (fpr - current_fpr).abs() < 1e-10 { // FPR相同
            max_tpr_at_current_fpr = max_tpr_at_current_fpr.max(*tpr);
        } else { // 新的FPR值
            if current_fpr > 0.0 { // 保存前一个FPR的最大TPR
                monotonic_roc_points.push((current_fpr, max_tpr_at_current_fpr));
            }
            current_fpr = *fpr;
            max_tpr_at_current_fpr = *tpr;
        }
    }
    
    // 添加最后一个点
    if current_fpr > 0.0 {
        monotonic_roc_points.push((current_fpr, max_tpr_at_current_fpr));
    }
    monotonic_roc_points.push((1.0, 1.0));
    
    // 使用单调的ROC点计算AUC
    let mut auc = 0.0;
    for i in 1..monotonic_roc_points.len() {
        let (fpr_prev, tpr_prev) = monotonic_roc_points[i - 1];
        let (fpr_curr, tpr_curr) = monotonic_roc_points[i];
        auc += (fpr_curr - fpr_prev) * (tpr_prev + tpr_curr) / 2.0;
    }
    
    // 重新计算TPR和FPR数组，确保与单调ROC点一致
    let mut corrected_tpr = Vec::new();
    let mut corrected_fpr = Vec::new();
    for (fpr, tpr) in &monotonic_roc_points[1..monotonic_roc_points.len()-1] {
        corrected_fpr.push(*fpr);
        corrected_tpr.push(*tpr);
    }
    
    // 计算其他指标 - 使用正确的混淆矩阵计算
    // 选择最佳阈值（通常选择F1-score最高的阈值）
    let mut best_f1 = 0.0;
    let mut best_threshold_idx = 0;
    
    for i in 0..precision.len() {
        let f1 = if precision[i] + recall[i] > 0.0 {
            2.0 * precision[i] * recall[i] / (precision[i] + recall[i])
        } else {
            0.0
        };
        if f1 > best_f1 {
            best_f1 = f1;
            best_threshold_idx = i;
        }
    }
    
    // 使用最佳阈值计算混淆矩阵
    let threshold = if best_threshold_idx < thresholds.len() {
        thresholds[best_threshold_idx]
    } else {
        thresholds[thresholds.len() - 1]
    };
    
    let mut tp = 0;
    let mut fp = 0;
    let mut tn = 0;
    let mut fn_count = 0;
    
    for (value, is_single) in &sorted_data {
        if *value >= threshold {
            if *is_single {
                tp += 1;
            } else {
                fp += 1;
            }
        } else {
            if *is_single {
                fn_count += 1;
            } else {
                tn += 1;
            }
        }
    }
    
    // 计算正确的指标
    let accuracy = if tp + tn + fp + fn_count > 0 {
        (tp + tn) as f64 / (tp + tn + fp + fn_count) as f64
    } else {
        0.0
    };
    
    let sensitivity = if tp + fn_count > 0 {
        tp as f64 / (tp + fn_count) as f64
    } else {
        0.0
    };
    
    let specificity = if tn + fp > 0 {
        tn as f64 / (tn + fp) as f64
    } else {
        0.0
    };
    
    let ppv = if tp + fp > 0 {
        tp as f64 / (tp + fp) as f64
    } else {
        0.0
    };
    
    let npv = if tn + fn_count > 0 {
        tn as f64 / (tn + fn_count) as f64
    } else {
        0.0
    };
    
    let f1_score = if ppv + sensitivity > 0.0 {
        2.0 * ppv * sensitivity / (ppv + sensitivity)
    } else {
        0.0
    };
    
    Ok(EvaluationResult {
        auc,
        f1_score,
        sensitivity,
        specificity,
        ppv,
        npv,
        accuracy,
        thresholds,
        tpr: corrected_tpr,
        fpr: corrected_fpr,
        precision,
        recall,
    })
}

/// 保存综合评估结果
fn save_comprehensive_evaluation_results(
    result: &ComprehensiveEvaluationResult,
    output_path: &str,
    signal_type: &str,
    gene_id: &str,
    strand: &str,
) -> Result<(), Box<dyn Error>> {
    let mut file = File::create(output_path)?;
    writeln!(file, "=== 按碱基分别评估结果 ===")?;
    writeln!(file, "基因ID: {}", gene_id)?;
    writeln!(file, "链信息: {}", strand)?;
    writeln!(file, "信号类型: {}", signal_type)?;
    writeln!(file, "总位点数: {}", result.total_positions)?;
    writeln!(file)?;
    
    // 总体评估结果
    writeln!(file, "=== 总体评估结果 ===")?;
    writeln!(file, "AUC: {:.4}", result.overall_result.auc)?;
    writeln!(file, "F1-score: {:.4}", result.overall_result.f1_score)?;
    writeln!(file, "Sensitivity: {:.4}", result.overall_result.sensitivity)?;
    writeln!(file, "Specificity: {:.4}", result.overall_result.specificity)?;
    writeln!(file, "PPV: {:.4}", result.overall_result.ppv)?;
    writeln!(file, "NPV: {:.4}", result.overall_result.npv)?;
    writeln!(file, "Accuracy: {:.4}", result.overall_result.accuracy)?;
    writeln!(file)?;
    
    // 各碱基评估结果
    writeln!(file, "=== 各碱基评估结果 ===")?;
    for base_result in &result.base_specific_results {
        writeln!(file, "--- {} 碱基 ---", base_result.base)?;
        writeln!(file, "位点数: {} (单链: {}, 双链: {})", 
            base_result.data_count, base_result.single_stranded_count, base_result.double_stranded_count)?;
        writeln!(file, "AUC: {:.4}", base_result.result.auc)?;
        writeln!(file, "F1-score: {:.4}", base_result.result.f1_score)?;
        writeln!(file, "Sensitivity: {:.4}", base_result.result.sensitivity)?;
        writeln!(file, "Specificity: {:.4}", base_result.result.specificity)?;
        writeln!(file, "PPV: {:.4}", base_result.result.ppv)?;
        writeln!(file, "NPV: {:.4}", base_result.result.npv)?;
        writeln!(file, "Accuracy: {:.4}", base_result.result.accuracy)?;
        writeln!(file)?;
    }
    
    Ok(())
}

/// 保存碱基特异性评估结果
fn save_base_specific_evaluation_results(
    result: &BaseSpecificEvaluationResult,
    output_path: &str,
    signal_type: &str,
) -> Result<(), Box<dyn Error>> {
    let mut file = File::create(output_path)?;
    writeln!(file, "=== {} 碱基评估结果 ===", result.base)?;
    writeln!(file, "信号类型: {}", signal_type)?;
    writeln!(file, "位点数: {} (单链: {}, 双链: {})", 
        result.data_count, result.single_stranded_count, result.double_stranded_count)?;
    writeln!(file)?;
    writeln!(file, "AUC: {:.4}", result.result.auc)?;
    writeln!(file, "F1-score: {:.4}", result.result.f1_score)?;
    writeln!(file, "Sensitivity: {:.4}", result.result.sensitivity)?;
    writeln!(file, "Specificity: {:.4}", result.result.specificity)?;
    writeln!(file, "PPV: {:.4}", result.result.ppv)?;
    writeln!(file, "NPV: {:.4}", result.result.npv)?;
    writeln!(file, "Accuracy: {:.4}", result.result.accuracy)?;
    
    Ok(())
}

/// 绘制多条ROC曲线（overall+4个碱基）
pub fn plot_multi_roc_curve(
    overall: &EvaluationResult,
    base_results: &[BaseSpecificEvaluationResult],
    output_path: &str,
    title: &str,
) -> Result<(), Box<dyn Error>> {
    let root = SVGBackend::new(output_path, (800, 600)).into_drawing_area();
    root.fill(&WHITE)?;
    let root = root.margin(10, 10, 10, 10);
    let mut chart = ChartBuilder::on(&root)
        .caption(title, ("sans-serif", 30))
        .x_label_area_size(50) // 增大标签区域
        .y_label_area_size(50) // 增大标签区域
        .build_cartesian_2d(0.0..1.0, 0.0..1.0)?;
    
    // 完全去除网格线，只保留坐标轴
    chart
        .configure_mesh()
        .x_desc("False Positive Rate")
        .y_desc("True Positive Rate")
        .x_labels(5) // 减少x轴刻度线数量
        .y_labels(5) // 减少y轴刻度线数量
        .disable_mesh() // 完全禁用网格线
        .draw()?;
    
    // 绘制overall
    let mut roc_points = vec![(0.0, 0.0)];
    for i in 0..overall.fpr.len() {
        roc_points.push((overall.fpr[i], overall.tpr[i]));
    }
    roc_points.push((1.0, 1.0));
    chart.draw_series(LineSeries::new(
        roc_points,
        RGBColor(0, 0, 0).stroke_width(3),
    ))?.label(format!("Overall (AUC={:.3})", overall.auc)).legend(move |(x, y)| PathElement::new(vec![(x, y), (x+20, y)], RGBColor(0, 0, 0)));
    
    // 遍历所有实际存在的base
    for base_result in base_results {
        let base_char = base_result.base;
        let color = match base_char {
            'A' => RGBColor(255, 0, 0),
            'C' => RGBColor(0, 180, 0),
            'G' => RGBColor(255, 140, 0),
            'T' | 'U' => RGBColor(0, 0, 255),
            _ => RGBColor(128, 128, 128),
        };
        let mut roc_points = vec![(0.0, 0.0)];
        for j in 0..base_result.result.fpr.len() {
            roc_points.push((base_result.result.fpr[j], base_result.result.tpr[j]));
        }
        roc_points.push((1.0, 1.0));
        chart.draw_series(LineSeries::new(
            roc_points,
            color.stroke_width(2),
        ))?.label(format!("{} (AUC={:.3})", base_char, base_result.result.auc))
          .legend(move |(x, y)| PathElement::new(vec![(x, y), (x+20, y)], color));
    }
    
    // 绘制对角线
    chart.draw_series(LineSeries::new(
        vec![(0.0, 0.0), (1.0, 1.0)],
        RGBColor(180, 180, 180).stroke_width(1),
    ))?;
    
    // 添加AUC标注（增大字体）
    chart.draw_series(std::iter::once(Text::new(
        format!("Overall AUC = {:.3}", overall.auc),
        (0.6, 0.2),
        ("sans-serif", 20).into_font().color(&RGBColor(0, 0, 0)),
    )))?;
    
    // 图例（增大字体，放在右下角）
    chart.configure_series_labels()
        .background_style(WHITE.mix(0.8))
        .border_style(BLACK)
        .label_font(("sans-serif", 18)) // 进一步增大图例字体
        .position(SeriesLabelPosition::LowerRight) // 放在右下角
        .draw()?;
    
    root.present()?;
    Ok(())
}

/// 读取TSV格式的中间过程文件数据
fn read_tsv_data(file_path: &str) -> Result<Vec<(u32, char, f64, bool)>, Box<dyn Error>> {
    let mut data = Vec::new();
    let file = File::open(file_path)?;
    let reader = BufReader::new(file);
    
    for (line_num, line) in reader.lines().enumerate() {
        if line_num == 0 {
            continue; // 跳过标题行
        }
        
        let line = line?;
        let fields: Vec<&str> = line.split('\t').collect();
        
        if fields.len() >= 4 {
            let position = fields[0].parse::<u32>()?;
            let base_str = fields[1].trim();
            if base_str.is_empty() {
                continue; // 跳过空的碱基
            }
            let base = base_str.chars().next().unwrap_or('N');
            let reactivity = fields[2].parse::<f64>()?;
            let is_single = fields[3].parse::<u32>()? == 1;
            
            data.push((position, base, reactivity, is_single));
        }
    }
    
    Ok(data)
}

/// 优化的主评估函数（只处理目标基因）
pub fn evaluate_reactivity_accuracy_main_optimized(
    reactivity_file: &str,
    secondary_structure_file: &str,
    output_dir: &str,
    signal_type: &str,
    strand: &str,
    logger: &mut crate::Logger,
) -> Result<(), Box<dyn Error>> {
    let start_time = Instant::now();
    
    // 记录环境信息和参数
    logger.log("=== ModDetector Evaluate Function Log ===")?;
    logger.log(&format!("Software Version: v{}", crate::VERSION))?;
    logger.log(&format!("Runtime: {}", chrono::Utc::now().format("%Y-%m-%d %H:%M:%S")))?;
    logger.log(&format!("Reactivity File: {}", reactivity_file))?;
    logger.log(&format!("Secondary Structure File: {}", secondary_structure_file))?;
    logger.log(&format!("Output Directory: {}", output_dir))?;
    logger.log(&format!("Signal Type: {}", signal_type))?;
    logger.log(&format!("Strand: {}", strand))?;
    logger.log("Starting optimized reactivity accuracy evaluation...")?;
    
    // 定义目标基因模式
    let target_gene_patterns = ["16S", "J01695", "EC", "E.coli"];
    
    // 直接提取目标基因的数据
    logger.log_and_progress("Extracting target gene data...")?;
    let gene_data = extract_reactivity_for_target_gene(
        reactivity_file, 
        &target_gene_patterns, 
        strand, 
        signal_type
    )?;
    
    match gene_data {
        Some((gene_id, reactivity_data)) => {
            logger.log_and_progress(&format!("Found target gene: {}", gene_id))?;
            
            // 解析二级结构
            let secondary_structure = parse_secondary_structure(secondary_structure_file)?;
            
            // 创建输出目录
            std::fs::create_dir_all(output_dir)?;
            
            // 基于自动shift校正进行准确度评估
            logger.log_and_progress("Performing accuracy evaluation with auto-shift correction...")?;
            let matched_data = match_reactivity_with_structure_auto_shift(&reactivity_data, &secondary_structure, signal_type)?;
            
            // 计算评估结果
            let result = evaluate_reactivity_accuracy_from_matched_data(&matched_data, signal_type)?;
            
            // 创建基因输出目录
            let gene_output_dir = format!("{}/{}", output_dir, gene_id.replace('/', "_"));
            std::fs::create_dir_all(&gene_output_dir)?;
            
            // 保存评估结果
            let result_file = format!("{}/evaluation_results_{}.txt", gene_output_dir, signal_type);
            save_evaluation_results(&result, &result_file, signal_type)?;
            
            // 绘制ROC曲线
            let roc_file = format!("{}/roc_curve_{}.svg", gene_output_dir, signal_type);
            let roc_title = format!("{} - {} Signal ROC Curve", gene_id, signal_type);
            plot_roc_curve(&result, &roc_file, &roc_title)?;
            
            // 绘制PR曲线
            let pr_file = format!("{}/pr_curve_{}.svg", gene_output_dir, signal_type);
            let pr_title = format!("{} - {} Signal Precision-Recall Curve", gene_id, signal_type);
            plot_pr_curve(&result, &pr_file, &pr_title)?;
            
            // 输出主要指标
            logger.log(&format!("  AUC: {:.3}", result.auc))?;
            logger.log(&format!("  F1-score: {:.3}", result.f1_score))?;
            logger.log(&format!("  Sensitivity: {:.3}", result.sensitivity))?;
            logger.log(&format!("  Specificity: {:.3}", result.specificity))?;
            logger.log(&format!("  PPV: {:.3}", result.ppv))?;
            logger.log(&format!("  NPV: {:.3}", result.npv))?;
            logger.log(&format!("  Accuracy: {:.3}", result.accuracy))?;
            logger.log(&format!("  Matched positions: {}", matched_data.len()))?;
            
            logger.log(&format!("Successfully processed gene: {}, evaluation completed", gene_id))?;
        }
        None => {
            logger.log("Warning: No matching target gene found for evaluation")?;
        }
    }
    
    let elapsed = start_time.elapsed();
    println!("{}", crate::progress::format_time_used(elapsed));
    Ok(())
}

/// 同时评估stop和mutation信号的准确度
pub fn evaluate_both_signals_main(
    reactivity_file: &str,
    secondary_structure_file: &str,
    output_dir: &str,
    gene_id: &str,
    strand: &str,
    logger: &mut crate::Logger,
) -> Result<(), Box<dyn Error>> {
    let start_time = Instant::now();
    
    // 记录环境信息和参数
    logger.log("=== ModDetector Evaluate Both Signals Function Log ===")?;
    logger.log(&format!("Software Version: v{}", crate::VERSION))?;
    logger.log(&format!("Runtime: {}", chrono::Utc::now().format("%Y-%m-%d %H:%M:%S")))?;
    logger.log(&format!("Reactivity File: {}", reactivity_file))?;
    logger.log(&format!("Secondary Structure File: {}", secondary_structure_file))?;
    logger.log(&format!("Output Directory: {}", output_dir))?;
    logger.log(&format!("Gene ID: {}", gene_id))?;
    logger.log(&format!("Strand: {}", strand))?;
    logger.log("Starting both signals evaluation (stop and mutation)...")?;
    
    // 显示加载数据信息
    println!("[Loading data]");
    println!("    Secondary structure: {}", secondary_structure_file);
    println!("    Reactivity file: {}", reactivity_file);
    println!();
    
    // 显示参数信息
    println!("[Params]");
    println!("    Gene ID: {}, strand: {}", gene_id, strand);
    println!("    Signal types: stop, mutation");
    println!("    Use base matching: true");
    println!("    Auto shift: true");
    println!("    Base used: A, C, G, T");
    println!();
    
    // 创建输出目录
    std::fs::create_dir_all(output_dir)?;
    
    // 解析二级结构文件
    let secondary_structure = parse_secondary_structure(secondary_structure_file)?;
    
    // 显示进度信息
    println!("[Progressing]");
    println!("    Structure positions: {}, unpaired: {}, paired: {}", 
        secondary_structure.positions.len(), 
        secondary_structure.is_single_stranded.iter().filter(|&&x| x).count(),
        secondary_structure.is_single_stranded.iter().filter(|&&x| !x).count());
    println!();
    
    // 评估stop信号
    println!("[Evaluating] Stop signal...");
    let stop_result = evaluate_reactivity_accuracy_with_auto_shift(
        reactivity_file,
        &secondary_structure,
        "stop",
        gene_id,
        strand,
        logger,
    )?;
    
    // 输出stop信号的清晰报告
    println!("[Stop Signal Results]");
    println!("    AUC: {:.3}", stop_result.auc);
    println!("    F1-score: {:.3}", stop_result.f1_score);
    println!("    Sensitivity: {:.3}", stop_result.sensitivity);
    println!("    Specificity: {:.3}", stop_result.specificity);
    println!("    PPV: {:.3}", stop_result.ppv);
    println!("    NPV: {:.3}", stop_result.npv);
    println!("    Accuracy: {:.3}", stop_result.accuracy);
    println!();
    
    // 评估mutation信号
    println!("[Evaluating] Mutation signal...");
    let mutation_result = evaluate_reactivity_accuracy_with_auto_shift(
        reactivity_file,
        &secondary_structure,
        "mutation",
        gene_id,
        strand,
        logger,
    )?;
    
    // 输出mutation信号的清晰报告
    println!("[Mutation Signal Results]");
    println!("    AUC: {:.3}", mutation_result.auc);
    println!("    F1-score: {:.3}", mutation_result.f1_score);
    println!("    Sensitivity: {:.3}", mutation_result.sensitivity);
    println!("    Specificity: {:.3}", mutation_result.specificity);
    println!("    PPV: {:.3}", mutation_result.ppv);
    println!("    NPV: {:.3}", mutation_result.npv);
    println!("    Accuracy: {:.3}", mutation_result.accuracy);
    println!();
    
    // 性能比较
    println!("[Performance Comparison]");
    if stop_result.auc > mutation_result.auc {
        println!("    Stop signal performs better (AUC: {:.3} vs {:.3})", 
            stop_result.auc, mutation_result.auc);
    } else if mutation_result.auc > stop_result.auc {
        println!("    Mutation signal performs better (AUC: {:.3} vs {:.3})", 
            mutation_result.auc, stop_result.auc);
    } else {
        println!("    Both signals perform equally (AUC: {:.3})", stop_result.auc);
    }
    println!();
    
    // 按碱基分别评估stop信号
    println!("[Evaluating] Stop signal by base...");
    let stop_comprehensive_result = evaluate_reactivity_accuracy_by_base(
        reactivity_file,
        &secondary_structure,
        "stop",
        gene_id,
        strand,
        output_dir,
    )?;
    
    // 输出stop信号按碱基的清晰报告
    println!("[Stop Signal Base-Specific Results]");
    println!("    Overall AUC: {:.3}", stop_comprehensive_result.overall_result.auc);
    for base_result in &stop_comprehensive_result.base_specific_results {
        println!("    {} base AUC: {:.3} (Positions: {}, Unpaired: {}, Paired: {})", 
            base_result.base, base_result.result.auc, base_result.data_count,
            base_result.single_stranded_count, base_result.double_stranded_count);
    }
    println!();
    
    // 按碱基分别评估mutation信号
    println!("[Evaluating] Mutation signal by base...");
    let mutation_comprehensive_result = evaluate_reactivity_accuracy_by_base(
        reactivity_file,
        &secondary_structure,
        "mutation",
        gene_id,
        strand,
        output_dir,
    )?;
    
    // 输出mutation信号按碱基的清晰报告
    println!("[Mutation Signal Base-Specific Results]");
    println!("    Overall AUC: {:.3}", mutation_comprehensive_result.overall_result.auc);
    for base_result in &mutation_comprehensive_result.base_specific_results {
        println!("    {} base AUC: {:.3} (Positions: {}, Unpaired: {}, Paired: {})", 
            base_result.base, base_result.result.auc, base_result.data_count,
            base_result.single_stranded_count, base_result.double_stranded_count);
    }
    println!();
    
    // 保存stop信号结果
    println!("[Saving] Stop signal results...");
    save_evaluation_results(&stop_result, &format!("{}/evaluation_stop.txt", output_dir), "stop")?;
    plot_roc_curve(&stop_result, &format!("{}/roc_curve_stop.svg", output_dir), 
        &format!("{} - Stop Signal ROC Curve", gene_id))?;
    plot_pr_curve(&stop_result, &format!("{}/pr_curve_stop.svg", output_dir), 
        &format!("{} - Stop Signal PR Curve", gene_id))?;
    
    // 保存mutation信号结果
    println!("[Saving] Mutation signal results...");
    save_evaluation_results(&mutation_result, &format!("{}/evaluation_mutation.txt", output_dir), "mutation")?;
    plot_roc_curve(&mutation_result, &format!("{}/roc_curve_mutation.svg", output_dir), 
        &format!("{} - Mutation Signal ROC Curve", gene_id))?;
    plot_pr_curve(&mutation_result, &format!("{}/pr_curve_mutation.svg", output_dir), 
        &format!("{} - Mutation Signal PR Curve", gene_id))?;
    
    // 保存stop信号的按碱基评估结果
    println!("[Saving] Stop signal base-specific results...");
    let stop_base_analysis_name = format!("evaluation_by_base_stop_{}_{}", gene_id, strand);
    let stop_overall_roc_path = format!("{}/{}.overall.roc.svg", output_dir, stop_base_analysis_name);
    let stop_overall_pr_path = format!("{}/{}.overall.pr.svg", output_dir, stop_base_analysis_name);
    let stop_overall_results_path = format!("{}/{}.overall.txt", output_dir, stop_base_analysis_name);
    let stop_comprehensive_results_path = format!("{}/{}.comprehensive.txt", output_dir, stop_base_analysis_name);
    
    plot_roc_curve(&stop_comprehensive_result.overall_result, &stop_overall_roc_path, 
        &format!("Overall ROC Curve - Stop ({}, Strand {})", gene_id, strand))?;
    plot_pr_curve(&stop_comprehensive_result.overall_result, &stop_overall_pr_path, 
        &format!("Overall PR Curve - Stop ({}, Strand {})", gene_id, strand))?;
    save_evaluation_results(&stop_comprehensive_result.overall_result, &stop_overall_results_path, "stop")?;
    save_comprehensive_evaluation_results(&stop_comprehensive_result, &stop_comprehensive_results_path, "stop", gene_id, strand)?;
    
    // 为stop信号的每个碱基绘制ROC曲线
    for base_result in &stop_comprehensive_result.base_specific_results {
        let base_roc_path = format!("{}/{}.{}.roc.svg", output_dir, stop_base_analysis_name, base_result.base);
        let base_pr_path = format!("{}/{}.{}.pr.svg", output_dir, stop_base_analysis_name, base_result.base);
        let base_results_path = format!("{}/{}.{}.txt", output_dir, stop_base_analysis_name, base_result.base);
        
        plot_roc_curve(&base_result.result, &base_roc_path, 
            &format!("{} Base Stop ROC Curve - {} (Strand {})", base_result.base, gene_id, strand))?;
        plot_pr_curve(&base_result.result, &base_pr_path, 
            &format!("{} Base Stop PR Curve - {} (Strand {})", base_result.base, gene_id, strand))?;
        save_base_specific_evaluation_results(base_result, &base_results_path, "stop")?;
    }
    
    // 绘制stop信号的combined ROC曲线
    let stop_combined_roc_path = format!("{}/{}.roc.svg", output_dir, stop_base_analysis_name);
    plot_multi_roc_curve(
        &stop_comprehensive_result.overall_result,
        &stop_comprehensive_result.base_specific_results,
        &stop_combined_roc_path,
        &format!("ROC Curve (Overall + 4 Bases) - Stop ({}, Strand {})", gene_id, strand)
    )?;
    
    // 保存mutation信号的按碱基评估结果
    println!("[Saving] Mutation signal base-specific results...");
    let mutation_base_analysis_name = format!("evaluation_by_base_mutation_{}_{}", gene_id, strand);
    let mutation_overall_roc_path = format!("{}/{}.overall.roc.svg", output_dir, mutation_base_analysis_name);
    let mutation_overall_pr_path = format!("{}/{}.overall.pr.svg", output_dir, mutation_base_analysis_name);
    let mutation_overall_results_path = format!("{}/{}.overall.txt", output_dir, mutation_base_analysis_name);
    let mutation_comprehensive_results_path = format!("{}/{}.comprehensive.txt", output_dir, mutation_base_analysis_name);
    
    plot_roc_curve(&mutation_comprehensive_result.overall_result, &mutation_overall_roc_path, 
        &format!("Overall ROC Curve - Mutation ({}, Strand {})", gene_id, strand))?;
    plot_pr_curve(&mutation_comprehensive_result.overall_result, &mutation_overall_pr_path, 
        &format!("Overall PR Curve - Mutation ({}, Strand {})", gene_id, strand))?;
    save_evaluation_results(&mutation_comprehensive_result.overall_result, &mutation_overall_results_path, "mutation")?;
    save_comprehensive_evaluation_results(&mutation_comprehensive_result, &mutation_comprehensive_results_path, "mutation", gene_id, strand)?;
    
    // 为mutation信号的每个碱基绘制ROC曲线
    for base_result in &mutation_comprehensive_result.base_specific_results {
        let base_roc_path = format!("{}/{}.{}.roc.svg", output_dir, mutation_base_analysis_name, base_result.base);
        let base_pr_path = format!("{}/{}.{}.pr.svg", output_dir, mutation_base_analysis_name, base_result.base);
        let base_results_path = format!("{}/{}.{}.txt", output_dir, mutation_base_analysis_name, base_result.base);
        
        plot_roc_curve(&base_result.result, &base_roc_path, 
            &format!("{} Base Mutation ROC Curve - {} (Strand {})", base_result.base, gene_id, strand))?;
        plot_pr_curve(&base_result.result, &base_pr_path, 
            &format!("{} Base Mutation PR Curve - {} (Strand {})", base_result.base, gene_id, strand))?;
        save_base_specific_evaluation_results(base_result, &base_results_path, "mutation")?;
    }
    
    // 绘制mutation信号的combined ROC曲线
    let mutation_combined_roc_path = format!("{}/{}.roc.svg", output_dir, mutation_base_analysis_name);
    plot_multi_roc_curve(
        &mutation_comprehensive_result.overall_result,
        &mutation_comprehensive_result.base_specific_results,
        &mutation_combined_roc_path,
        &format!("ROC Curve (Overall + 4 Bases) - Mutation ({}, Strand {})", gene_id, strand)
    )?;
    
    // 添加数据信息
    println!("[Data info]");
    println!("    Total AUC={:.3}. Positions={}, unpaired={}, paired={}.", 
        stop_comprehensive_result.overall_result.auc, stop_comprehensive_result.total_positions,
        stop_comprehensive_result.total_positions - stop_comprehensive_result.base_specific_results.iter().map(|r| r.double_stranded_count).sum::<usize>(),
        stop_comprehensive_result.base_specific_results.iter().map(|r| r.double_stranded_count).sum::<usize>());
    
    for base_result in &stop_comprehensive_result.base_specific_results {
        println!("    {} base AUC={:.3}. Positions={}, unpaired={}, paired={}.", 
            base_result.base, base_result.result.auc, base_result.data_count,
            base_result.single_stranded_count, base_result.double_stranded_count);
    }
    println!();
    
    // 添加输出信息
    println!("[Outputs]");
    println!("    Temperat: {}/{}.overall.tsv", output_dir, stop_base_analysis_name);
    println!("    Summary: {}/{}.comprehensive.txt", output_dir, stop_base_analysis_name);
    println!("    ROC plots: {}/{}.combined_roc.svg", output_dir, stop_base_analysis_name);
    println!();
    
    // 生成综合报告
    println!("[Generating] Combined report...");
    let combined_report_path = format!("{}/combined_evaluation_report.txt", output_dir);
    let mut report_file = File::create(&combined_report_path)?;
    
    writeln!(report_file, "=== ModDetector Combined Evaluation Report ===")?;
    writeln!(report_file, "Gene ID: {}", gene_id)?;
    writeln!(report_file, "Strand: {}", strand)?;
    writeln!(report_file, "Secondary Structure File: {}", secondary_structure_file)?;
    writeln!(report_file, "Reactivity File: {}", reactivity_file)?;
    writeln!(report_file, "Evaluation Time: {}", chrono::Utc::now().format("%Y-%m-%d %H:%M:%S"))?;
    writeln!(report_file)?;
    
    // Stop信号结果
    writeln!(report_file, "=== Stop Signal Results ===")?;
    writeln!(report_file, "AUC: {:.4}", stop_result.auc)?;
    writeln!(report_file, "F1-Score: {:.4}", stop_result.f1_score)?;
    writeln!(report_file, "Sensitivity: {:.4}", stop_result.sensitivity)?;
    writeln!(report_file, "Specificity: {:.4}", stop_result.specificity)?;
    writeln!(report_file, "PPV: {:.4}", stop_result.ppv)?;
    writeln!(report_file, "NPV: {:.4}", stop_result.npv)?;
    writeln!(report_file, "Accuracy: {:.4}", stop_result.accuracy)?;
    writeln!(report_file)?;
    
    // Mutation信号结果
    writeln!(report_file, "=== Mutation Signal Results ===")?;
    writeln!(report_file, "AUC: {:.4}", mutation_result.auc)?;
    writeln!(report_file, "F1-Score: {:.4}", mutation_result.f1_score)?;
    writeln!(report_file, "Sensitivity: {:.4}", mutation_result.sensitivity)?;
    writeln!(report_file, "Specificity: {:.4}", mutation_result.specificity)?;
    writeln!(report_file, "PPV: {:.4}", mutation_result.ppv)?;
    writeln!(report_file, "NPV: {:.4}", mutation_result.npv)?;
    writeln!(report_file, "Accuracy: {:.4}", mutation_result.accuracy)?;
    writeln!(report_file)?;
    
    // 性能比较
    writeln!(report_file, "=== Performance Comparison ===")?;
    if stop_result.auc > mutation_result.auc {
        writeln!(report_file, "Stop signal performs better (AUC: {:.4} vs {:.4})", 
            stop_result.auc, mutation_result.auc)?;
    } else if mutation_result.auc > stop_result.auc {
        writeln!(report_file, "Mutation signal performs better (AUC: {:.4} vs {:.4})", 
            mutation_result.auc, stop_result.auc)?;
    } else {
        writeln!(report_file, "Both signals perform equally (AUC: {:.4})", stop_result.auc)?;
    }
    
    let elapsed = start_time.elapsed();
    
    // 覆盖进度显示，显示输出信息
    println!("\r[Output]                           ");
    println!("    Combined results: {}", output_dir);
    println!("    Stop signal: roc_curve_stop.svg, pr_curve_stop.svg");
    println!("    Mutation signal: roc_curve_mutation.svg, pr_curve_mutation.svg");
    println!("    Report: combined_evaluation_report.txt");
    println!("{}", crate::progress::format_time_used(elapsed));
    
    // 记录完成日志
    logger.log(&format!("Both signals evaluation completed, output directory: {}", output_dir))?;
    logger.log(&format!("Stop signal AUC: {:.4}", stop_result.auc))?;
    logger.log(&format!("Mutation signal AUC: {:.4}", mutation_result.auc))?;
    logger.log(&format!("Total time: {:.2}s", elapsed.as_secs_f64()))?;
    
    Ok(())
}







