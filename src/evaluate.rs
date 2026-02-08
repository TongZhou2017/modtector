use crate::plot::ReactivityData;
use chrono;
use plotters::prelude::*;
use std::collections::HashMap;
use std::error::Error;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::time::Instant;

/// Parse f64 value from string, supporting "NaN" string
fn parse_f64_allow_nan(s: &str) -> Option<f64> {
    let trimmed = s.trim();
    if trimmed.eq_ignore_ascii_case("nan") {
        Some(f64::NAN)
    } else {
        trimmed.parse::<f64>().ok()
    }
}

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
        else if found_sequence && !found_structure && trimmed.chars().all(|c| "ACGU".contains(c))
        {
            sequence.push_str(trimmed);
        }
        // Check if it's a structure line (contains .() etc.)
        else if found_sequence
            && !found_structure
            && trimmed.chars().all(|c| ".()[]{}<> aA ".contains(c))
        {
            structure.push_str(trimmed);
            found_structure = true;
        }
        // Continue reading structure lines
        else if found_sequence
            && found_structure
            && trimmed.chars().all(|c| ".()[]{}<> aA ".contains(c))
        {
            structure.push_str(trimmed);
        }
    }

    // Filter out spaces, only keep valid structure characters
    structure = structure.chars().filter(|&c| c != ' ').collect();

    if sequence.is_empty() || structure.is_empty() {
        return Err("Could not find sequence or structure in file".into());
    }

    if sequence.len() != structure.len() {
        return Err(format!(
            "Sequence length ({}) and structure length ({}) do not match",
            sequence.len(),
            structure.len()
        )
        .into());
    }

    let mut positions = Vec::new();
    let mut is_single_stranded = Vec::new();

    for (i, (seq_char, struct_char)) in sequence.chars().zip(structure.chars()).enumerate() {
        // Only process standard bases
        if "ACGU".contains(seq_char) {
            positions.push(i as u32 + 1); // 1-based position
                                          // Single-stranded: . Others are considered double-stranded
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

/// Extract base information from reactivity CSV file
pub fn extract_reactivity_with_bases(
    csv_file: &str,
    gene_id: &str,
    strand: &str,
    signal_type: &str,
) -> Result<Vec<(u32, char, f64)>, Box<dyn Error>> {
    let mut reactivity_data: Vec<(u32, char, f64)> = Vec::new();
    let file = File::open(csv_file)?;
    let reader = BufReader::new(file);

    // First calculate total line count for progress display
    let total_lines = std::fs::read_to_string(csv_file)?.lines().count();

    for (line_num, line) in reader.lines().enumerate() {
        if line_num == 0 {
            continue; // Skip header row
        }

        // Update progress display
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

            // Determine which column to read reactivity value based on signal_type
            let reactivity_value = if signal_type == "mutation" {
                // mutation signal uses column 6 (reactivity_mutation)
                parse_f64_allow_nan(fields[5]).ok_or_else(|| format!("Invalid reactivity_mutation value: {}", fields[5]))?
            } else {
                // stop signal uses column 5 (reactivity_stop)
                parse_f64_allow_nan(fields[4]).ok_or_else(|| format!("Invalid reactivity_stop value: {}", fields[4]))?
            };

            reactivity_data.push((position, base, reactivity_value));
        }
    }

    // Sort by position
    reactivity_data.sort_by_key(|(pos, _, _)| *pos);

    Ok(reactivity_data)
}

/// Extract base information for specified gene from reactivity CSV file (optimized version)
pub fn extract_reactivity_for_target_gene(
    csv_file: &str,
    target_gene_patterns: &[&str],
    strand: &str,
    signal_type: &str,
) -> Result<Option<(String, Vec<(u32, char, f64)>)>, Box<dyn Error>> {
    let file = File::open(csv_file)?;
    let reader = BufReader::new(file);
    let mut reactivity_data: Vec<(u32, char, f64)> = Vec::new();
    let mut found_gene_id: Option<String> = None;

    for (line_num, line) in reader.lines().enumerate() {
        if line_num == 0 {
            continue; // Skip header row
        }

        let line = line?;
        let fields: Vec<&str> = line.split(',').collect();

        if fields.len() >= 6 {
            let chr_id = fields[0].trim_matches('"');
            let strand_field = fields[1].trim();

            // Check if matches target gene pattern
            let is_target_gene = target_gene_patterns
                .iter()
                .any(|pattern| chr_id.contains(pattern));

            if is_target_gene && strand_field == strand {
                // Found target gene, record gene ID
                if found_gene_id.is_none() {
                    found_gene_id = Some(chr_id.to_string());
                }

                let position = fields[2].parse::<u32>()?;
                let base = fields[3].chars().next().unwrap_or('N');

                // Determine which column to read reactivity value based on signal_type
                let reactivity_value = if signal_type == "mutation" {
                    // mutation signal uses column 6 (reactivity_mutation)
                    fields[5].parse::<f64>()?
                } else {
                    // stop signal uses column 5 (reactivity_stop)
                    fields[4].parse::<f64>()?
                };

                reactivity_data.push((position, base, reactivity_value));
            }
        }
    }

    // Sort by position
    reactivity_data.sort_by_key(|(pos, _, _)| *pos);

    if let Some(gene_id) = found_gene_id {
        println!(
            "Found target gene: {} (strand {}): {} sites",
            gene_id,
            strand,
            reactivity_data.len()
        );
        Ok(Some((gene_id, reactivity_data)))
    } else {
        println!("No matching target gene found");
        Ok(None)
    }
}

/// Match reactivity data with secondary structure data based on sequence bases (supports automatic shift correction)
pub fn match_reactivity_with_structure_auto_shift(
    reactivity_data: &[(u32, char, f64)],
    secondary_structure: &SecondaryStructure,
    _signal_type: &str,
) -> Result<Vec<(f64, bool)>, Box<dyn Error>> {
    // Create base mapping for secondary structure sequence
    let mut structure_bases: Vec<(u32, char, bool)> = Vec::new();
    let mut pos_idx = 0;

    for (_i, (seq_char, struct_char)) in secondary_structure
        .sequence
        .chars()
        .zip(secondary_structure.structure.chars())
        .enumerate()
    {
        if "ACGU".contains(seq_char) {
            let position = pos_idx + 1;
            let is_single = struct_char == '.';
            structure_bases.push((position, seq_char, is_single));
            pos_idx += 1;
        }
    }

    // Calculate max_shift based on length difference
    let length_diff = reactivity_data.len() as i32 - structure_bases.len() as i32;
    let max_shift = (length_diff.abs() + 5).min(20); // Maximum 20, minimum 5

    // Find best shift
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
                    // Handle T and U equivalence
                    let base_match = *react_base == struct_base
                        || (*react_base == 'T' && struct_base == 'U')
                        || (*react_base == 'U' && struct_base == 'T');
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

    // Use best shift for matching, only keep data with corresponding secondary structure sites
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
                // Handle T and U equivalence
                let base_match = *react_base == struct_base
                    || (*react_base == 'T' && struct_base == 'U')
                    || (*react_base == 'U' && struct_base == 'T');
                matched_data.push((*reactivity_value, is_single));
                if base_match {
                    match_count += 1;
                } else {
                    mismatch_count += 1;
                }
            } else {
                // Sites beyond secondary structure range are skipped
                skipped_count += 1;
            }
        } else {
            // Sites beyond secondary structure range are skipped
            skipped_count += 1;
        }
    }

    if matched_data.is_empty() {
        return Err("No matched data found".into());
    }

    Ok(matched_data)
}

/// Calculate accuracy evaluation metrics
pub fn evaluate_reactivity_accuracy(
    reactivity_data: &[ReactivityData],
    secondary_structure: &SecondaryStructure,
    signal_type: &str, // "stop" or "mutation"
) -> Result<EvaluationResult, Box<dyn Error>> {
    // Match reactivity data with secondary structure data
    let mut matched_data: Vec<(f64, bool)> = Vec::new();

    for reactivity in reactivity_data {
        if let Some(pos_idx) = secondary_structure
            .positions
            .iter()
            .position(|&p| p == reactivity.position)
        {
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

    // Sort data, NaN values will be placed last (but will be filtered out during statistics)
    matched_data.sort_by(|a, b| {
        match (a.0.is_nan(), b.0.is_nan()) {
            (true, true) => std::cmp::Ordering::Equal,
            (true, false) => std::cmp::Ordering::Greater, // NaN placed last
            (false, true) => std::cmp::Ordering::Less,
            (false, false) => a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal),
        }
    });

    // Calculate ROC curve data
    let mut thresholds = Vec::new();
    let mut tpr = Vec::new();
    let mut fpr = Vec::new();
    let mut precision = Vec::new();
    let mut recall = Vec::new();

    // Calculate total positive and negative sample counts, filter out NaN values
    let total_positive: usize = matched_data
        .iter()
        .filter(|(value, is_single)| !value.is_nan() && *is_single)
        .count();
    let total_negative: usize = matched_data
        .iter()
        .filter(|(value, is_single)| !value.is_nan() && !*is_single)
        .count();

    if total_positive == 0 || total_negative == 0 {
        return Err("No positive or negative samples found".into());
    }

    // Generate thresholds, filter out NaN values
    let valid_data: Vec<f64> = matched_data.iter()
        .map(|(val, _)| *val)
        .filter(|v| !v.is_nan())
        .collect();
    
    if valid_data.is_empty() {
        return Err("No valid reactivity data (all values are NaN)".into());
    }
    
    let min_val = valid_data.iter().fold(valid_data[0], |a, &b| a.min(b));
    let max_val = valid_data.iter().fold(valid_data[0], |a, &b| a.max(b));
    let step = (max_val - min_val) / 100.0;

    for i in 0..=100 {
        let threshold = min_val + i as f64 * step;
        thresholds.push(threshold);

        let mut tp = 0;
        let mut fp = 0;
        let mut tn = 0;
        let mut fn_count = 0;

        for (value, is_single) in &matched_data {
            // Skip NaN values, do not participate in statistical calculation
            if value.is_nan() {
                continue;
            }
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

        let tpr_val = if tp + fn_count > 0 {
            tp as f64 / (tp + fn_count) as f64
        } else {
            0.0
        };
        let fpr_val = if fp + tn > 0 {
            fp as f64 / (fp + tn) as f64
        } else {
            0.0
        };
        let precision_val = if tp + fp > 0 {
            tp as f64 / (tp + fp) as f64
        } else {
            0.0
        };
        let recall_val = tpr_val;

        tpr.push(tpr_val);
        fpr.push(fpr_val);
        precision.push(precision_val);
        recall.push(recall_val);
    }

    // Calculate AUC (using trapezoidal rule)
    // Ensure ROC curve is correct: FPR should go from 0 to 1, TPR should go from 0 to 1
    let mut auc = 0.0;

    // Add start point (0,0) and end point (1,1)
    let mut roc_points = vec![(0.0, 0.0)];
    for i in 0..fpr.len() {
        roc_points.push((fpr[i], tpr[i]));
    }
    roc_points.push((1.0, 1.0));

    // Calculate AUC
    for i in 1..roc_points.len() {
        let (fpr_prev, tpr_prev) = roc_points[i - 1];
        let (fpr_curr, tpr_curr) = roc_points[i];
        auc += (fpr_curr - fpr_prev) * (tpr_curr + tpr_prev) / 2.0;
    }

    // Calculate other metrics, filter out NaN values
    let tp = matched_data
        .iter()
        .filter(|(value, is_single)| !value.is_nan() && *value > 0.0 && *is_single)
        .count();
    let fp = matched_data
        .iter()
        .filter(|(value, is_single)| !value.is_nan() && *value > 0.0 && !*is_single)
        .count();
    let tn = matched_data
        .iter()
        .filter(|(value, is_single)| !value.is_nan() && *value <= 0.0 && !*is_single)
        .count();
    let fn_count = matched_data
        .iter()
        .filter(|(value, is_single)| !value.is_nan() && *value <= 0.0 && *is_single)
        .count();

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
    let accuracy = (tp + tn) as f64 / (tp + tn + fp + fn_count) as f64;
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
        tpr,
        fpr,
        precision,
        recall,
    })
}

/// Accuracy evaluation based on base matching
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

    // Extract reactivity data and base information
    let reactivity_data =
        extract_reactivity_with_bases(reactivity_file, gene_id, strand, signal_type)?;

    // Match data based on bases
    let matched_data = match_reactivity_with_structure_auto_shift(
        &reactivity_data,
        secondary_structure,
        signal_type,
    )?;

    // Sort data, NaN values will be placed last (but will be filtered out during statistics)
    let mut sorted_data = matched_data;
    sorted_data.sort_by(|a, b| {
        match (a.0.is_nan(), b.0.is_nan()) {
            (true, true) => std::cmp::Ordering::Equal,
            (true, false) => std::cmp::Ordering::Greater, // NaN placed last
            (false, true) => std::cmp::Ordering::Less,
            (false, false) => a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal),
        }
    });

    // Calculate ROC curve data
    let mut thresholds = Vec::new();
    let mut tpr = Vec::new();
    let mut fpr = Vec::new();
    let mut precision = Vec::new();
    let mut recall = Vec::new();

    // Calculate total positive and negative sample counts, filter out NaN values
    let total_positive: usize = sorted_data
        .iter()
        .filter(|(value, is_single)| !value.is_nan() && *is_single)
        .count();
    let total_negative: usize = sorted_data
        .iter()
        .filter(|(value, is_single)| !value.is_nan() && !*is_single)
        .count();

    if total_positive == 0 || total_negative == 0 {
        return Err("No positive or negative samples found".into());
    }

    logger.log(&format!(
        "Sample Statistics: Single-stranded sites={}, Double-stranded sites={}",
        total_positive, total_negative
    ))?;

    // Generate ROC curve data (based on actual data points)
    let mut roc_points = Vec::new();

    // Add start point (0,0)
    roc_points.push((0.0, 0.0));

    // Calculate ROC points for each unique data value, filter out NaN values
    let mut unique_values: Vec<f64> = sorted_data.iter()
        .map(|(val, _)| *val)
        .filter(|v| !v.is_nan())
        .collect();
    unique_values.sort_by(|a, b| {
        match (a.is_nan(), b.is_nan()) {
            (true, true) => std::cmp::Ordering::Equal,
            (true, false) => std::cmp::Ordering::Greater, // NaN placed last
            (false, true) => std::cmp::Ordering::Less,
            (false, false) => a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal),
        }
    });
    unique_values.dedup(); // Remove duplicate values

    // Generate threshold list (for output) - must be set before calculating precision/recall
    thresholds = unique_values.clone();

    // Calculate precision, recall, and FPR for each threshold (in threshold order)
    for threshold in &thresholds {
        let mut tp = 0;
        let mut fp = 0;
        let mut tn = 0;
        let mut fn_count = 0;

        for (value, is_single) in &sorted_data {
            // Skip NaN values, do not participate in statistical calculation
            if value.is_nan() {
                continue;
            }
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

        let tpr_val = if tp + fn_count > 0 {
            tp as f64 / (tp + fn_count) as f64
        } else {
            0.0
        };
        let fpr_val = if fp + tn > 0 {
            fp as f64 / (fp + tn) as f64
        } else {
            0.0
        };
        let precision_val = if tp + fp > 0 {
            tp as f64 / (tp + fp) as f64
        } else {
            0.0
        };
        let recall_val = if tp + fn_count > 0 {
            tp as f64 / (tp + fn_count) as f64
        } else {
            0.0
        };

        // Store FPR, TPR, precision, and recall in the same order as thresholds
        fpr.push(fpr_val);
        tpr.push(tpr_val);
        precision.push(precision_val);
        recall.push(recall_val);
        
        // Also add to roc_points for AUC calculation
        roc_points.push((fpr_val, tpr_val));
    }

    // Add end point (1,1)
    roc_points.push((1.0, 1.0));

    // Ensure ROC points are monotonic (FPR increasing) for AUC calculation
    roc_points.sort_by(|a, b| {
        // FPR and TPR should not be NaN, but handle it for safety
        match (a.0.is_nan(), b.0.is_nan()) {
            (true, true) => std::cmp::Ordering::Equal,
            (true, false) => std::cmp::Ordering::Greater,
            (false, true) => std::cmp::Ordering::Less,
            (false, false) => a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal),
        }
    });

    // Calculate AUC (using trapezoidal rule)
    let mut auc = 0.0;
    for i in 1..roc_points.len() {
        let (fpr_prev, tpr_prev) = roc_points[i - 1];
        let (fpr_curr, tpr_curr) = roc_points[i];
        auc += (fpr_curr - fpr_prev) * (tpr_curr + tpr_prev) / 2.0;
    }

    // Calculate other metrics - use correct confusion matrix calculation
    // Select best threshold using Youden's J statistic (Sensitivity + Specificity - 1)
    // This avoids selecting thresholds that lead to extreme results (Sensitivity=1.0, Specificity=0.0)
    let mut best_youden_j = -1.0;
    let mut best_threshold_idx: Option<usize> = None;

    for i in 0..precision.len() {
        // Calculate sensitivity and specificity for this threshold
        let sensitivity = recall[i]; // recall is TPR (sensitivity)
        let specificity = if fpr[i] < 1.0 {
            1.0 - fpr[i] // specificity = 1 - FPR
        } else {
            0.0
        };
        
        // Calculate Youden's J statistic
        let youden_j = sensitivity + specificity - 1.0;
        
        // Skip thresholds that lead to extreme results
        if (sensitivity == 1.0 && specificity == 0.0) || (sensitivity == 0.0 && specificity == 1.0) {
            continue;
        }
        
        if youden_j > best_youden_j {
            best_youden_j = youden_j;
            best_threshold_idx = Some(i);
        }
    }

    // Calculate confusion matrix using best threshold
    // If no valid threshold found (all were extreme), fall back to F1-score based selection
    let threshold = if let Some(idx) = best_threshold_idx {
        if idx < thresholds.len() {
            thresholds[idx]
        } else if !thresholds.is_empty() {
            thresholds[thresholds.len() - 1]
        } else {
            0.0
        }
    } else {
        // Fallback: if all thresholds were extreme, use F1-score to find best threshold
        let mut best_f1 = 0.0;
        let mut fallback_idx = 0;
        for i in 0..precision.len() {
            let f1 = if precision[i] + recall[i] > 0.0 {
                2.0 * precision[i] * recall[i] / (precision[i] + recall[i])
            } else {
                0.0
            };
            if f1 > best_f1 {
                best_f1 = f1;
                fallback_idx = i;
            }
        }
        if fallback_idx < thresholds.len() {
            thresholds[fallback_idx]
        } else if !thresholds.is_empty() {
            thresholds[thresholds.len() - 1]
        } else {
            0.0
        }
    };

    let mut tp = 0;
    let mut fp = 0;
    let mut tn = 0;
    let mut fn_count = 0;

    for (value, is_single) in &sorted_data {
        // Skip NaN values, do not participate in statistical calculation
        if value.is_nan() {
            continue;
        }
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

    // Calculate correct metrics
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

    logger.log(&format!(
        "Evaluation completed: AUC={:.4}, F1={:.4}, Sensitivity={:.4}, Specificity={:.4}",
        auc, f1_score, sensitivity, specificity
    ))?;

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

/// Plot ROC curve
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
        .x_label_area_size(50) // Increase label area
        .y_label_area_size(50) // Increase label area
        .build_cartesian_2d(0.0..1.0, 0.0..1.0)?;

    // Completely remove grid lines, only keep axes
    chart
        .configure_mesh()
        .x_desc("False Positive Rate")
        .y_desc("True Positive Rate")
        .x_labels(5) // Reduce x-axis tick count
        .y_labels(5) // Reduce y-axis tick count
        .disable_mesh() // Completely disable grid lines
        .draw()?;

    // Plot ROC curve (including start point (0,0) and end point (1,1))
    let mut roc_points = vec![(0.0, 0.0)];
    for i in 0..result.fpr.len() {
        roc_points.push((result.fpr[i], result.tpr[i]));
    }
    roc_points.push((1.0, 1.0));

    chart.draw_series(LineSeries::new(
        roc_points,
        RGBColor(255, 0, 0).stroke_width(2),
    ))?;

    // Plot diagonal line
    chart.draw_series(LineSeries::new(
        vec![(0.0, 0.0), (1.0, 1.0)],
        RGBColor(0, 0, 0).stroke_width(1),
    ))?;

    // Add AUC annotation (increase font size)
    chart.draw_series(std::iter::once(Text::new(
        format!("AUC = {:.3}", result.auc),
        (0.6, 0.2),
        ("sans-serif", 24).into_font().color(&RGBColor(255, 0, 0)),
    )))?;

    root.present()?;
    Ok(())
}

/// Plot Precision-Recall curve
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

    // Plot PR curve
    chart.draw_series(LineSeries::new(
        result
            .recall
            .iter()
            .zip(result.precision.iter())
            .map(|(&x, &y)| (x, y)),
        RGBColor(0, 0, 255).stroke_width(2),
    ))?;

    // Add F1-score annotation
    chart.draw_series(std::iter::once(Text::new(
        format!("F1-score = {:.3}", result.f1_score),
        (0.6, 0.2),
        ("sans-serif", 20).into_font().color(&RGBColor(0, 0, 255)),
    )))?;

    root.present()?;
    Ok(())
}

/// Output evaluation results to file
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

    // Use shortest array length to avoid out of bounds
    let min_len = result
        .thresholds
        .len()
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

/// Main evaluation function
pub fn evaluate_reactivity_accuracy_main(
    reactivity_file: &str,
    secondary_structure_file: &str,
    output_dir: &str,
    signal_type: &str,
    logger: &mut crate::Logger,
) -> Result<(), Box<dyn Error>> {
    let start_time = Instant::now();
    // Load reactivity data
    let reactivity_data = crate::plot::load_reactivity_data(reactivity_file)?;

    // Parse secondary structure
    let secondary_structure = parse_secondary_structure(secondary_structure_file)?;

    // Create output directory
    std::fs::create_dir_all(output_dir)?;

    // Find matching genes, process immediately when found and return
    for (gene_id, gene_reactivity) in reactivity_data {
        logger.log_and_progress(&format!("Checking gene: {}", gene_id))?;

        // Check if gene ID matches secondary structure (assume 16S corresponds to specific gene ID)
        // This can be adjusted according to actual needs
        if gene_id.contains("16S")
            || gene_id.contains("J01695")
            || gene_id.contains("EC")
            || gene_id.contains("E.coli")
        {
            logger.log_and_progress(&format!("Found matching gene: {}", gene_id))?;

            // Evaluate accuracy
            let result =
                evaluate_reactivity_accuracy(&gene_reactivity, &secondary_structure, signal_type)?;

            // Create gene output directory
            let gene_output_dir = format!("{}/{}", output_dir, gene_id.replace('/', "_"));
            std::fs::create_dir_all(&gene_output_dir)?;

            // Save evaluation results
            let result_file = format!("{}/evaluation_results_{}.txt", gene_output_dir, signal_type);
            save_evaluation_results(&result, &result_file, signal_type)?;

            // Plot ROC curve
            let roc_file = format!("{}/roc_curve_{}.svg", gene_output_dir, signal_type);
            let roc_title = format!("{} - {} Signal ROC Curve", gene_id, signal_type);
            plot_roc_curve(&result, &roc_file, &roc_title)?;

            // Plot PR curve
            let pr_file = format!("{}/pr_curve_{}.svg", gene_output_dir, signal_type);
            let pr_title = format!(
                "{} - {} Signal Precision-Recall Curve",
                gene_id, signal_type
            );
            plot_pr_curve(&result, &pr_file, &pr_title)?;

            // Output main metrics
            logger.log(&format!("  AUC: {:.3}", result.auc))?;
            logger.log(&format!("  F1-score: {:.3}", result.f1_score))?;
            logger.log(&format!("  Sensitivity: {:.3}", result.sensitivity))?;
            logger.log(&format!("  Specificity: {:.3}", result.specificity))?;
            logger.log(&format!("  PPV: {:.3}", result.ppv))?;
            logger.log(&format!("  NPV: {:.3}", result.npv))?;
            logger.log(&format!("  Accuracy: {:.3}", result.accuracy))?;
            logger.log(&format!("  Matched positions: {}", gene_reactivity.len()))?;

            // Return immediately after finding matching gene, no longer continue traversal
            logger.log(&format!("Successfully processed gene: {}, evaluation completed", gene_id))?;
            let elapsed = start_time.elapsed();
            logger.log(&format!("[evaluate] Total time: {:.9}s", elapsed.as_secs_f64()))?;
            return Ok(());
        } else {
            logger.log(&format!("Skipping gene: {} (does not match 16S structure)", gene_id))?;
        }
    }

    // If no matching gene found
    logger.log("Warning: No matching gene found for evaluation")?;
    let elapsed = start_time.elapsed();
    println!("{}", crate::progress::format_time_used(elapsed));
    Ok(())
}

/// Main function for accuracy evaluation based on base matching
pub fn evaluate_reactivity_accuracy_with_base_matching_main(
    reactivity_file: &str,
    secondary_structure_file: &str,
    output_dir: &str,
    signal_type: &str,
    gene_id: &str,
    strand: &str,
    reactive_bases: &str,
    logger: &mut crate::Logger,
) -> Result<(), Box<dyn Error>> {
    let start_time = Instant::now();

    // Record environment information and parameters
    logger.log("=== ModDetector Evaluate Function Log ===")?;
    logger.log(&format!("Software Version: v{}", crate::VERSION))?;
    logger.log(&format!(
        "Runtime: {}",
        chrono::Utc::now().format("%Y-%m-%d %H:%M:%S")
    ))?;
    logger.log(&format!("Reactivity File: {}", reactivity_file))?;
    logger.log(&format!(
        "Secondary Structure File: {}",
        secondary_structure_file
    ))?;
    logger.log(&format!("Output Directory: {}", output_dir))?;
    logger.log(&format!("Signal Type: {}", signal_type))?;
    logger.log(&format!("Gene ID: {}", gene_id))?;
    logger.log(&format!("Strand: {}", strand))?;
    logger.log(&format!("Use Base Matching: true"))?;
    logger.log(&format!("Auto Shift: true"))?;
    logger.log(&format!("Reactive bases (for Overall/ROC): {}", reactive_bases))?;
    logger.log("Starting reactivity accuracy evaluation with base matching...")?;

    // Display data loading information
    println!("[Loading data]");
    println!("    Secondary structure: {}", secondary_structure_file);
    println!("    Reactivity file: {}", reactivity_file);
    println!();

    // Display parameter information
    println!("[Params]");
    println!("    Gene ID: {}, strand: {}", gene_id, strand);
    println!("    Signal type: {}", signal_type);
    println!("    Use base matching: true");
    println!("    Auto shift: true");
    println!("    Reactive bases (Overall/ROC): {}", reactive_bases);
    println!();

    // Create output directory
    std::fs::create_dir_all(output_dir)?;

    // Parse secondary structure file
    let secondary_structure = parse_secondary_structure(secondary_structure_file)?;

    // Extract reactivity data first to get actual count
    let reactivity_data = extract_reactivity_with_bases(
        reactivity_file,
        gene_id,
        strand,
        signal_type,
    )?;

    // Display progress information
    println!("[Progressing]");
    println!(
        "    Structure positions: {}, unpaired: {}, paired: {}",
        secondary_structure.positions.len(),
        secondary_structure
            .is_single_stranded
            .iter()
            .filter(|&&x| x)
            .count(),
        secondary_structure
            .is_single_stranded
            .iter()
            .filter(|&&x| !x)
            .count()
    );
    println!("    Reactivity positions: {}", reactivity_data.len());
    let length_diff = reactivity_data.len() as i32 - secondary_structure.positions.len() as i32;
    println!("    Length difference: {}, Max shift: {}", length_diff, (length_diff.abs() + 5).min(20));
    println!();

    // Perform accuracy evaluation based on automatic shift correction

    // 1. Execute original base matching evaluation
    let original_result = evaluate_reactivity_accuracy_with_base_matching(
        reactivity_file,
        &secondary_structure,
        signal_type,
        gene_id,
        strand,
        logger,
    )?;

    // Generate output file names for original evaluation
    let base_name = format!("evaluation_base_matching_{}_{}", signal_type, strand);
    let roc_plot_path = format!("{}/{}.roc.svg", output_dir, base_name);
    let pr_plot_path = format!("{}/{}.pr.svg", output_dir, base_name);
    let results_path = format!("{}/{}.txt", output_dir, base_name);

    // Plot ROC curve for original evaluation
    plot_roc_curve(
        &original_result,
        &roc_plot_path,
        &format!(
            "ROC Curve - {} (Base Matching, Strand {})",
            signal_type, strand
        ),
    )?;

    // Plot PR curve for original evaluation
    plot_pr_curve(
        &original_result,
        &pr_plot_path,
        &format!(
            "PR Curve - {} (Base Matching, Strand {})",
            signal_type, strand
        ),
    )?;

    // Save original evaluation results
    save_evaluation_results(&original_result, &results_path, signal_type)?;

    // 2. Execute evaluation by base separately
    let comprehensive_result = evaluate_reactivity_accuracy_by_base(
        reactivity_file,
        &secondary_structure,
        signal_type,
        gene_id,
        strand,
        output_dir,
        reactive_bases,
    )?;

    // Generate output file names for base-specific evaluation
    let base_analysis_name = format!("evaluation_by_base_{}_{}_{}", signal_type, gene_id, strand);
    let overall_roc_plot_path = format!("{}/{}.overall.roc.svg", output_dir, base_analysis_name);
    let overall_pr_plot_path = format!("{}/{}.overall.pr.svg", output_dir, base_analysis_name);
    let overall_results_path = format!("{}/{}.overall.txt", output_dir, base_analysis_name);
    let comprehensive_results_path =
        format!("{}/{}.comprehensive.txt", output_dir, base_analysis_name);

    // Plot overall ROC curve
    plot_roc_curve(
        &comprehensive_result.overall_result,
        &overall_roc_plot_path,
        &format!(
            "Overall ROC Curve - {} ({}, Strand {})",
            signal_type, gene_id, strand
        ),
    )?;

    // Plot overall PR curve
    plot_pr_curve(
        &comprehensive_result.overall_result,
        &overall_pr_plot_path,
        &format!(
            "Overall PR Curve - {} ({}, Strand {})",
            signal_type, gene_id, strand
        ),
    )?;

    // Save overall evaluation results
    save_evaluation_results(
        &comprehensive_result.overall_result,
        &overall_results_path,
        signal_type,
    )?;

    // Save comprehensive evaluation results (including base-specific evaluations)
    save_comprehensive_evaluation_results(
        &comprehensive_result,
        &comprehensive_results_path,
        signal_type,
        gene_id,
        strand,
    )?;

    // Plot ROC curve for each base
    for base_result in &comprehensive_result.base_specific_results {
        let base_roc_plot_path = format!(
            "{}/{}.{}.roc.svg",
            output_dir, base_analysis_name, base_result.base
        );
        let base_pr_plot_path = format!(
            "{}/{}.{}.pr.svg",
            output_dir, base_analysis_name, base_result.base
        );
        let base_results_path = format!(
            "{}/{}.{}.txt",
            output_dir, base_analysis_name, base_result.base
        );

        // Plot base-specific ROC curve
        plot_roc_curve(
            &base_result.result,
            &base_roc_plot_path,
            &format!(
                "{} Base {} ROC Curve - {} (Strand {})",
                base_result.base, signal_type, gene_id, strand
            ),
        )?;

        // Plot base-specific PR curve
        plot_pr_curve(
            &base_result.result,
            &base_pr_plot_path,
            &format!(
                "{} Base {} PR Curve - {} (Strand {})",
                base_result.base, signal_type, gene_id, strand
            ),
        )?;

        // Save base-specific evaluation results
        save_base_specific_evaluation_results(base_result, &base_results_path, signal_type)?;
    }

    // Plot overall+4 bases ROC curve in one figure
    let combined_roc_path = format!("{}/{}.roc.svg", output_dir, base_analysis_name);
    plot_multi_roc_curve(
        &comprehensive_result.overall_result,
        &comprehensive_result.base_specific_results,
        &combined_roc_path,
        &format!(
            "ROC Curve (Overall + 4 Bases) - {} ({}, Strand {})",
            signal_type, gene_id, strand
        ),
    )?;

    // Add data information
    println!("[Data info]");
    println!(
        "    Total AUC={:.3}. Positions={}, unpaired={}, paired={}.",
        comprehensive_result.overall_result.auc,
        comprehensive_result.total_positions,
        comprehensive_result.total_positions
            - comprehensive_result
                .base_specific_results
                .iter()
                .map(|r| r.double_stranded_count)
                .sum::<usize>(),
        comprehensive_result
            .base_specific_results
            .iter()
            .map(|r| r.double_stranded_count)
            .sum::<usize>()
    );

    for base_result in &comprehensive_result.base_specific_results {
        println!(
            "    {} base AUC={:.3}. Positions={}, unpaired={}, paired={}.",
            base_result.base,
            base_result.result.auc,
            base_result.data_count,
            base_result.single_stranded_count,
            base_result.double_stranded_count
        );
    }
    println!();

    // Add output information
    println!("[Outputs]");
    println!(
        "    Temperat: {}/{}.overall.tsv",
        output_dir, base_analysis_name
    );
    println!(
        "    Summary: {}/{}.comprehensive.txt",
        output_dir, base_analysis_name
    );
    println!(
        "    ROC plots: {}/{}.combined_roc.svg",
        output_dir, base_analysis_name
    );
    println!();

    logger.finish_progress()?;
    let elapsed = start_time.elapsed();
    println!("{}", crate::progress::format_time_used(elapsed));

    Ok(())
}

/// Main function for accuracy evaluation based on automatic shift correction
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

    // Record environment information and parameters
    logger.log("=== ModDetector Evaluate Function Log ===")?;
    logger.log(&format!("Software Version: v{}", crate::VERSION))?;
    logger.log(&format!(
        "Runtime: {}",
        chrono::Utc::now().format("%Y-%m-%d %H:%M:%S")
    ))?;
    logger.log(&format!("Reactivity File: {}", reactivity_file))?;
    logger.log(&format!(
        "Secondary Structure File: {}",
        secondary_structure_file
    ))?;
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
    plot_roc_curve(
        &result,
        &roc_plot_path,
        &format!(
            "ROC Curve - {} (Auto Shift, Strand {})",
            signal_type, strand
        ),
    )?;
    logger.log_and_progress("Plotting PR curve...")?;
    plot_pr_curve(
        &result,
        &pr_plot_path,
        &format!("PR Curve - {} (Auto Shift, Strand {})", signal_type, strand),
    )?;
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

/// Accuracy evaluation based on automatic shift correction
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

    // Extract reactivity data and base information
    let reactivity_data =
        extract_reactivity_with_bases(reactivity_file, gene_id, strand, signal_type)?;

    // Match data based on automatic shift correction
    let matched_data = match_reactivity_with_structure_auto_shift(
        &reactivity_data,
        secondary_structure,
        signal_type,
    )?;

    // Sort data, NaN values will be placed last (but will be filtered out during statistics)
    let mut sorted_data = matched_data;
    sorted_data.sort_by(|a, b| {
        match (a.0.is_nan(), b.0.is_nan()) {
            (true, true) => std::cmp::Ordering::Equal,
            (true, false) => std::cmp::Ordering::Greater, // NaN placed last
            (false, true) => std::cmp::Ordering::Less,
            (false, false) => a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal),
        }
    });

    // Calculate ROC curve data
    let mut thresholds = Vec::new();
    let mut tpr = Vec::new();
    let mut fpr = Vec::new();
    let mut precision = Vec::new();
    let mut recall = Vec::new();

    // Calculate total positive and negative sample counts, filter out NaN values
    let total_positive: usize = sorted_data
        .iter()
        .filter(|(value, is_single)| !value.is_nan() && *is_single)
        .count();
    let total_negative: usize = sorted_data
        .iter()
        .filter(|(value, is_single)| !value.is_nan() && !*is_single)
        .count();

    if total_positive == 0 || total_negative == 0 {
        return Err("No positive or negative samples found".into());
    }

    logger.log(&format!(
        "Sample Statistics: Single-stranded sites={}, Double-stranded sites={}",
        total_positive, total_negative
    ))?;

    // Generate ROC curve data (based on actual data points)
    let mut roc_points = Vec::new();

    // Add start point (0,0)
    roc_points.push((0.0, 0.0));

    // Calculate ROC points for each unique data value, filter out NaN values
    let mut unique_values: Vec<f64> = sorted_data.iter()
        .map(|(val, _)| *val)
        .filter(|v| !v.is_nan())
        .collect();
    unique_values.sort_by(|a, b| {
        match (a.is_nan(), b.is_nan()) {
            (true, true) => std::cmp::Ordering::Equal,
            (true, false) => std::cmp::Ordering::Greater, // NaN placed last
            (false, true) => std::cmp::Ordering::Less,
            (false, false) => a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal),
        }
    });
    unique_values.dedup(); // Remove duplicate values

    for threshold in unique_values {
        let mut tp = 0;
        let mut fp = 0;
        let mut tn = 0;
        let mut fn_count = 0;

        for (value, is_single) in &sorted_data {
            // Skip NaN values, do not participate in statistical calculation
            if value.is_nan() {
                continue;
            }
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

        let tpr_val = if total_positive > 0 {
            tp as f64 / total_positive as f64
        } else {
            0.0
        };
        let fpr_val = if total_negative > 0 {
            fp as f64 / total_negative as f64
        } else {
            0.0
        };
        let precision_val = if (tp + fp) > 0 {
            tp as f64 / (tp + fp) as f64
        } else {
            0.0
        };
        let recall_val = if total_positive > 0 {
            tp as f64 / total_positive as f64
        } else {
            0.0
        };

        thresholds.push(threshold);
        tpr.push(tpr_val);
        fpr.push(fpr_val);
        precision.push(precision_val);
        recall.push(recall_val);

        roc_points.push((fpr_val, tpr_val));
    }

    // Add end point (1,1)
    roc_points.push((1.0, 1.0));

    // Directly calculate AUC using original order of roc_points
    let mut auc = 0.0;
    for i in 1..roc_points.len() {
        let (fpr_prev, tpr_prev) = roc_points[i - 1];
        let (fpr_curr, tpr_curr) = roc_points[i];
        auc += (fpr_curr - fpr_prev) * (tpr_prev + tpr_curr) / 2.0;
    }

    // Calculate other metrics - use correct confusion matrix calculation
    // Select best threshold using Youden's J statistic (Sensitivity + Specificity - 1)
    // This avoids selecting thresholds that lead to extreme results (Sensitivity=1.0, Specificity=0.0)
    let mut best_youden_j = -1.0;
    let mut best_threshold_idx: Option<usize> = None;

    for i in 0..precision.len() {
        // Calculate sensitivity and specificity for this threshold
        let sensitivity = recall[i]; // recall is TPR (sensitivity)
        let specificity = if fpr[i] < 1.0 {
            1.0 - fpr[i] // specificity = 1 - FPR
        } else {
            0.0
        };
        
        // Calculate Youden's J statistic
        let youden_j = sensitivity + specificity - 1.0;
        
        // Skip thresholds that lead to extreme results
        if (sensitivity == 1.0 && specificity == 0.0) || (sensitivity == 0.0 && specificity == 1.0) {
            continue;
        }
        
        if youden_j > best_youden_j {
            best_youden_j = youden_j;
            best_threshold_idx = Some(i);
        }
    }

    // Calculate confusion matrix using best threshold
    // If no valid threshold found (all were extreme), fall back to F1-score based selection
    let threshold = if let Some(idx) = best_threshold_idx {
        if idx < thresholds.len() {
            thresholds[idx]
        } else {
            thresholds[thresholds.len() - 1]
        }
    } else {
        // Fallback: if all thresholds were extreme, use F1-score to find best threshold
        let mut best_f1 = 0.0;
        let mut fallback_idx = 0;
        for i in 0..precision.len() {
            let f1 = if precision[i] + recall[i] > 0.0 {
                2.0 * precision[i] * recall[i] / (precision[i] + recall[i])
            } else {
                0.0
            };
            if f1 > best_f1 {
                best_f1 = f1;
                fallback_idx = i;
            }
        }
        if fallback_idx < thresholds.len() {
            thresholds[fallback_idx]
        } else {
            thresholds[thresholds.len() - 1]
        }
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

    // Calculate correct metrics
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

    logger.log(&format!(
        "Evaluation completed: AUC={:.4}, F1={:.4}, Sensitivity={:.4}, Specificity={:.4}",
        auc, f1_score, sensitivity, specificity
    ))?;

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

/// Evaluate reactivity accuracy by base separately.
/// `reactive_bases`: only these bases are used for overall AUC and combined ROC (e.g. "AC" for DMS).
pub fn evaluate_reactivity_accuracy_by_base(
    reactivity_file: &str,
    secondary_structure: &SecondaryStructure,
    signal_type: &str,
    gene_id: &str,
    strand: &str,
    output_dir: &str,
    reactive_bases: &str,
) -> Result<ComprehensiveEvaluationResult, Box<dyn Error>> {
    // Extract reactivity data and base information

    // Extract reactivity data and base information
    let reactivity_data =
        extract_reactivity_with_bases(reactivity_file, gene_id, strand, signal_type)?;

    // Match data based on bases
    let matched_data = match_reactivity_with_structure_auto_shift(
        &reactivity_data,
        secondary_structure,
        signal_type,
    )?;

    // Output intermediate process file - overall
    let base_analysis_name = format!("evaluation_by_base_{}_{}_{}", signal_type, gene_id, strand);
    let overall_tsv_path = format!("{}/{}.overall.tsv", output_dir, base_analysis_name);

    let mut overall_file = File::create(&overall_tsv_path)?;
    writeln!(
        overall_file,
        "position\tbase\treactivity\tis_single_stranded"
    )?;

    for ((pos, base, _), (reactivity, is_single)) in reactivity_data.iter().zip(matched_data.iter())
    {
        if reactive_bases.contains(*base) {
            writeln!(
                overall_file,
                "{}\t{}\t{:.6}\t{}",
                pos,
                base,
                reactivity,
                if *is_single { 1 } else { 0 }
            )?;
        }
    }

    // Group by base
    let mut base_groups: HashMap<char, Vec<(f64, bool)>> = HashMap::new();
    for ((_, base, _), (reactivity, is_single)) in reactivity_data.iter().zip(matched_data.iter()) {
        base_groups
            .entry(*base)
            .or_insert_with(Vec::new)
            .push((*reactivity, *is_single));
    }

    // Output intermediate process files for each base
    for (base, data) in &base_groups {
        let base_tsv_path = format!("{}/{}.{}.tsv", output_dir, base_analysis_name, base);

        let mut base_file = File::create(&base_tsv_path)?;
        writeln!(base_file, "position\tbase\treactivity\tis_single_stranded")?;

        // Re-traverse reactivity_data, only output sites for this base
        for ((pos, base_char, _), (reactivity, is_single)) in
            reactivity_data.iter().zip(matched_data.iter())
        {
            if *base_char == *base {
                writeln!(
                    base_file,
                    "{}\t{}\t{:.6}\t{}",
                    pos,
                    base_char,
                    reactivity,
                    if *is_single { 1 } else { 0 }
                )?;
            }
        }
    }

    // Calculate AUC based on intermediate file data

    // Read overall data and calculate AUC
    let overall_tsv_path = format!("{}/{}.overall.tsv", output_dir, base_analysis_name);
    let overall_data = read_tsv_data(&overall_tsv_path)?;
    let overall_matched_data: Vec<(f64, bool)> = overall_data
        .iter()
        .map(|(_, _, reactivity, is_single)| (*reactivity, *is_single))
        .collect();
    let overall_result =
        evaluate_reactivity_accuracy_from_matched_data(&overall_matched_data, signal_type)?;

    // Read each base data and calculate AUC (only for bases in reactive_bases)
    let mut base_specific_results = Vec::new();
    for base in ['A', 'C', 'G', 'T', 'U'] {
        if !reactive_bases.contains(base) {
            continue;
        }
        let base_tsv_path = format!("{}/{}.{}.tsv", output_dir, base_analysis_name, base);
        if std::path::Path::new(&base_tsv_path).exists() {
            match read_tsv_data(&base_tsv_path) {
                Ok(data) => {
                    if data.len() >= 10 {
                        // At least 10 sites required
                        let base_matched_data: Vec<(f64, bool)> = data
                            .iter()
                            .map(|(_, _, reactivity, is_single)| (*reactivity, *is_single))
                            .collect();
                        match evaluate_reactivity_accuracy_from_matched_data(
                            &base_matched_data,
                            signal_type,
                        ) {
                            Ok(result) => {
                                let single_count = data
                                    .iter()
                                    .filter(|(_, _, _, is_single)| *is_single)
                                    .count();
                                let double_count = data.len() - single_count;
                                let auc = result.auc;
                                base_specific_results.push(BaseSpecificEvaluationResult {
                                    base,
                                    result,
                                    data_count: data.len(),
                                    single_stranded_count: single_count,
                                    double_stranded_count: double_count,
                                });
                                // Temporarily not output, wait for unified output later
                            }
                            Err(e) => {
                                println!("Warning: {} base AUC calculation failed: {}", base, e);
                            }
                        }
                    } else {
                        println!("Warning: {} base site count insufficient ({} < 10), skipping", base, data.len());
                    }
                }
                Err(e) => {
                    println!("Warning: failed to read {} base file: {}", base, e);
                }
            }
        } else {
            // println!("Info: {} base file does not exist: {}", base, base_tsv_path);
        }
    }

    // Plot overall+4 bases 5-line ROC curve
    let combined_roc_path = format!("{}/{}.combined_roc.svg", output_dir, base_analysis_name);
    plot_multi_roc_curve(
        &overall_result,
        &base_specific_results,
        &combined_roc_path,
        &format!(
            "ROC Curve (Overall + 4 Bases) - {} ({}, Strand {})",
            signal_type, gene_id, strand
        ),
    )?;

    // Sort by base
    base_specific_results.sort_by_key(|r| r.base);

    Ok(ComprehensiveEvaluationResult {
        overall_result,
        base_specific_results,
        total_positions: matched_data.len(),
    })
}

/// Calculate evaluation results from matched data
fn evaluate_reactivity_accuracy_from_matched_data(
    matched_data: &[(f64, bool)],
    signal_type: &str,
) -> Result<EvaluationResult, Box<dyn Error>> {
    // Sort data, NaN values will be placed last (but will be filtered out during statistics)
    let mut sorted_data = matched_data.to_vec();
    sorted_data.sort_by(|a, b| {
        match (a.0.is_nan(), b.0.is_nan()) {
            (true, true) => std::cmp::Ordering::Equal,
            (true, false) => std::cmp::Ordering::Greater, // NaN placed last
            (false, true) => std::cmp::Ordering::Less,
            (false, false) => a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal),
        }
    });

    // Calculate ROC curve data
    let mut thresholds = Vec::new();
    let mut tpr = Vec::new();
    let mut fpr = Vec::new();
    let mut precision = Vec::new();
    let mut recall = Vec::new();

    // Calculate total positive and negative sample counts, filter out NaN values
    let total_positive: usize = sorted_data
        .iter()
        .filter(|(value, is_single)| !value.is_nan() && *is_single)
        .count();
    let total_negative: usize = sorted_data
        .iter()
        .filter(|(value, is_single)| !value.is_nan() && !*is_single)
        .count();

    if total_positive == 0 || total_negative == 0 {
        return Err("No positive or negative samples found".into());
    }

    // Sample statistics will be displayed in final summary

    // Generate ROC curve data (based on actual data points)
    let mut roc_points = Vec::new();

    // Add starting point (0,0)
    roc_points.push((0.0, 0.0));

    // Calculate ROC points for each unique data value, sorted by threshold from high to low, filter out NaN values
    let mut unique_values: Vec<f64> = sorted_data
        .iter()
        .map(|(val, _)| *val)
        .filter(|v| !v.is_nan())
        .collect();
    unique_values.sort_by(|a, b| {
        match (a.is_nan(), b.is_nan()) {
            (true, true) => std::cmp::Ordering::Equal,
            (true, false) => std::cmp::Ordering::Greater, // NaN placed at the end
            (false, true) => std::cmp::Ordering::Less,
            (false, false) => b.partial_cmp(a).unwrap_or(std::cmp::Ordering::Equal), // Sort from high to low
        }
    });
    unique_values.dedup(); // Remove duplicate values

    for threshold in unique_values {
        let mut tp = 0;
        let mut fp = 0;
        let mut tn = 0;
        let mut fn_count = 0;

        for (value, is_single) in &sorted_data {
            // Skip NaN values
            if value.is_nan() {
                continue;
            }
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

        let tpr_val = if total_positive > 0 {
            tp as f64 / total_positive as f64
        } else {
            0.0
        };
        let fpr_val = if total_negative > 0 {
            fp as f64 / total_negative as f64
        } else {
            0.0
        };
        let precision_val = if (tp + fp) > 0 {
            tp as f64 / (tp + fp) as f64
        } else {
            0.0
        };
        let recall_val = if total_positive > 0 {
            tp as f64 / total_positive as f64
        } else {
            0.0
        };

        thresholds.push(threshold);
        tpr.push(tpr_val);
        fpr.push(fpr_val);
        precision.push(precision_val);
        recall.push(recall_val);

        roc_points.push((fpr_val, tpr_val));
    }

    // Add ending point (1,1)
    roc_points.push((1.0, 1.0));

    // Correct ROC curve generation: sort by FPR, ensure monotonicity
    let mut sorted_roc_points = roc_points.clone();
    sorted_roc_points.sort_by(|a, b| {
        // FPR and TPR should not be NaN, but handle it for safety
        match (a.0.is_nan(), b.0.is_nan()) {
            (true, true) => std::cmp::Ordering::Equal,
            (true, false) => std::cmp::Ordering::Greater,
            (false, true) => std::cmp::Ordering::Less,
            (false, false) => a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal),
        }
    }); // Sort by FPR

    // Ensure monotonicity of ROC points: if FPR is the same, take maximum TPR
    let mut monotonic_roc_points = Vec::new();
    monotonic_roc_points.push((0.0, 0.0));

    let mut current_fpr = 0.0;
    let mut max_tpr_at_current_fpr: f64 = 0.0;

    for (fpr, tpr) in sorted_roc_points.iter().skip(1) {
        // Skip starting point (0,0)
        if (fpr - current_fpr).abs() < 1e-10 {
            // Same FPR
            max_tpr_at_current_fpr = max_tpr_at_current_fpr.max(*tpr);
        } else {
            // New FPR value
            if current_fpr > 0.0 {
                // Save maximum TPR for previous FPR
                monotonic_roc_points.push((current_fpr, max_tpr_at_current_fpr));
            }
            current_fpr = *fpr;
            max_tpr_at_current_fpr = *tpr;
        }
    }

    // Add last point
    if current_fpr > 0.0 {
        monotonic_roc_points.push((current_fpr, max_tpr_at_current_fpr));
    }
    monotonic_roc_points.push((1.0, 1.0));

    // Calculate AUC using monotonic ROC points
    let mut auc = 0.0;
    for i in 1..monotonic_roc_points.len() {
        let (fpr_prev, tpr_prev) = monotonic_roc_points[i - 1];
        let (fpr_curr, tpr_curr) = monotonic_roc_points[i];
        auc += (fpr_curr - fpr_prev) * (tpr_prev + tpr_curr) / 2.0;
    }

    // Recalculate TPR and FPR arrays to ensure consistency with monotonic ROC points
    let mut corrected_tpr = Vec::new();
    let mut corrected_fpr = Vec::new();
    for (fpr, tpr) in &monotonic_roc_points[1..monotonic_roc_points.len() - 1] {
        corrected_fpr.push(*fpr);
        corrected_tpr.push(*tpr);
    }

    // Calculate other metrics - use correct confusion matrix calculation
    // Select best threshold using Youden's J statistic (Sensitivity + Specificity - 1)
    // This avoids selecting thresholds that lead to extreme results (Sensitivity=1.0, Specificity=0.0)
    let mut best_youden_j = -1.0;
    let mut best_threshold_idx: Option<usize> = None;

    for i in 0..precision.len() {
        // Calculate sensitivity and specificity for this threshold
        let sensitivity = recall[i]; // recall is TPR (sensitivity)
        let specificity = if fpr.len() > i && fpr[i] < 1.0 {
            1.0 - fpr[i] // specificity = 1 - FPR
        } else {
            0.0
        };
        
        // Calculate Youden's J statistic
        let youden_j = sensitivity + specificity - 1.0;
        
        // Skip thresholds that lead to extreme results
        if (sensitivity == 1.0 && specificity == 0.0) || (sensitivity == 0.0 && specificity == 1.0) {
            continue;
        }
        
        if youden_j > best_youden_j {
            best_youden_j = youden_j;
            best_threshold_idx = Some(i);
        }
    }

    // Calculate confusion matrix using best threshold
    // If no valid threshold found (all were extreme), fall back to F1-score based selection
    let threshold = if let Some(idx) = best_threshold_idx {
        if idx < thresholds.len() {
            thresholds[idx]
        } else {
            thresholds[thresholds.len() - 1]
        }
    } else {
        // Fallback: if all thresholds were extreme, use F1-score to find best threshold
        let mut best_f1 = 0.0;
        let mut fallback_idx = 0;
        for i in 0..precision.len() {
            let f1 = if precision[i] + recall[i] > 0.0 {
                2.0 * precision[i] * recall[i] / (precision[i] + recall[i])
            } else {
                0.0
            };
            if f1 > best_f1 {
                best_f1 = f1;
                fallback_idx = i;
            }
        }
        if fallback_idx < thresholds.len() {
            thresholds[fallback_idx]
        } else {
            thresholds[thresholds.len() - 1]
        }
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

    // Calculate correct metrics
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

/// Save comprehensive evaluation results
fn save_comprehensive_evaluation_results(
    result: &ComprehensiveEvaluationResult,
    output_path: &str,
    signal_type: &str,
    gene_id: &str,
    strand: &str,
) -> Result<(), Box<dyn Error>> {
    let mut file = File::create(output_path)?;
    writeln!(file, "=== Evaluation Results by Base ===")?;
    writeln!(file, "Gene ID: {}", gene_id)?;
    writeln!(file, "Strand: {}", strand)?;
    writeln!(file, "Signal Type: {}", signal_type)?;
    writeln!(file, "Total Positions: {}", result.total_positions)?;
    writeln!(file)?;

    // Overall evaluation results
    writeln!(file, "=== Overall Evaluation Results ===")?;
    writeln!(file, "AUC: {:.4}", result.overall_result.auc)?;
    writeln!(file, "F1-score: {:.4}", result.overall_result.f1_score)?;
    writeln!(
        file,
        "Sensitivity: {:.4}",
        result.overall_result.sensitivity
    )?;
    writeln!(
        file,
        "Specificity: {:.4}",
        result.overall_result.specificity
    )?;
    writeln!(file, "PPV: {:.4}", result.overall_result.ppv)?;
    writeln!(file, "NPV: {:.4}", result.overall_result.npv)?;
    writeln!(file, "Accuracy: {:.4}", result.overall_result.accuracy)?;
    writeln!(file)?;

    // Base-specific evaluation results
    writeln!(file, "=== Base-Specific Evaluation Results ===")?;
    for base_result in &result.base_specific_results {
        writeln!(file, "--- {} Base ---", base_result.base)?;
        writeln!(
            file,
            "Position Count: {} (Single-stranded: {}, Double-stranded: {})",
            base_result.data_count,
            base_result.single_stranded_count,
            base_result.double_stranded_count
        )?;
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

/// Save base-specific evaluation results
fn save_base_specific_evaluation_results(
    result: &BaseSpecificEvaluationResult,
    output_path: &str,
    signal_type: &str,
) -> Result<(), Box<dyn Error>> {
    let mut file = File::create(output_path)?;
    writeln!(file, "=== {} Base Evaluation Results ===", result.base)?;
    writeln!(file, "Signal Type: {}", signal_type)?;
    writeln!(
        file,
        "Position Count: {} (Single-stranded: {}, Double-stranded: {})",
        result.data_count, result.single_stranded_count, result.double_stranded_count
    )?;
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

/// Plot multiple ROC curves (overall + 4 bases)
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
        .x_label_area_size(50) // Increase label area
        .y_label_area_size(50) // Increase label area
        .build_cartesian_2d(0.0..1.0, 0.0..1.0)?;

    // Completely remove grid lines, only keep axes
    chart
        .configure_mesh()
        .x_desc("False Positive Rate")
        .y_desc("True Positive Rate")
        .x_labels(5) // Reduce number of x-axis tick marks
        .y_labels(5) // Reduce number of y-axis tick marks
        .disable_mesh() // Completely disable grid lines
        .draw()?;

    // Plot overall
    let mut roc_points = vec![(0.0, 0.0)];
    for i in 0..overall.fpr.len() {
        roc_points.push((overall.fpr[i], overall.tpr[i]));
    }
    roc_points.push((1.0, 1.0));
    chart
        .draw_series(LineSeries::new(
            roc_points,
            RGBColor(0, 0, 0).stroke_width(3),
        ))?
        .label(format!("Overall (AUC={:.3})", overall.auc))
        .legend(move |(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], RGBColor(0, 0, 0)));

    // Iterate through all existing bases
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
        chart
            .draw_series(LineSeries::new(roc_points, color.stroke_width(2)))?
            .label(format!("{} (AUC={:.3})", base_char, base_result.result.auc))
            .legend(move |(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], color));
    }

    // Draw diagonal line
    chart.draw_series(LineSeries::new(
        vec![(0.0, 0.0), (1.0, 1.0)],
        RGBColor(180, 180, 180).stroke_width(1),
    ))?;

    // Add AUC annotation (increase font size)
    chart.draw_series(std::iter::once(Text::new(
        format!("Overall AUC = {:.3}", overall.auc),
        (0.6, 0.2),
        ("sans-serif", 20).into_font().color(&RGBColor(0, 0, 0)),
    )))?;

    // Legend (increase font size, place in lower right corner)
    chart
        .configure_series_labels()
        .background_style(WHITE.mix(0.8))
        .border_style(BLACK)
        .label_font(("sans-serif", 18)) // Further increase legend font size
        .position(SeriesLabelPosition::LowerRight) // Place in lower right corner
        .draw()?;

    root.present()?;
    Ok(())
}

/// Read TSV format intermediate process file data
fn read_tsv_data(file_path: &str) -> Result<Vec<(u32, char, f64, bool)>, Box<dyn Error>> {
    let mut data = Vec::new();
    let file = File::open(file_path)?;
    let reader = BufReader::new(file);

    for (line_num, line) in reader.lines().enumerate() {
        if line_num == 0 {
            continue; // Skip header row
        }

        let line = line?;
        let fields: Vec<&str> = line.split('\t').collect();

        if fields.len() >= 4 {
            let position = fields[0].parse::<u32>()?;
            let base_str = fields[1].trim();
            if base_str.is_empty() {
                continue; // Skip empty base
            }
            let base = base_str.chars().next().unwrap_or('N');
            let reactivity = parse_f64_allow_nan(fields[2]).ok_or_else(|| format!("Invalid reactivity value: {}", fields[2]))?;
            let is_single = fields[3].parse::<u32>()? == 1;

            data.push((position, base, reactivity, is_single));
        }
    }

    Ok(data)
}

/// Optimized main evaluation function (only processes target genes)
pub fn evaluate_reactivity_accuracy_main_optimized(
    reactivity_file: &str,
    secondary_structure_file: &str,
    output_dir: &str,
    signal_type: &str,
    strand: &str,
    logger: &mut crate::Logger,
) -> Result<(), Box<dyn Error>> {
    let start_time = Instant::now();

    // Record environment information and parameters
    logger.log("=== ModDetector Evaluate Function Log ===")?;
    logger.log(&format!("Software Version: v{}", crate::VERSION))?;
    logger.log(&format!(
        "Runtime: {}",
        chrono::Utc::now().format("%Y-%m-%d %H:%M:%S")
    ))?;
    logger.log(&format!("Reactivity File: {}", reactivity_file))?;
    logger.log(&format!(
        "Secondary Structure File: {}",
        secondary_structure_file
    ))?;
    logger.log(&format!("Output Directory: {}", output_dir))?;
    logger.log(&format!("Signal Type: {}", signal_type))?;
    logger.log(&format!("Strand: {}", strand))?;
    logger.log("Starting optimized reactivity accuracy evaluation...")?;

    // Define target gene patterns
    let target_gene_patterns = ["16S", "J01695", "EC", "E.coli"];

    // Directly extract target gene data
    logger.log_and_progress("Extracting target gene data...")?;
    let gene_data = extract_reactivity_for_target_gene(
        reactivity_file,
        &target_gene_patterns,
        strand,
        signal_type,
    )?;

    match gene_data {
        Some((gene_id, reactivity_data)) => {
            logger.log_and_progress(&format!("Found target gene: {}", gene_id))?;

            // Parse secondary structure
            let secondary_structure = parse_secondary_structure(secondary_structure_file)?;

            // Create output directory
            std::fs::create_dir_all(output_dir)?;

            // Perform accuracy evaluation with auto-shift correction
            logger
                .log_and_progress("Performing accuracy evaluation with auto-shift correction...")?;
            let matched_data = match_reactivity_with_structure_auto_shift(
                &reactivity_data,
                &secondary_structure,
                signal_type,
            )?;

            // Calculate evaluation results
            let result =
                evaluate_reactivity_accuracy_from_matched_data(&matched_data, signal_type)?;

            // Create gene output directory
            let gene_output_dir = format!("{}/{}", output_dir, gene_id.replace('/', "_"));
            std::fs::create_dir_all(&gene_output_dir)?;

            // Save evaluation results
            let result_file = format!("{}/evaluation_results_{}.txt", gene_output_dir, signal_type);
            save_evaluation_results(&result, &result_file, signal_type)?;

            // Plot ROC curve
            let roc_file = format!("{}/roc_curve_{}.svg", gene_output_dir, signal_type);
            let roc_title = format!("{} - {} Signal ROC Curve", gene_id, signal_type);
            plot_roc_curve(&result, &roc_file, &roc_title)?;

            // Plot PR curve
            let pr_file = format!("{}/pr_curve_{}.svg", gene_output_dir, signal_type);
            let pr_title = format!(
                "{} - {} Signal Precision-Recall Curve",
                gene_id, signal_type
            );
            plot_pr_curve(&result, &pr_file, &pr_title)?;

            // Output main metrics
            logger.log(&format!("  AUC: {:.3}", result.auc))?;
            logger.log(&format!("  F1-score: {:.3}", result.f1_score))?;
            logger.log(&format!("  Sensitivity: {:.3}", result.sensitivity))?;
            logger.log(&format!("  Specificity: {:.3}", result.specificity))?;
            logger.log(&format!("  PPV: {:.3}", result.ppv))?;
            logger.log(&format!("  NPV: {:.3}", result.npv))?;
            logger.log(&format!("  Accuracy: {:.3}", result.accuracy))?;
            logger.log(&format!("  Matched positions: {}", matched_data.len()))?;

            logger.log(&format!(
                "Successfully processed gene: {}, evaluation completed",
                gene_id
            ))?;
        }
        None => {
            logger.log("Warning: No matching target gene found for evaluation")?;
        }
    }

    let elapsed = start_time.elapsed();
    println!("{}", crate::progress::format_time_used(elapsed));
    Ok(())
}

/// Evaluate accuracy of both stop and mutation signals simultaneously
pub fn evaluate_both_signals_main(
    reactivity_file: &str,
    secondary_structure_file: &str,
    output_dir: &str,
    gene_id: &str,
    strand: &str,
    reactive_bases: &str,
    logger: &mut crate::Logger,
) -> Result<(), Box<dyn Error>> {
    let start_time = Instant::now();

    // Record environment information and parameters
    logger.log("=== ModDetector Evaluate Both Signals Function Log ===")?;
    logger.log(&format!("Software Version: v{}", crate::VERSION))?;
    logger.log(&format!(
        "Runtime: {}",
        chrono::Utc::now().format("%Y-%m-%d %H:%M:%S")
    ))?;
    logger.log(&format!("Reactivity File: {}", reactivity_file))?;
    logger.log(&format!(
        "Secondary Structure File: {}",
        secondary_structure_file
    ))?;
    logger.log(&format!("Output Directory: {}", output_dir))?;
    logger.log(&format!("Gene ID: {}", gene_id))?;
    logger.log(&format!("Strand: {}", strand))?;
    logger.log("Starting both signals evaluation (stop and mutation)...")?;

    // Display data loading information
    println!("[Loading data]");
    println!("    Secondary structure: {}", secondary_structure_file);
    println!("    Reactivity file: {}", reactivity_file);
    println!();

    // Display parameter information
    println!("[Params]");
    println!("    Gene ID: {}, strand: {}", gene_id, strand);
    println!("    Signal types: stop, mutation");
    println!("    Use base matching: true");
    println!("    Auto shift: true");
    println!("    Reactive bases (Overall/ROC): {}", reactive_bases);
    println!();

    // Create output directory
    std::fs::create_dir_all(output_dir)?;

    // Parse secondary structure file
    let secondary_structure = parse_secondary_structure(secondary_structure_file)?;

    // Display progress information
    println!("[Progressing]");
    println!(
        "    Structure positions: {}, unpaired: {}, paired: {}",
        secondary_structure.positions.len(),
        secondary_structure
            .is_single_stranded
            .iter()
            .filter(|&&x| x)
            .count(),
        secondary_structure
            .is_single_stranded
            .iter()
            .filter(|&&x| !x)
            .count()
    );
    println!();

    // Evaluate stop signal
    println!("[Evaluating] Stop signal...");
    let stop_result = evaluate_reactivity_accuracy_with_auto_shift(
        reactivity_file,
        &secondary_structure,
        "stop",
        gene_id,
        strand,
        logger,
    )?;

    // Output clear report for stop signal
    println!("[Stop Signal Results]");
    println!("    AUC: {:.3}", stop_result.auc);
    println!("    F1-score: {:.3}", stop_result.f1_score);
    println!("    Sensitivity: {:.3}", stop_result.sensitivity);
    println!("    Specificity: {:.3}", stop_result.specificity);
    println!("    PPV: {:.3}", stop_result.ppv);
    println!("    NPV: {:.3}", stop_result.npv);
    println!("    Accuracy: {:.3}", stop_result.accuracy);
    println!();

    // Evaluate mutation signal
    println!("[Evaluating] Mutation signal...");
    let mutation_result = evaluate_reactivity_accuracy_with_auto_shift(
        reactivity_file,
        &secondary_structure,
        "mutation",
        gene_id,
        strand,
        logger,
    )?;

    // Output clear report for mutation signal
    println!("[Mutation Signal Results]");
    println!("    AUC: {:.3}", mutation_result.auc);
    println!("    F1-score: {:.3}", mutation_result.f1_score);
    println!("    Sensitivity: {:.3}", mutation_result.sensitivity);
    println!("    Specificity: {:.3}", mutation_result.specificity);
    println!("    PPV: {:.3}", mutation_result.ppv);
    println!("    NPV: {:.3}", mutation_result.npv);
    println!("    Accuracy: {:.3}", mutation_result.accuracy);
    println!();

    // Performance comparison
    println!("[Performance Comparison]");
    if stop_result.auc > mutation_result.auc {
        println!(
            "    Stop signal performs better (AUC: {:.3} vs {:.3})",
            stop_result.auc, mutation_result.auc
        );
    } else if mutation_result.auc > stop_result.auc {
        println!(
            "    Mutation signal performs better (AUC: {:.3} vs {:.3})",
            mutation_result.auc, stop_result.auc
        );
    } else {
        println!(
            "    Both signals perform equally (AUC: {:.3})",
            stop_result.auc
        );
    }
    println!();

    // Evaluate stop signal by base separately
    println!("[Evaluating] Stop signal by base...");
    let stop_comprehensive_result = evaluate_reactivity_accuracy_by_base(
        reactivity_file,
        &secondary_structure,
        "stop",
        gene_id,
        strand,
        output_dir,
        reactive_bases,
    )?;

    // Output clear report for stop signal by base
    println!("[Stop Signal Base-Specific Results]");
    println!(
        "    Overall AUC: {:.3}",
        stop_comprehensive_result.overall_result.auc
    );
    for base_result in &stop_comprehensive_result.base_specific_results {
        println!(
            "    {} base AUC: {:.3} (Positions: {}, Unpaired: {}, Paired: {})",
            base_result.base,
            base_result.result.auc,
            base_result.data_count,
            base_result.single_stranded_count,
            base_result.double_stranded_count
        );
    }
    println!();

    // Evaluate mutation signal by base separately
    println!("[Evaluating] Mutation signal by base...");
    let mutation_comprehensive_result = evaluate_reactivity_accuracy_by_base(
        reactivity_file,
        &secondary_structure,
        "mutation",
        gene_id,
        strand,
        output_dir,
        reactive_bases,
    )?;

    // Output clear report for mutation signal by base
    println!("[Mutation Signal Base-Specific Results]");
    println!(
        "    Overall AUC: {:.3}",
        mutation_comprehensive_result.overall_result.auc
    );
    for base_result in &mutation_comprehensive_result.base_specific_results {
        println!(
            "    {} base AUC: {:.3} (Positions: {}, Unpaired: {}, Paired: {})",
            base_result.base,
            base_result.result.auc,
            base_result.data_count,
            base_result.single_stranded_count,
            base_result.double_stranded_count
        );
    }
    println!();

    // Save stop signal results
    println!("[Saving] Stop signal results...");
    save_evaluation_results(
        &stop_result,
        &format!("{}/evaluation_stop.txt", output_dir),
        "stop",
    )?;
    plot_roc_curve(
        &stop_result,
        &format!("{}/roc_curve_stop.svg", output_dir),
        &format!("{} - Stop Signal ROC Curve", gene_id),
    )?;
    plot_pr_curve(
        &stop_result,
        &format!("{}/pr_curve_stop.svg", output_dir),
        &format!("{} - Stop Signal PR Curve", gene_id),
    )?;

    // Save mutation signal results
    println!("[Saving] Mutation signal results...");
    save_evaluation_results(
        &mutation_result,
        &format!("{}/evaluation_mutation.txt", output_dir),
        "mutation",
    )?;
    plot_roc_curve(
        &mutation_result,
        &format!("{}/roc_curve_mutation.svg", output_dir),
        &format!("{} - Mutation Signal ROC Curve", gene_id),
    )?;
    plot_pr_curve(
        &mutation_result,
        &format!("{}/pr_curve_mutation.svg", output_dir),
        &format!("{} - Mutation Signal PR Curve", gene_id),
    )?;

    // Save base-specific evaluation results for stop signal
    println!("[Saving] Stop signal base-specific results...");
    let stop_base_analysis_name = format!("evaluation_by_base_stop_{}_{}", gene_id, strand);
    let stop_overall_roc_path =
        format!("{}/{}.overall.roc.svg", output_dir, stop_base_analysis_name);
    let stop_overall_pr_path = format!("{}/{}.overall.pr.svg", output_dir, stop_base_analysis_name);
    let stop_overall_results_path =
        format!("{}/{}.overall.txt", output_dir, stop_base_analysis_name);
    let stop_comprehensive_results_path = format!(
        "{}/{}.comprehensive.txt",
        output_dir, stop_base_analysis_name
    );

    plot_roc_curve(
        &stop_comprehensive_result.overall_result,
        &stop_overall_roc_path,
        &format!("Overall ROC Curve - Stop ({}, Strand {})", gene_id, strand),
    )?;
    plot_pr_curve(
        &stop_comprehensive_result.overall_result,
        &stop_overall_pr_path,
        &format!("Overall PR Curve - Stop ({}, Strand {})", gene_id, strand),
    )?;
    save_evaluation_results(
        &stop_comprehensive_result.overall_result,
        &stop_overall_results_path,
        "stop",
    )?;
    save_comprehensive_evaluation_results(
        &stop_comprehensive_result,
        &stop_comprehensive_results_path,
        "stop",
        gene_id,
        strand,
    )?;

    // Plot ROC curve for each base of stop signal
    for base_result in &stop_comprehensive_result.base_specific_results {
        let base_roc_path = format!(
            "{}/{}.{}.roc.svg",
            output_dir, stop_base_analysis_name, base_result.base
        );
        let base_pr_path = format!(
            "{}/{}.{}.pr.svg",
            output_dir, stop_base_analysis_name, base_result.base
        );
        let base_results_path = format!(
            "{}/{}.{}.txt",
            output_dir, stop_base_analysis_name, base_result.base
        );

        plot_roc_curve(
            &base_result.result,
            &base_roc_path,
            &format!(
                "{} Base Stop ROC Curve - {} (Strand {})",
                base_result.base, gene_id, strand
            ),
        )?;
        plot_pr_curve(
            &base_result.result,
            &base_pr_path,
            &format!(
                "{} Base Stop PR Curve - {} (Strand {})",
                base_result.base, gene_id, strand
            ),
        )?;
        save_base_specific_evaluation_results(base_result, &base_results_path, "stop")?;
    }

    // Plot combined ROC curve for stop signal
    let stop_combined_roc_path = format!("{}/{}.roc.svg", output_dir, stop_base_analysis_name);
    plot_multi_roc_curve(
        &stop_comprehensive_result.overall_result,
        &stop_comprehensive_result.base_specific_results,
        &stop_combined_roc_path,
        &format!(
            "ROC Curve (Overall + 4 Bases) - Stop ({}, Strand {})",
            gene_id, strand
        ),
    )?;

    // Save base-specific evaluation results for mutation signal
    println!("[Saving] Mutation signal base-specific results...");
    let mutation_base_analysis_name = format!("evaluation_by_base_mutation_{}_{}", gene_id, strand);
    let mutation_overall_roc_path = format!(
        "{}/{}.overall.roc.svg",
        output_dir, mutation_base_analysis_name
    );
    let mutation_overall_pr_path = format!(
        "{}/{}.overall.pr.svg",
        output_dir, mutation_base_analysis_name
    );
    let mutation_overall_results_path =
        format!("{}/{}.overall.txt", output_dir, mutation_base_analysis_name);
    let mutation_comprehensive_results_path = format!(
        "{}/{}.comprehensive.txt",
        output_dir, mutation_base_analysis_name
    );

    plot_roc_curve(
        &mutation_comprehensive_result.overall_result,
        &mutation_overall_roc_path,
        &format!(
            "Overall ROC Curve - Mutation ({}, Strand {})",
            gene_id, strand
        ),
    )?;
    plot_pr_curve(
        &mutation_comprehensive_result.overall_result,
        &mutation_overall_pr_path,
        &format!(
            "Overall PR Curve - Mutation ({}, Strand {})",
            gene_id, strand
        ),
    )?;
    save_evaluation_results(
        &mutation_comprehensive_result.overall_result,
        &mutation_overall_results_path,
        "mutation",
    )?;
    save_comprehensive_evaluation_results(
        &mutation_comprehensive_result,
        &mutation_comprehensive_results_path,
        "mutation",
        gene_id,
        strand,
    )?;

    // Plot ROC curve for each base of mutation signal
    for base_result in &mutation_comprehensive_result.base_specific_results {
        let base_roc_path = format!(
            "{}/{}.{}.roc.svg",
            output_dir, mutation_base_analysis_name, base_result.base
        );
        let base_pr_path = format!(
            "{}/{}.{}.pr.svg",
            output_dir, mutation_base_analysis_name, base_result.base
        );
        let base_results_path = format!(
            "{}/{}.{}.txt",
            output_dir, mutation_base_analysis_name, base_result.base
        );

        plot_roc_curve(
            &base_result.result,
            &base_roc_path,
            &format!(
                "{} Base Mutation ROC Curve - {} (Strand {})",
                base_result.base, gene_id, strand
            ),
        )?;
        plot_pr_curve(
            &base_result.result,
            &base_pr_path,
            &format!(
                "{} Base Mutation PR Curve - {} (Strand {})",
                base_result.base, gene_id, strand
            ),
        )?;
        save_base_specific_evaluation_results(base_result, &base_results_path, "mutation")?;
    }

    // Plot combined ROC curve for mutation signal
    let mutation_combined_roc_path =
        format!("{}/{}.roc.svg", output_dir, mutation_base_analysis_name);
    plot_multi_roc_curve(
        &mutation_comprehensive_result.overall_result,
        &mutation_comprehensive_result.base_specific_results,
        &mutation_combined_roc_path,
        &format!(
            "ROC Curve (Overall + 4 Bases) - Mutation ({}, Strand {})",
            gene_id, strand
        ),
    )?;

    // Add data information
    println!("[Data info]");
    println!(
        "    Total AUC={:.3}. Positions={}, unpaired={}, paired={}.",
        stop_comprehensive_result.overall_result.auc,
        stop_comprehensive_result.total_positions,
        stop_comprehensive_result.total_positions
            - stop_comprehensive_result
                .base_specific_results
                .iter()
                .map(|r| r.double_stranded_count)
                .sum::<usize>(),
        stop_comprehensive_result
            .base_specific_results
            .iter()
            .map(|r| r.double_stranded_count)
            .sum::<usize>()
    );

    for base_result in &stop_comprehensive_result.base_specific_results {
        println!(
            "    {} base AUC={:.3}. Positions={}, unpaired={}, paired={}.",
            base_result.base,
            base_result.result.auc,
            base_result.data_count,
            base_result.single_stranded_count,
            base_result.double_stranded_count
        );
    }
    println!();

    // Add output information
    println!("[Outputs]");
    println!(
        "    Temperat: {}/{}.overall.tsv",
        output_dir, stop_base_analysis_name
    );
    println!(
        "    Summary: {}/{}.comprehensive.txt",
        output_dir, stop_base_analysis_name
    );
    println!(
        "    ROC plots: {}/{}.combined_roc.svg",
        output_dir, stop_base_analysis_name
    );
    println!();

    // Generate comprehensive report
    println!("[Generating] Combined report...");
    let combined_report_path = format!("{}/combined_evaluation_report.txt", output_dir);
    let mut report_file = File::create(&combined_report_path)?;

    writeln!(
        report_file,
        "=== ModDetector Combined Evaluation Report ==="
    )?;
    writeln!(report_file, "Gene ID: {}", gene_id)?;
    writeln!(report_file, "Strand: {}", strand)?;
    writeln!(
        report_file,
        "Secondary Structure File: {}",
        secondary_structure_file
    )?;
    writeln!(report_file, "Reactivity File: {}", reactivity_file)?;
    writeln!(
        report_file,
        "Evaluation Time: {}",
        chrono::Utc::now().format("%Y-%m-%d %H:%M:%S")
    )?;
    writeln!(report_file)?;

    // Stop signal results
    writeln!(report_file, "=== Stop Signal Results ===")?;
    writeln!(report_file, "AUC: {:.4}", stop_result.auc)?;
    writeln!(report_file, "F1-Score: {:.4}", stop_result.f1_score)?;
    writeln!(report_file, "Sensitivity: {:.4}", stop_result.sensitivity)?;
    writeln!(report_file, "Specificity: {:.4}", stop_result.specificity)?;
    writeln!(report_file, "PPV: {:.4}", stop_result.ppv)?;
    writeln!(report_file, "NPV: {:.4}", stop_result.npv)?;
    writeln!(report_file, "Accuracy: {:.4}", stop_result.accuracy)?;
    writeln!(report_file)?;

    // Mutation signal results
    writeln!(report_file, "=== Mutation Signal Results ===")?;
    writeln!(report_file, "AUC: {:.4}", mutation_result.auc)?;
    writeln!(report_file, "F1-Score: {:.4}", mutation_result.f1_score)?;
    writeln!(
        report_file,
        "Sensitivity: {:.4}",
        mutation_result.sensitivity
    )?;
    writeln!(
        report_file,
        "Specificity: {:.4}",
        mutation_result.specificity
    )?;
    writeln!(report_file, "PPV: {:.4}", mutation_result.ppv)?;
    writeln!(report_file, "NPV: {:.4}", mutation_result.npv)?;
    writeln!(report_file, "Accuracy: {:.4}", mutation_result.accuracy)?;
    writeln!(report_file)?;

    // Performance comparison
    writeln!(report_file, "=== Performance Comparison ===")?;
    if stop_result.auc > mutation_result.auc {
        writeln!(
            report_file,
            "Stop signal performs better (AUC: {:.4} vs {:.4})",
            stop_result.auc, mutation_result.auc
        )?;
    } else if mutation_result.auc > stop_result.auc {
        writeln!(
            report_file,
            "Mutation signal performs better (AUC: {:.4} vs {:.4})",
            mutation_result.auc, stop_result.auc
        )?;
    } else {
        writeln!(
            report_file,
            "Both signals perform equally (AUC: {:.4})",
            stop_result.auc
        )?;
    }

    let elapsed = start_time.elapsed();

    // Overwrite progress display, show output information
    println!("\r[Output]                           ");
    println!("    Combined results: {}", output_dir);
    println!("    Stop signal: roc_curve_stop.svg, pr_curve_stop.svg");
    println!("    Mutation signal: roc_curve_mutation.svg, pr_curve_mutation.svg");
    println!("    Report: combined_evaluation_report.txt");
    println!("{}", crate::progress::format_time_used(elapsed));

    // Log completion
    logger.log(&format!(
        "Both signals evaluation completed, output directory: {}",
        output_dir
    ))?;
    logger.log(&format!("Stop signal AUC: {:.4}", stop_result.auc))?;
    logger.log(&format!("Mutation signal AUC: {:.4}", mutation_result.auc))?;
    logger.log(&format!("Total time: {:.2}s", elapsed.as_secs_f64()))?;

    Ok(())
}
