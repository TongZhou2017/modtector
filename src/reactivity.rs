use chrono;
use clap::Args;
use rayon::prelude::*;
use std::collections::HashMap;
use std::error::Error;
use std::fs::File;
use std::io::{BufRead, BufReader, Write, Read};
use std::path::Path;
use std::time::Instant;

// Import SecondaryStructure for recursive k-factor calculation
use crate::evaluate::{SecondaryStructure, parse_secondary_structure};

/// Validate reactivity command parameters
fn validate_reactivity_args(args: &ReactivityArgs) -> Result<(), Box<dyn Error>> {
    // Validate modified sample file
    if args.mod_csv.trim().is_empty() {
        return Err("Error: Modified sample file path cannot be empty".into());
    }
    if !Path::new(&args.mod_csv).exists() {
        return Err(format!(
            "Error: Modified sample file does not exist: {}",
            args.mod_csv
        )
        .into());
    }
    if !args.mod_csv.ends_with(".csv") {
        return Err(format!(
            "Error: Modified sample file path must end with .csv: {}",
            args.mod_csv
        )
        .into());
    }

    // Validate unmodified sample file (optional)
    if let Some(unmod_path) = args.unmod_csv.as_ref() {
        if unmod_path.trim().is_empty() {
            return Err("Error: Unmodified sample file path cannot be empty".into());
        }
        if !Path::new(unmod_path).exists() {
            return Err(format!(
                "Error: Unmodified sample file does not exist: {}",
                unmod_path
            )
            .into());
        }
        if !unmod_path.ends_with(".csv") {
            return Err(format!(
                "Error: Unmodified sample file path must end with .csv: {}",
                unmod_path
            )
            .into());
        }
    }

    // Validate output file
    if args.output.trim().is_empty() {
        return Err("Error: reactivity output file path cannot be empty".into());
    }
    if !args.output.ends_with(".csv") {
        return Err(format!(
            "Error: reactivity output file path must end with .csv: {}",
            args.output
        )
        .into());
    }

    // Validate stop method
    match args.stop_method.to_lowercase().as_str() {
        "kfactor" | "ding" | "rouskin" => {}
        _ => {
            return Err(format!(
                "Error: unknown stop signal method: {}. Supported methods: kfactor, ding, rouskin",
                args.stop_method
            )
            .into())
        }
    }

    // Validate mutation method
    match args.mutation_method.to_lowercase().as_str() {
        "kfactor" | "siegfried" | "zubradt" => {}
        _ => {
            return Err(format!(
                "Error: unknown mutation signal method: {}. Supported methods: kfactor, siegfried, zubradt",
                args.mutation_method
            )
            .into())
        }
    }

    // Validate Ding method related parameters
    if args.stop_method.to_lowercase() == "ding" {
        if args.pseudocount <= 0.0 {
            return Err("Error: Ding method pseudocount parameter must be greater than 0".into());
        }
        if args.maxscore <= 0.0 {
            return Err("Error: Ding method maxscore parameter must be greater than 0".into());
        }
        if args.maxscore > 100.0 {
            return Err(format!(
                "Error: Ding method maxscore parameter cannot exceed 100, current: {}",
                args.maxscore
            )
            .into());
        }
    }

    // Validate SNP cutoff parameter
    if args.snp_cutoff < 0.0 || args.snp_cutoff > 1.0 {
        return Err(format!(
            "Error: SNP cutoff must be between 0.0 and 1.0, current: {}",
            args.snp_cutoff
        )
        .into());
    }

    // Validate k prediction method
    match args.k_prediction_method.to_lowercase().as_str() {
        "background" | "distribution" | "recursive" => {}
        _ => {
            return Err(format!(
                "Error: unknown k prediction method: {}. Supported methods: background, distribution, recursive",
                args.k_prediction_method
            )
            .into())
        }
    }

    // Validate structure file for recursive method
    if args.k_prediction_method.to_lowercase() == "recursive" {
        if args.structure_file.is_none() {
            return Err("Error: recursive k prediction method requires --structure-file parameter".into());
        }
        if let Some(ref struct_file) = args.structure_file {
            if !Path::new(struct_file).exists() {
                return Err(format!(
                    "Error: structure file does not exist: {}",
                    struct_file
                )
                .into());
            }
        }
    }

    // Validate thread count
    if let Some(threads) = args.threads {
        if threads == 0 {
            return Err("Error: thread count cannot be 0".into());
        }
        // Remove hard limit, only provide suggestion
        if threads > 128 {
            eprintln!("Warning: thread count {} may exceed optimal range, suggest using 64-128 threads for best performance", threads);
        }
    }

    Ok(())
}

#[derive(Args, Debug)]
pub struct ReactivityArgs {
    // Input files
    /// Modified sample CSV
    #[arg(short = 'M', long = "mod")]
    pub mod_csv: String,
    /// Unmodified sample CSV
    #[arg(short = 'U', long = "unmod")]
    pub unmod_csv: Option<String>,

    // Output files
    /// reactivity output file (containing stop and mutation signals)
    #[arg(short = 'O', long = "output")]
    pub output: String,

    // Method selection
    /// stop signal reactivity method: kfactor, ding, rouskin
    #[arg(short = 's', long = "stop-method", default_value = "kfactor")]
    pub stop_method: String,
    /// mutation signal reactivity method: kfactor, siegfried, zubradt
    #[arg(short = 'm', long = "mutation-method", default_value = "kfactor")]
    pub mutation_method: String,

    // Method-specific parameters (only used with specific methods)
    /// Pseudocount parameter to avoid zero-value issues in logarithmic calculations (Ding method)
    #[arg(long = "pseudocount", default_value_t = 1.0)]
    pub pseudocount: f64,
    /// Maximum score limit parameter to limit the upper bound of reactivity values (Ding method)
    #[arg(long = "maxscore", default_value_t = 10.0)]
    pub maxscore: f64,
    /// SNP threshold for filtering positions (filter positions where unmod sample mutation rate >= cutoff)
    /// This matches the strategy used by RNA Framework (--max-untreated-mut) and ShapeMapper2 (--max-bg)
    /// Unmod samples should have low mutation rates; high rates indicate SNPs or sequencing errors
    #[arg(long = "snp-cutoff", default_value_t = 0.25)]
    pub snp_cutoff: f64,
    /// k-factor prediction method: background (default, from background regions), distribution (based on statistical distribution), recursive (requires --structure-file)
    #[arg(long = "k-prediction-method", default_value = "background")]
    pub k_prediction_method: String,
    /// Reference secondary structure file (required for recursive method, format: Vienna/BPSEQ with sequence and structure)
    #[arg(long = "structure-file")]
    pub structure_file: Option<String>,
    /// Gene ID to use for k-factor calculation (only positions with this gene_id will be used as background regions)
    #[arg(long = "k-background-gene-id")]
    pub k_background_gene_id: Option<String>,

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
    Kfactor, // Kfactor method: difference
    Ding,    // Ding method
    Rouskin, // Rouskin method: use only mod sample stop signal
}

impl std::str::FromStr for StopReactivityMethod {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "kfactor" => Ok(StopReactivityMethod::Kfactor),
            "ding" => Ok(StopReactivityMethod::Ding),
            "rouskin" => Ok(StopReactivityMethod::Rouskin),
            _ => Err(format!("Unknown stop signal reactivity method: {}", s)),
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub enum MutationReactivityMethod {
    Kfactor,   // Kfactor method
    Siegfried, // Siegfried method
    Zubradt,   // Zubradt method: use only mod sample mutation signal
}

impl std::str::FromStr for MutationReactivityMethod {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "kfactor" => Ok(MutationReactivityMethod::Kfactor),
            "siegfried" => Ok(MutationReactivityMethod::Siegfried),
            "zubradt" => Ok(MutationReactivityMethod::Zubradt),
            _ => Err(format!("Unknown mutation signal reactivity method: {}", s)),
        }
    }
}

/// Calculate reactivity with k-factor correction
/// Formula: Rᵢ = Tᵢ - k * Uᵢ
/// where k is a correction factor estimated from the signal intensity ratio in background regions
/// 
/// Improved version: allows negative values (similar to Siegfried method)
/// Negative values indicate that background mutation rate is higher than modified sample
fn calculate_k_corrected_reactivity(mod_signal: f64, unmod_signal: f64, k_factor: f64) -> f64 {
    mod_signal - k_factor * unmod_signal
    // Allow negative values - no clipping
}

/// Distribution-based k-factor prediction targets (based on optimal k analysis from 4 samples)
/// These values are derived from statistical distribution analysis of optimal k values
struct DistributionTargets {
    vs_low_p10: f64,
    vs_low_p25: f64,
    vs_low_p50: f64,
    vs_low_p75: f64,
    vs_low_p90: f64,
    vs_low_iqr: f64,
    vs_low_mad: f64,
    vs_overall_p10: f64,
    vs_overall_p25: f64,
    vs_overall_p50: f64,
    vs_overall_p75: f64,
    vs_overall_p90: f64,
    vs_overall_iqr: f64,
    vs_overall_mad: f64,
    correction_p10: f64,
    correction_p25: f64,
    correction_p50: f64,
    correction_p75: f64,
    correction_p90: f64,
    correction_iqr: f64,
    correction_mad: f64,
}

impl DistributionTargets {
    fn new() -> Self {
        // Values from feature_distribution_stats.csv analysis
        Self {
            vs_low_p10: 7.8414,
            vs_low_p25: 8.9416,
            vs_low_p50: 11.7483,
            vs_low_p75: 14.3412,
            vs_low_p90: 15.0563,
            vs_low_iqr: 5.3996,
            vs_low_mad: 2.9901,
            vs_overall_p10: 4.4704,
            vs_overall_p25: 4.6323,
            vs_overall_p50: 5.0081,
            vs_overall_p75: 5.4478,
            vs_overall_p90: 5.7248,
            vs_overall_iqr: 0.8155,
            vs_overall_mad: 0.4657,
            correction_p10: 0.3297,
            correction_p25: 0.4243,
            correction_p50: 0.6202,
            correction_p75: 0.7732,
            correction_p90: 0.7907,
            correction_iqr: 0.3489,
            correction_mad: 0.1628,
        }
    }
}

/// Calculate features for a given k value
fn calculate_features_for_k(
    mod_rates: &[f64],
    unmod_rates: &[f64],
    k: f64,
    high_mask: &[bool],
    low_mask: &[bool],
) -> Option<(f64, f64, f64)> {
    if mod_rates.len() != unmod_rates.len()
        || mod_rates.len() != high_mask.len()
        || mod_rates.len() != low_mask.len()
    {
        return None;
    }

    let mut high_reactivity = Vec::new();
    let mut low_reactivity = Vec::new();
    let mut overall_reactivity = Vec::new();

    for i in 0..mod_rates.len() {
        let reactivity = (mod_rates[i] - k * unmod_rates[i]).max(0.0);
        overall_reactivity.push(reactivity);
        if high_mask[i] {
            high_reactivity.push(reactivity);
        }
        if low_mask[i] {
            low_reactivity.push(reactivity);
        }
    }

    if high_reactivity.is_empty() || low_reactivity.is_empty() {
        return None;
    }

    let mean_high: f64 = high_reactivity.iter().sum::<f64>() / high_reactivity.len() as f64;
    let mean_low: f64 = low_reactivity.iter().sum::<f64>() / low_reactivity.len() as f64;
    let mean_overall: f64 = overall_reactivity.iter().sum::<f64>() / overall_reactivity.len() as f64;

    if mean_low <= 0.0 || mean_overall <= 0.0 {
        return None;
    }

    let vs_low_ratio = mean_high / mean_low;
    let vs_overall_ratio = mean_high / mean_overall;

    let mean_high_mut_rate: f64 = mod_rates
        .iter()
        .enumerate()
        .filter(|(i, _)| high_mask[*i])
        .map(|(_, &rate)| rate)
        .sum::<f64>()
        / high_reactivity.len() as f64;
    let correction_ratio = if mean_high_mut_rate > 0.0 {
        mean_high / mean_high_mut_rate
    } else {
        0.0
    };

    Some((vs_low_ratio, vs_overall_ratio, correction_ratio))
}

/// Calculate loss for a single feature value
fn feature_loss(value: f64, target: &DistributionTargets, feature_type: &str) -> f64 {
    let (p10, p25, p75, p90, iqr, mad) = match feature_type {
        "vs_low" => (
            target.vs_low_p10,
            target.vs_low_p25,
            target.vs_low_p75,
            target.vs_low_p90,
            target.vs_low_iqr,
            target.vs_low_mad,
        ),
        "vs_overall" => (
            target.vs_overall_p10,
            target.vs_overall_p25,
            target.vs_overall_p75,
            target.vs_overall_p90,
            target.vs_overall_iqr,
            target.vs_overall_mad,
        ),
        "correction" => (
            target.correction_p10,
            target.correction_p25,
            target.correction_p75,
            target.correction_p90,
            target.correction_iqr,
            target.correction_mad,
        ),
        _ => return f64::INFINITY,
    };

    // If value is in IQR range (P25-P75), loss is 0
    if value >= p25 && value <= p75 {
        return 0.0;
    }

    // If value is in P10-P90 range, use IQR-normalized distance
    if value >= p10 && value <= p90 {
        let dist = if value < p25 {
            (p25 - value) / iqr
        } else {
            (value - p75) / iqr
        };
        return dist * 0.5; // Smaller penalty
    }

    // Beyond P10-P90 range, use MAD-normalized distance
    let dist = if value < p10 {
        (p10 - value) / mad
    } else {
        (value - p90) / mad
    };
    dist * 2.0 // Larger penalty
}

/// Calculate distribution-based loss for a given k value
fn calculate_distribution_loss(
    k: f64,
    mod_rates: &[f64],
    unmod_rates: &[f64],
    high_mask: &[bool],
    low_mask: &[bool],
    target: &DistributionTargets,
) -> f64 {
    if let Some((vs_low_ratio, vs_overall_ratio, correction_ratio)) =
        calculate_features_for_k(mod_rates, unmod_rates, k, high_mask, low_mask)
    {
        let loss_vs_low = feature_loss(vs_low_ratio, target, "vs_low");
        let loss_vs_overall = feature_loss(vs_overall_ratio, target, "vs_overall");
        let loss_correction = feature_loss(correction_ratio, target, "correction");

        // Weighted combination
        2.5 * loss_vs_low + 2.0 * loss_vs_overall + 1.0 * loss_correction
    } else {
        f64::INFINITY
    }
}

/// Predict k value using distribution-based method
/// Uses percentile thresholds (top 10% for high, 40-60% for low) and loss minimization
fn predict_k_by_distribution(
    mod_rates: &[f64],
    unmod_rates: &[f64],
    k_range: (f64, f64),
    step: f64,
) -> Option<f64> {
    if mod_rates.is_empty() || mod_rates.len() != unmod_rates.len() {
        return None;
    }

    // Use percentile thresholds: top 10% for high, 40-60% for low
    let mut sorted_rates: Vec<(usize, f64)> = mod_rates
        .iter()
        .enumerate()
        .map(|(i, &rate)| (i, rate))
        .collect();
    sorted_rates.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));

    let high_percentile = 90;
    let low_percentile_low = 40;
    let low_percentile_high = 60;

    let high_threshold_idx = (sorted_rates.len() * high_percentile / 100).min(sorted_rates.len() - 1);
    let low_threshold_low_idx = (sorted_rates.len() * low_percentile_low / 100).min(sorted_rates.len() - 1);
    let low_threshold_high_idx = (sorted_rates.len() * low_percentile_high / 100).min(sorted_rates.len() - 1);

    let high_threshold = sorted_rates[high_threshold_idx].1;
    let low_threshold_low = sorted_rates[low_threshold_low_idx].1;
    let low_threshold_high = sorted_rates[low_threshold_high_idx].1;

    let mut high_mask = vec![false; mod_rates.len()];
    let mut low_mask = vec![false; mod_rates.len()];

    for (i, &rate) in mod_rates.iter().enumerate() {
        high_mask[i] = rate >= high_threshold;
        low_mask[i] = rate >= low_threshold_low && rate < low_threshold_high;
    }

    if high_mask.iter().filter(|&&x| x).count() < 10
        || low_mask.iter().filter(|&&x| x).count() < 10
    {
        // Try more relaxed thresholds
        let high_percentile = 85;
        let low_percentile_low = 30;
        let low_percentile_high = 50;

        let high_threshold_idx = (sorted_rates.len() * high_percentile / 100).min(sorted_rates.len() - 1);
        let low_threshold_low_idx = (sorted_rates.len() * low_percentile_low / 100).min(sorted_rates.len() - 1);
        let low_threshold_high_idx = (sorted_rates.len() * low_percentile_high / 100).min(sorted_rates.len() - 1);

        let high_threshold = sorted_rates[high_threshold_idx].1;
        let low_threshold_low = sorted_rates[low_threshold_low_idx].1;
        let low_threshold_high = sorted_rates[low_threshold_high_idx].1;

        for (i, &rate) in mod_rates.iter().enumerate() {
            high_mask[i] = rate >= high_threshold;
            low_mask[i] = rate >= low_threshold_low && rate < low_threshold_high;
        }

        if high_mask.iter().filter(|&&x| x).count() < 10
            || low_mask.iter().filter(|&&x| x).count() < 10
        {
            return None;
        }
    }

    let target = DistributionTargets::new();
    let mut best_k = None;
    let mut best_loss = f64::INFINITY;

    let mut k = k_range.0;
    while k <= k_range.1 {
        let loss = calculate_distribution_loss(k, mod_rates, unmod_rates, &high_mask, &low_mask, &target);
        if loss < best_loss {
            best_loss = loss;
            best_k = Some(k);
        }
        k += step;
    }

    best_k
}

/// Metrics for k value evaluation
#[derive(Debug, Clone)]
struct KValueMetrics {
    auc: f64,
    sensitivity: f64,
    specificity: f64,
    negative_ratio: f64, // Ratio of negative reactivity values
    zero_count_double_stranded: usize, // Count of zero values in double-stranded regions
    zero_count_single_stranded: usize, // Count of zero values in single-stranded regions
    zero_ratio_double_stranded: f64, // Ratio of zero values in double-stranded regions
    zero_ratio_single_stranded: f64, // Ratio of zero values in single-stranded regions
}

/// Calculate metrics (AUC, sensitivity, specificity) for a given k value using reactivity and structure information
/// Returns metrics or None if calculation fails
/// Allows negative reactivity values (similar to Siegfried method)
fn calculate_metrics_for_k(
    k: f64,
    mod_rates: &[f64],
    unmod_rates: &[f64],
    is_single_stranded: &[bool],
) -> Option<KValueMetrics> {
    if mod_rates.len() != unmod_rates.len() || mod_rates.len() != is_single_stranded.len() {
        return None;
    }

    // Calculate reactivity for each position (allow negative values)
    let mut reactivity_and_label: Vec<(f64, bool)> = mod_rates
        .iter()
        .zip(unmod_rates.iter())
        .zip(is_single_stranded.iter())
        .map(|((&m_rate, &u_rate), &is_single)| {
            let reactivity = m_rate - k * u_rate; // Allow negative values
            (reactivity, is_single)
        })
        .filter(|(r, _)| !r.is_nan())
        .collect();

    if reactivity_and_label.is_empty() {
        return None;
    }

    // Calculate negative ratio (proportion of negative reactivity values)
    let negative_count = reactivity_and_label.iter().filter(|(r, _)| *r < 0.0).count();
    let negative_ratio = negative_count as f64 / reactivity_and_label.len() as f64;
    
    // Calculate zero value statistics by structure type
    // Zero values: reactivity <= threshold (considering small positive values near zero as effectively zero)
    let zero_threshold = 1e-6;
    let zero_count_double_stranded = reactivity_and_label.iter()
        .filter(|(r, is_single)| *r <= zero_threshold && !*is_single)
        .count();
    let zero_count_single_stranded = reactivity_and_label.iter()
        .filter(|(r, is_single)| *r <= zero_threshold && *is_single)
        .count();
    
    // Calculate zero ratios
    let total_double_stranded = reactivity_and_label.iter().filter(|(_, is_single)| !*is_single).count();
    let total_single_stranded = reactivity_and_label.iter().filter(|(_, is_single)| *is_single).count();
    
    let zero_ratio_double_stranded = if total_double_stranded > 0 {
        zero_count_double_stranded as f64 / total_double_stranded as f64
    } else {
        0.0
    };
    
    let zero_ratio_single_stranded = if total_single_stranded > 0 {
        zero_count_single_stranded as f64 / total_single_stranded as f64
    } else {
        0.0
    };

    // Sort by reactivity (descending)
    reactivity_and_label.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap_or(std::cmp::Ordering::Equal));

    // Calculate total positive and negative counts
    let total_positive = reactivity_and_label.iter().filter(|(_, is_single)| *is_single).count();
    let total_negative = reactivity_and_label.len() - total_positive;

    if total_positive == 0 || total_negative == 0 {
        return None;
    }

    // Find optimal threshold (Youden's J statistic: maximizes sensitivity + specificity - 1)
    let mut best_threshold = 0.0;
    let mut best_sensitivity = 0.0;
    let mut best_specificity = 0.0;
    let mut best_youden = -1.0;

    // Try different thresholds (use reactivity values as thresholds)
    let mut thresholds: Vec<f64> = reactivity_and_label.iter().map(|(r, _)| *r).collect();
    thresholds.sort_by(|a, b| b.partial_cmp(a).unwrap_or(std::cmp::Ordering::Equal));
    thresholds.dedup();

    for threshold in thresholds {
        let mut tp = 0;
        let mut fp = 0;
        let mut tn = 0;
        let mut fn_count = 0;

        for (reactivity, is_single) in &reactivity_and_label {
            let predicted_positive = *reactivity >= threshold;
            if *is_single {
                if predicted_positive {
                    tp += 1;
                } else {
                    fn_count += 1;
                }
            } else {
                if predicted_positive {
                    fp += 1;
                } else {
                    tn += 1;
                }
            }
        }

        let sensitivity = if tp + fn_count > 0 { tp as f64 / (tp + fn_count) as f64 } else { 0.0 };
        let specificity = if tn + fp > 0 { tn as f64 / (tn + fp) as f64 } else { 0.0 };
        let youden = sensitivity + specificity - 1.0;

        if youden > best_youden {
            best_youden = youden;
            best_threshold = threshold;
            best_sensitivity = sensitivity;
            best_specificity = specificity;
        }
    }

    // Calculate AUC using ROC curve
    let mut auc = 0.0;
    let mut tp = 0;
    let mut fp = 0;
    let mut prev_fpr = 0.0;
    let mut prev_tpr = 0.0;

    for (_, is_single) in &reactivity_and_label {
        if *is_single {
            tp += 1;
        } else {
            fp += 1;
        }

        let tpr = tp as f64 / total_positive as f64;
        let fpr = fp as f64 / total_negative as f64;

        // Calculate AUC using trapezoidal rule
        auc += (fpr - prev_fpr) * (prev_tpr + tpr) / 2.0;

        prev_fpr = fpr;
        prev_tpr = tpr;
    }

    // Add final point
    auc += (1.0 - prev_fpr) * (prev_tpr + 1.0) / 2.0;

    Some(KValueMetrics {
        auc,
        sensitivity: best_sensitivity,
        specificity: best_specificity,
        negative_ratio,
        zero_count_double_stranded,
        zero_count_single_stranded,
        zero_ratio_double_stranded,
        zero_ratio_single_stranded,
    })
}

/// Calculate AUC for a given k value (backward compatibility)
/// Returns AUC or None if calculation fails
fn calculate_auc_for_k(
    k: f64,
    mod_rates: &[f64],
    unmod_rates: &[f64],
    is_single_stranded: &[bool],
) -> Option<f64> {
    calculate_metrics_for_k(k, mod_rates, unmod_rates, is_single_stranded)
        .map(|metrics| metrics.auc)
}

/// Calculate composite score for k value selection
/// Balances AUC, sensitivity, and specificity, with penalty for high negative ratio
/// Also considers zero value distribution: penalizes when single-stranded zeros increase while double-stranded zeros don't
/// Formula: score = w_auc * AUC + w_sens * sensitivity + w_spec * specificity - penalty_neg - penalty_zero
/// where penalty_neg increases exponentially with negative_ratio
/// and penalty_zero penalizes inappropriate zero value distribution
fn calculate_composite_score(metrics: &KValueMetrics) -> f64 {
    let w_auc = 0.4;      // Weight for AUC
    let w_sens = 0.35;     // Weight for sensitivity (higher weight to avoid false negatives)
    let w_spec = 0.25;     // Weight for specificity
    
    // Penalty for negative ratio: exponential penalty to strongly discourage excessive false negatives
    // When negative_ratio > 0.5, penalty becomes very large
    let penalty_neg = if metrics.negative_ratio > 0.5 {
        // Strong penalty for high negative ratio (>50%)
        2.0 * (metrics.negative_ratio - 0.5).powi(2)
    } else if metrics.negative_ratio > 0.3 {
        // Moderate penalty for moderate negative ratio (30-50%)
        0.5 * (metrics.negative_ratio - 0.3) / 0.2
    } else {
        // Small penalty for low negative ratio (<30%)
        0.1 * metrics.negative_ratio / 0.3
    };
    
    // Penalty for zero value distribution
    // Penalize when single-stranded zero ratio is high while double-stranded zero ratio is not increasing
    // This indicates that k is too high, causing loss of signal in single-stranded regions
    let penalty_zero = if metrics.zero_ratio_single_stranded > 0.3 {
        // If single-stranded zero ratio > 30%, apply penalty
        // Penalty increases with the difference between single and double stranded zero ratios
        let zero_diff = metrics.zero_ratio_single_stranded - metrics.zero_ratio_double_stranded;
        if zero_diff > 0.1 {
            // Strong penalty when single-stranded zeros are much higher than double-stranded zeros
            1.0 * zero_diff
        } else {
            // Moderate penalty
            0.5 * metrics.zero_ratio_single_stranded
        }
    } else {
        0.0
    };
    
    w_auc * metrics.auc + w_sens * metrics.sensitivity + w_spec * metrics.specificity - penalty_neg - penalty_zero
}

/// Predict k value using improved recursive search strategy based on structure
/// Strategy:
/// 1. Start from k=1.0, explore ±3 steps of 0.1 (0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3)
/// 2. Determine search direction by comparing left vs right performance
/// 3. Continue coarse search (0.1 step) in chosen direction until AUC drops for 3 consecutive steps
/// 4. Find best coarse position, then perform fine search (0.01 step) around it
/// 5. Fine search: evaluate 18 new positions + 1 existing best = 19 positions total
/// Uses AUC as evaluation metric (requires structure information)
/// Allows negative reactivity values
fn predict_k_by_recursive_with_structure(
    mod_rates: &[f64],
    unmod_rates: &[f64],
    positions: &[u32],
    is_single_stranded: &[bool],
    _k_range: (f64, f64), // Not used in new strategy, kept for compatibility
    _coarse_step: f64,    // Not used in new strategy, kept for compatibility
    _fine_step: f64,      // Not used in new strategy, kept for compatibility
    _fine_range: f64,     // Not used in new strategy, kept for compatibility
) -> Option<f64> {
    if mod_rates.is_empty() 
        || mod_rates.len() != unmod_rates.len() 
        || mod_rates.len() != positions.len()
        || mod_rates.len() != is_single_stranded.len() {
        return None;
    }

    // Step 1: Initial exploration around k=1.0
    // Explore ±3 steps of 0.1: 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3
    let initial_k_values: Vec<f64> = (7..=13).map(|i| i as f64 / 10.0).collect();
    
    // Parallel evaluation of initial k values with comprehensive metrics
    let initial_results: Vec<(f64, Option<KValueMetrics>)> = initial_k_values
        .into_par_iter()
        .map(|k| {
            let metrics = calculate_metrics_for_k(k, mod_rates, unmod_rates, is_single_stranded);
            (k, metrics)
        })
        .collect();

    // Find best initial k using composite score
    let initial_scores: Vec<(f64, f64, KValueMetrics)> = initial_results
        .iter()
        .filter_map(|(k, metrics_opt)| {
            metrics_opt.as_ref().map(|metrics| {
                let score = calculate_composite_score(metrics);
                (*k, score, metrics.clone())
            })
        })
        .collect();

    if initial_scores.is_empty() {
        return None;
    }

    // Find best initial position based on composite score
    let (best_initial_k, best_initial_score, best_initial_metrics) = initial_scores
        .iter()
        .max_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal))
        .map(|(k, score, metrics)| (*k, *score, metrics.clone()))?;
    
    let best_initial_auc = best_initial_metrics.auc;

    // Determine search direction: consider both average composite score and best score location
    let left_scores: Vec<f64> = initial_scores
        .iter()
        .filter(|(k, _, _)| *k < 1.0)
        .map(|(_, score, _)| *score)
        .collect();
    
    let right_scores: Vec<f64> = initial_scores
        .iter()
        .filter(|(k, _, _)| *k > 1.0)
        .map(|(_, score, _)| *score)
        .collect();

    let left_avg_score: f64 = if !left_scores.is_empty() {
        left_scores.iter().sum::<f64>() / left_scores.len() as f64
    } else {
        -1.0
    };
    
    let right_avg_score: f64 = if !right_scores.is_empty() {
        right_scores.iter().sum::<f64>() / right_scores.len() as f64
    } else {
        -1.0
    };

    // Check which side has the best composite score
    let best_score_on_right = best_initial_k > 1.0;
    let best_score_on_left = best_initial_k < 1.0;
    let best_score_at_center = (best_initial_k - 1.0).abs() < 1e-6;
    
    // Determine search strategy: consider both best score location and average score
    // Strategy 1: If best score and average score are consistent (both on same side), search that side only
    // Strategy 2: If inconsistent (best score on one side, but average score higher on other side), search both sides
    let (search_right_only, search_left_only, search_both) = if best_score_at_center {
        // Best score is at k=1.0, use average score to decide single direction
        if right_avg_score > left_avg_score {
            (true, false, false)
        } else {
            (false, true, false)
        }
    } else if best_score_on_right {
        // Best score is on right side
        if right_avg_score >= left_avg_score {
            // Consistent: best score and average score both favor right
            (true, false, false)
        } else {
            // Inconsistent: best score on right, but average score higher on left
            (false, false, true)
        }
    } else {
        // Best score is on left side
        if left_avg_score >= right_avg_score {
            // Consistent: best score and average score both favor left
            (false, true, false)
        } else {
            // Inconsistent: best score on left, but average score higher on right
            (false, false, true)
        }
    };

    // Step 2: Coarse search based on strategy
    let mut best_coarse_k = best_initial_k;
    let mut best_coarse_score = best_initial_score;
    let mut best_coarse_auc = best_initial_auc;
    let coarse_step = 0.1;
    let max_iterations = 20; // Safety limit

    // Helper function to search in one direction using composite score
    // Also tracks zero value distribution to stop when single-stranded zeros increase while double-stranded zeros don't
    let search_in_direction = |start_k: f64, start_score: f64, direction_right: bool| -> (f64, f64, f64) {
        let mut current_k = start_k;
        let mut consecutive_drops = 0;
        let mut best_k = start_k;
        let mut best_score = start_score;
        let mut best_auc = best_initial_auc;
        let mut last_score = start_score;
        
        // Track zero value statistics from previous iteration
        let mut prev_zero_count_double = 0;
        let mut prev_zero_count_single = 0;
        
        // Get initial zero statistics
        if let Some(initial_metrics) = calculate_metrics_for_k(start_k, mod_rates, unmod_rates, is_single_stranded) {
            prev_zero_count_double = initial_metrics.zero_count_double_stranded;
            prev_zero_count_single = initial_metrics.zero_count_single_stranded;
        }

        for _ in 0..max_iterations {
            // Move in chosen direction
            if direction_right {
                current_k += coarse_step;
            } else {
                current_k -= coarse_step;
            }

            // Safety bounds
            if current_k < 0.1 || current_k > 3.0 {
                break;
            }

            if let Some(metrics) = calculate_metrics_for_k(current_k, mod_rates, unmod_rates, is_single_stranded) {
                // Check zero value distribution constraint
                // Stop if: double-stranded zeros don't increase AND single-stranded zeros increase
                // This indicates k is too high, causing loss of signal in single-stranded regions
                let zero_double_increased = metrics.zero_count_double_stranded > prev_zero_count_double;
                let zero_single_increased = metrics.zero_count_single_stranded > prev_zero_count_single;
                
                // Stop condition: double-stranded zeros not increasing while single-stranded zeros are increasing
                if !zero_double_increased && zero_single_increased && prev_zero_count_double > 0 {
                    // This indicates k is too high, stop searching in this direction
                    break;
                }
                
                // Update previous zero counts for next iteration
                prev_zero_count_double = metrics.zero_count_double_stranded;
                prev_zero_count_single = metrics.zero_count_single_stranded;
                
                let score = calculate_composite_score(&metrics);
                if score > best_score {
                    best_score = score;
                    best_k = current_k;
                    best_auc = metrics.auc;
                    consecutive_drops = 0;
                } else if score < last_score {
                    consecutive_drops += 1;
                    if consecutive_drops >= 3 {
                        break; // Stop when score drops for 3 consecutive steps
                    }
                } else {
                    consecutive_drops = 0; // Reset if score doesn't drop
                }
                
                last_score = score;
            } else {
                consecutive_drops += 1;
                if consecutive_drops >= 3 {
                    break;
                }
            }
        }
        (best_k, best_score, best_auc)
    };

    // Execute search based on strategy
    if search_right_only {
        // Search only right side
        let (k, score, auc) = search_in_direction(best_initial_k, best_initial_score, true);
        if score > best_coarse_score {
            best_coarse_k = k;
            best_coarse_score = score;
            best_coarse_auc = auc;
        }
    } else if search_left_only {
        // Search only left side
        let (k, score, auc) = search_in_direction(best_initial_k, best_initial_score, false);
        if score > best_coarse_score {
            best_coarse_k = k;
            best_coarse_score = score;
            best_coarse_auc = auc;
        }
    } else if search_both {
        // Search both sides (inconsistent case)
        let (k_right, score_right, auc_right) = search_in_direction(best_initial_k, best_initial_score, true);
        let (k_left, score_left, auc_left) = search_in_direction(best_initial_k, best_initial_score, false);
        
        // Choose the best result from both sides based on composite score
        if score_right > score_left && score_right > best_coarse_score {
            best_coarse_k = k_right;
            best_coarse_score = score_right;
            best_coarse_auc = auc_right;
        } else if score_left > best_coarse_score {
            best_coarse_k = k_left;
            best_coarse_score = score_left;
            best_coarse_auc = auc_left;
        }
    }

    // Step 3: Fine search around best coarse position
    // Fine search: 0.01 step, evaluate 18 new positions + 1 existing best = 19 positions total
    // Strategy: Generate 9 positions on left side, 9 positions on right side, plus best_coarse_k
    let fine_step = 0.01;
    
    let mut fine_k_values = Vec::new();
    
    // Generate 9 positions on the left side (excluding best_coarse_k)
    for i in 1..=9 {
        let k = best_coarse_k - (i as f64 * fine_step);
        if k >= 0.1 {
            let k_rounded = (k * 100.0).round() / 100.0;
            fine_k_values.push(k_rounded);
        }
    }
    
    // Add best_coarse_k itself (the +1 existing result)
    let best_coarse_rounded = (best_coarse_k * 100.0).round() / 100.0;
    fine_k_values.push(best_coarse_rounded);
    
    // Generate 9 positions on the right side (excluding best_coarse_k)
    for i in 1..=9 {
        let k = best_coarse_k + (i as f64 * fine_step);
        if k <= 3.0 {
            let k_rounded = (k * 100.0).round() / 100.0;
            fine_k_values.push(k_rounded);
        }
    }
    
    // Sort to ensure consistent order
    fine_k_values.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

    // Parallel fine search using composite score
    let fine_results: Vec<(f64, Option<KValueMetrics>)> = fine_k_values
        .into_par_iter()
        .map(|k| {
            let metrics = calculate_metrics_for_k(k, mod_rates, unmod_rates, is_single_stranded);
            (k, metrics)
        })
        .collect();

    // Find best k from fine search (highest composite score)
    let best_fine = fine_results
        .iter()
        .filter_map(|(k, metrics_opt)| {
            metrics_opt.as_ref().map(|metrics| {
                let score = calculate_composite_score(metrics);
                (*k, score)
            })
        })
        .max_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));

    match best_fine {
        Some((k, _)) => Some(k),
        None => Some(best_coarse_k), // Fallback to coarse result
    }
}

/// Predict k value using recursive method with parallel two-phase search
/// Phase 1: Coarse search with large step size to locate optimal region
/// Phase 2: Fine search with small step size in the optimal region
/// Uses distribution-based loss as evaluation metric
/// Note: Uses rayon's global thread pool for parallelization
fn predict_k_by_recursive(
    mod_rates: &[f64],
    unmod_rates: &[f64],
    k_range: (f64, f64),
    coarse_step: f64,
    fine_step: f64,
    fine_range: f64,
    _thread_count: usize, // Reserved for future use, rayon uses global thread pool
) -> Option<f64> {
    if mod_rates.is_empty() || mod_rates.len() != unmod_rates.len() {
        return None;
    }

    // Prepare masks for distribution-based evaluation (same as distribution method)
    let mut sorted_rates: Vec<(usize, f64)> = mod_rates
        .iter()
        .enumerate()
        .map(|(i, &rate)| (i, rate))
        .collect();
    sorted_rates.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));

    let high_percentile = 90;
    let low_percentile_low = 40;
    let low_percentile_high = 60;

    let high_threshold_idx = (sorted_rates.len() * high_percentile / 100).min(sorted_rates.len() - 1);
    let low_threshold_low_idx = (sorted_rates.len() * low_percentile_low / 100).min(sorted_rates.len() - 1);
    let low_threshold_high_idx = (sorted_rates.len() * low_percentile_high / 100).min(sorted_rates.len() - 1);

    let high_threshold = sorted_rates[high_threshold_idx].1;
    let low_threshold_low = sorted_rates[low_threshold_low_idx].1;
    let low_threshold_high = sorted_rates[low_threshold_high_idx].1;

    let mut high_mask = vec![false; mod_rates.len()];
    let mut low_mask = vec![false; mod_rates.len()];

    for (i, &rate) in mod_rates.iter().enumerate() {
        high_mask[i] = rate >= high_threshold;
        low_mask[i] = rate >= low_threshold_low && rate < low_threshold_high;
    }

    if high_mask.iter().filter(|&&x| x).count() < 10
        || low_mask.iter().filter(|&&x| x).count() < 10
    {
        // Try more relaxed thresholds
        let high_percentile = 85;
        let low_percentile_low = 30;
        let low_percentile_high = 50;

        let high_threshold_idx = (sorted_rates.len() * high_percentile / 100).min(sorted_rates.len() - 1);
        let low_threshold_low_idx = (sorted_rates.len() * low_percentile_low / 100).min(sorted_rates.len() - 1);
        let low_threshold_high_idx = (sorted_rates.len() * low_percentile_high / 100).min(sorted_rates.len() - 1);

        let high_threshold = sorted_rates[high_threshold_idx].1;
        let low_threshold_low = sorted_rates[low_threshold_low_idx].1;
        let low_threshold_high = sorted_rates[low_threshold_high_idx].1;

        for (i, &rate) in mod_rates.iter().enumerate() {
            high_mask[i] = rate >= high_threshold;
            low_mask[i] = rate >= low_threshold_low && rate < low_threshold_high;
        }

        if high_mask.iter().filter(|&&x| x).count() < 10
            || low_mask.iter().filter(|&&x| x).count() < 10
        {
            return None;
        }
    }

    let target = DistributionTargets::new();
    
    // Phase 1: Coarse search (parallel)
    let coarse_k_start = k_range.0;
    let coarse_k_end = k_range.1;
    let mut coarse_k_values = Vec::new();
    let mut k = coarse_k_start;
    while k <= coarse_k_end {
        coarse_k_values.push(k);
        k += coarse_step;
    }

    // Parallel coarse search
    let coarse_results: Vec<(f64, f64)> = coarse_k_values
        .into_par_iter()
        .map(|k| {
            let loss = calculate_distribution_loss(k, mod_rates, unmod_rates, &high_mask, &low_mask, &target);
            (k, loss)
        })
        .collect();

    // Find best k from coarse search
    let best_coarse = coarse_results
        .iter()
        .min_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));

    let best_coarse_k = match best_coarse {
        Some((k, _)) => *k,
        None => return None,
    };

    // Phase 2: Fine search (parallel) around best_coarse_k
    let fine_start = (best_coarse_k - fine_range).max(k_range.0);
    let fine_end = (best_coarse_k + fine_range).min(k_range.1);
    
    let mut fine_k_values = Vec::new();
    let mut k = fine_start;
    while k <= fine_end {
        fine_k_values.push(k);
        k += fine_step;
    }

    // Parallel fine search
    let fine_results: Vec<(f64, f64)> = fine_k_values
        .into_par_iter()
        .map(|k| {
            let loss = calculate_distribution_loss(k, mod_rates, unmod_rates, &high_mask, &low_mask, &target);
            (k, loss)
        })
        .collect();

    // Find best k from fine search
    let best_fine = fine_results
        .iter()
        .min_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));

    match best_fine {
        Some((k, _)) => Some(*k),
        None => Some(best_coarse_k), // Fallback to coarse result
    }
}

/// Ding et al., 2014 method reactivity calculation
/// Correct Ding method: directly calculate difference after log normalization
/// Formula: S_i = max(0, T_i - U_i)
/// where T_i = ln(n_Ti + p) / mean(ln(n_Tj + p)), U_i = ln(n_Ui + p) / mean(ln(n_Uj + p))
fn calculate_reactivity_ding_style(
    mod_signal: f64,
    unmod_signal: f64,
    _k_factor: f64, // Not using k-factor correction
    pseudocount: f64,
    _maxscore: f64, // Not using maxscore limit
) -> f64 {
    // Apply log transformation to raw signals (add pseudocount to avoid ln(0))
    let log_mod = (mod_signal + pseudocount).ln();
    let log_unmod = (unmod_signal + pseudocount).ln();

    // Directly calculate difference (not using k-factor correction)
    let reactivity = log_mod - log_unmod;

    // Return max with 0
    reactivity.max(0.0)
}

fn parse_signal(
    fields: &Vec<String>,
    rfc_idx: usize,
    pipc_idx: usize,
    depth_idx: usize,
    effective_depth_idx: Option<usize>,
) -> (f64, f64) {
    // Use effective_depth if available, otherwise fall back to depth
    let depth: f64 = if let Some(eff_idx) = effective_depth_idx {
        // Try to parse effective_depth, fall back to depth if parsing fails
        fields.get(eff_idx)
            .and_then(|s| s.parse::<f64>().ok())
            .unwrap_or_else(|| fields.get(depth_idx).and_then(|s| s.parse::<f64>().ok()).unwrap_or(0.0))
    } else {
        fields.get(depth_idx).and_then(|s| s.parse::<f64>().ok()).unwrap_or(0.0)
    };
    let rfc: f64 = fields.get(rfc_idx).and_then(|s| s.parse::<f64>().ok()).unwrap_or(0.0);
    let pipc: f64 = fields.get(pipc_idx).and_then(|s| s.parse::<f64>().ok()).unwrap_or(0.0);
    let mutation = if depth > 0.0 { rfc / depth } else { 0.0 };
    let stop = if depth > 0.0 { pipc / depth } else { 0.0 };
    (mutation, stop)
}

pub fn calc_reactivity(
    args: &ReactivityArgs,
    logger: &mut crate::Logger,
) -> Result<(), Box<dyn Error>> {
    // Validate reactivity command parameters
    validate_reactivity_args(args)?;
    let start_time = Instant::now();
    let has_unmod = args
        .unmod_csv
        .as_ref()
        .map(|p| !p.trim().is_empty())
        .unwrap_or(false);

    // Record environment information and parameters
    logger.log("=== ModDetector Reactivity Function Log ===")?;
    logger.log(&format!("Software Version: v{}", crate::VERSION))?;
    logger.log(&format!(
        "Runtime: {}",
        chrono::Utc::now().format("%Y-%m-%d %H:%M:%S")
    ))?;
    logger.log(&format!("Modified Sample File: {}", args.mod_csv))?;
    if let Some(unmod_path) = args.unmod_csv.as_ref() {
        logger.log(&format!("Unmodified Sample File: {}", unmod_path))?;
    } else {
        logger.log("Unmodified Sample File: <not provided> (mod-only mode)")?;
    }
    logger.log(&format!("Reactivity Output: {}", args.output))?;
    logger.log(&format!("Stop Method: {}", args.stop_method))?;
    logger.log(&format!("Mutation Method: {}", args.mutation_method))?;
    logger.log(&format!("Pseudocount Parameter: {}", args.pseudocount))?;
    logger.log(&format!("Max Score Limit: {}", args.maxscore))?;
    logger.log(&format!("SNP Cutoff: {} (filters positions where unmod sample mutation rate >= cutoff)", args.snp_cutoff))?;
    logger.log(&format!("K Prediction Method: {}", args.k_prediction_method))?;
    if let Some(threads) = args.threads {
        logger.log(&format!("Threads: {}", threads))?;
    } else {
        logger.log("Threads: Auto-detect")?;
    }
    logger.log("Starting reactivity calculation...")?;

    // Parse stop signal method
    let stop_method: StopReactivityMethod = args
        .stop_method
        .parse()
        .unwrap_or(StopReactivityMethod::Kfactor);

    // Parse mutation signal method
    let mutation_method: MutationReactivityMethod = args
        .mutation_method
        .parse()
        .unwrap_or(MutationReactivityMethod::Kfactor);

    // Set thread count
    let thread_count = args.threads.unwrap_or(8);
    rayon::ThreadPoolBuilder::new()
        .num_threads(thread_count)
        .build_global()
        .ok();

    println!("[Loading data] Threads: {}", thread_count);
    println!("    Mod={}", args.mod_csv);
    match args.unmod_csv.as_ref() {
        Some(unmod_path) => println!("    Unmod={}", unmod_path),
        None => println!("    Unmod=<not provided> (mod-only mode)"),
    }
    if !has_unmod {
        println!("    Mode: mod-only reactivity (no unmodified subtraction)");
    }
    println!();

    // Use chunked parallel data loading
    // Check if effective_depth column exists in the CSV header
    let mod_file = File::open(&args.mod_csv)?;
    let mut mod_header_reader = BufReader::new(mod_file);
    let mut header_line = String::new();
    mod_header_reader.read_line(&mut header_line)?;
    let header_fields: Vec<&str> = header_line.trim().split(',').collect();
    let effective_depth_idx = header_fields.iter().position(|&s| s == "effective_depth");
    if effective_depth_idx.is_some() {
        logger.log("Detected effective_depth column in mod CSV, will use it for mutation rate calculation")?;
        logger.log(&format!("effective_depth column index: {}", effective_depth_idx.unwrap()))?;
    } else {
        logger.log("No effective_depth column found in mod CSV, will use depth column for mutation rate calculation")?;
    }
    
    // Also check unmod CSV header if provided
    if has_unmod {
        let unmod_csv = args.unmod_csv.as_ref().unwrap();
        let unmod_file = File::open(unmod_csv)?;
        let mut unmod_header_reader = BufReader::new(unmod_file);
        let mut unmod_header_line = String::new();
        unmod_header_reader.read_line(&mut unmod_header_line)?;
        let unmod_header_fields: Vec<&str> = unmod_header_line.trim().split(',').collect();
        let unmod_effective_depth_idx = unmod_header_fields.iter().position(|&s| s == "effective_depth");
        if unmod_effective_depth_idx.is_some() {
            logger.log("Detected effective_depth column in unmod CSV, will use it for mutation rate calculation")?;
            logger.log(&format!("unmod effective_depth column index: {}", unmod_effective_depth_idx.unwrap()))?;
        } else {
            logger.log("No effective_depth column found in unmod CSV, will use depth column for mutation rate calculation")?;
        }
        // Use the effective_depth_idx from mod CSV (they should be the same)
        // If unmod doesn't have it, parse_signal will fall back to depth automatically
    }
    
    let (mod_records, unmod_records, unmod_key_map) = if has_unmod {
        let unmod_csv = args.unmod_csv.as_ref().unwrap();
        let (mod_records, unmod_records, _mod_key_map, unmod_key_map, _header) =
            load_data_in_chunks(&args.mod_csv, unmod_csv, thread_count)?;
        (mod_records, unmod_records, unmod_key_map)
    } else {
        println!("    Running in mod-only mode (no unmodified sample provided)");
        let mod_records = load_mod_data_only(&args.mod_csv, thread_count)?;
        (mod_records, Vec::new(), HashMap::new())
    };

    // 1. Collect background region signals in parallel, calculate k-factor
    let mut stop_mod_bg = Vec::new();
    let mut stop_unmod_bg = Vec::new();
    let mut mut_mod_bg = Vec::new();
    let mut mut_unmod_bg = Vec::new();
    let mut all_pairs: Vec<(usize, f64, f64, f64, f64)> = Vec::new(); // (mod_idx, m_mut, m_stop, u_mut, u_stop)

    // Extract SNP cutoff parameter for closure
    let snp_cutoff = args.snp_cutoff;

    let total_mod = mod_records.len();
    let mut progress = crate::progress::SimpleProgress::new(total_mod);

    // Display processing progress
    print!("\r[Processing] ");
    std::io::stdout().flush()?;

    // Process data collection in parallel
    let chunk_size = (total_mod + thread_count - 1) / thread_count;
    let chunks: Vec<Vec<usize>> = (0..total_mod)
        .collect::<Vec<usize>>()
        .chunks(chunk_size)
        .map(|chunk| chunk.to_vec())
        .collect();

    // Capture effective_depth_idx for use in closure
    let effective_depth_idx_clone = effective_depth_idx;

    let results: Vec<(
        Vec<f64>,
        Vec<f64>,
        Vec<f64>,
        Vec<f64>,
        Vec<(usize, f64, f64, f64, f64)>,
    )> = chunks
        .into_par_iter()
        .map(|chunk| {
            let effective_depth_idx = effective_depth_idx_clone;
            let mut local_stop_mod_bg = Vec::new();
            let mut local_stop_unmod_bg = Vec::new();
            let mut local_mut_mod_bg = Vec::new();
            let mut local_mut_unmod_bg = Vec::new();
            let mut local_pairs = Vec::new();

            for &i in &chunk {
                let m = &mod_records[i];
                let key = format!("{}|{}|{}", m[0], m[1], m[2]);
                
                // Calculate mod and unmod sample mutation rates
                // Use effective_depth if available (index 7), otherwise use depth (index 6)
                let (m_mut, m_stop) = parse_signal(m, 4, 5, 6, effective_depth_idx);
                
                // SNP filtering: filter positions where unmod sample mutation rate >= cutoff
                // This matches the strategy used by RNA Framework and ShapeMapper2
                // Unmod samples should have low mutation rates; high rates indicate SNPs or sequencing errors
                let is_filtered = if has_unmod {
                    if let Some(&u_idx) = unmod_key_map.get(&key) {
                        let u = &unmod_records[u_idx];
                        let (u_mut, _u_stop) = parse_signal(u, 4, 5, 6, effective_depth_idx);
                        u_mut >= snp_cutoff
                    } else {
                        false // If no matching unmod record, don't filter
                    }
                } else {
                    false // If no unmod sample provided, don't filter (mod-only mode)
                };
                
                if has_unmod {
                    if let Some(&u_idx) = unmod_key_map.get(&key) {
                        let u = &unmod_records[u_idx];
                        let (u_mut, u_stop) = parse_signal(u, 4, 5, 6, effective_depth_idx);
                        
                        // Only unfiltered positions are used for background region k-factor calculation
                        if !is_filtered {
                            let is_bg = if let Some(ref gene_id) = args.k_background_gene_id {
                                // Use specified gene_id for background regions
                                &m[0] == gene_id
                            } else {
                                // Default background region detection
                                m[0].contains("rRNA")
                                    || m[0].contains("16S")
                                    || m[0].contains("18S")
                                    || args.mod_csv.contains("16S")
                                    || args.mod_csv.contains("18S")
                            };
                            if is_bg {
                                local_stop_mod_bg.push(m_stop);
                                local_stop_unmod_bg.push(u_stop);
                                local_mut_mod_bg.push(m_mut);
                                local_mut_unmod_bg.push(u_mut);
                            }
                        }
                        
                        // Use f64::NAN to mark filtered positions, but keep position information
                        if is_filtered {
                            local_pairs.push((i, f64::NAN, f64::NAN, f64::NAN, f64::NAN));
                        } else {
                            local_pairs.push((i, m_mut, m_stop, u_mut, u_stop));
                        }
                    }
                } else {
                    // Mod-only mode: no unmod sample, so no filtering applied
                    // Only unfiltered positions are used for background region k-factor calculation
                    if !is_filtered {
                        let is_bg = if let Some(ref gene_id) = args.k_background_gene_id {
                            // Use specified gene_id for background regions
                            &m[0] == gene_id
                        } else {
                            // Default background region detection
                            m[0].contains("rRNA")
                                || m[0].contains("16S")
                                || m[0].contains("18S")
                                || args.mod_csv.contains("16S")
                                || args.mod_csv.contains("18S")
                        };
                        if is_bg {
                            local_stop_mod_bg.push(m_stop);
                            local_mut_mod_bg.push(m_mut);
                        }
                    }
                    
                    // Use f64::NAN to mark filtered positions, but keep position information
                    if is_filtered {
                        local_pairs.push((i, f64::NAN, f64::NAN, 0.0, 0.0));
                    } else {
                        local_pairs.push((i, m_mut, m_stop, 0.0, 0.0));
                    }
                }
            }

            (
                local_stop_mod_bg,
                local_stop_unmod_bg,
                local_mut_mod_bg,
                local_mut_unmod_bg,
                local_pairs,
            )
        })
        .collect();

    // Merge results
    for (local_stop_mod, local_stop_unmod, local_mut_mod, local_mut_unmod, local_pairs) in results {
        stop_mod_bg.extend(local_stop_mod);
        stop_unmod_bg.extend(local_stop_unmod);
        mut_mod_bg.extend(local_mut_mod);
        mut_unmod_bg.extend(local_mut_unmod);
        all_pairs.extend(local_pairs);
    }

    progress.finish()?;
    println!("\r[Data info]                      ");
    println!(
        "    mod records: {}, unmod records: {}",
        mod_records.len(),
        unmod_records.len()
    );
    println!(
        "    matched positions: {} (after SNP filtering, cutoff={:.2}), background regions: {}",
        all_pairs.len(),
        args.snp_cutoff,
        stop_mod_bg.len()
    );

    // Calculate k-factor
    print!("\r[Calculating] Computing k coefficients...");
    std::io::stdout().flush()?;
    let mut k_stop = 1.0;
    let mut k_mut = 1.0;

    let use_distribution_method = args.k_prediction_method.to_lowercase() == "distribution";
    let use_recursive_method = args.k_prediction_method.to_lowercase() == "recursive";

    if use_distribution_method && has_unmod {
        // Use distribution-based k prediction
        // Collect all mutation rates for distribution analysis
        let mut mod_mut_rates = Vec::new();
        let mut unmod_mut_rates = Vec::new();
        let mut mod_stop_rates = Vec::new();
        let mut unmod_stop_rates = Vec::new();

        for (_, m_mut, m_stop, u_mut, u_stop) in &all_pairs {
            if !m_mut.is_nan() && !u_mut.is_nan() {
                mod_mut_rates.push(*m_mut);
                unmod_mut_rates.push(*u_mut);
            }
            if !m_stop.is_nan() && !u_stop.is_nan() {
                mod_stop_rates.push(*m_stop);
                unmod_stop_rates.push(*u_stop);
            }
        }

        // Predict k_mut using distribution method
        if let Some(predicted_k_mut) = predict_k_by_distribution(&mod_mut_rates, &unmod_mut_rates, (0.0, 2.0), 0.01) {
            k_mut = predicted_k_mut;
            logger.log(&format!("k_mut predicted using distribution method: {:.4}", k_mut))?;
            } else {
                // Fallback to background method if distribution prediction fails
                if !mut_mod_bg.is_empty() && !mut_unmod_bg.is_empty() {
                    let mod_mean = mut_mod_bg.iter().sum::<f64>() / mut_mod_bg.len() as f64;
                    let unmod_mean = mut_unmod_bg.iter().sum::<f64>() / mut_unmod_bg.len() as f64;
                    if unmod_mean > 0.0 {
                        k_mut = mod_mean / unmod_mean;
                        logger.log(&format!("k_mut distribution prediction failed, using background method: {:.4}", k_mut))?;
                    } else {
                        logger.log("k_mut distribution prediction failed, using background method (k_mut=1.0 default)")?;
                    }
                } else {
                    logger.log("k_mut distribution prediction failed, using background method (k_mut=1.0 default)")?;
                }
            }

        // Predict k_stop using distribution method
        if let Some(predicted_k_stop) = predict_k_by_distribution(&mod_stop_rates, &unmod_stop_rates, (0.0, 2.0), 0.01) {
            k_stop = predicted_k_stop;
            logger.log(&format!("k_stop predicted using distribution method: {:.4}", k_stop))?;
            } else {
                // Fallback to background method if distribution prediction fails
                if !stop_mod_bg.is_empty() && !stop_unmod_bg.is_empty() {
                    let mod_mean = stop_mod_bg.iter().sum::<f64>() / stop_mod_bg.len() as f64;
                    let unmod_mean = stop_unmod_bg.iter().sum::<f64>() / stop_unmod_bg.len() as f64;
                    if unmod_mean > 0.0 {
                        k_stop = mod_mean / unmod_mean;
                        logger.log(&format!("k_stop distribution prediction failed, using background method: {:.4}", k_stop))?;
                    } else {
                        logger.log("k_stop distribution prediction failed, using background method (k_stop=1.0 default)")?;
                    }
                } else {
                    logger.log("k_stop distribution prediction failed, using background method (k_stop=1.0 default)")?;
                }
            }
    } else if use_recursive_method && has_unmod {
        // Use recursive k prediction with parallel two-phase search
        // Check if structure file is provided for structure-based recursive method
        if let Some(ref structure_file) = args.structure_file {
            // Load secondary structure
            let secondary_structure: SecondaryStructure = parse_secondary_structure(structure_file)?;
            logger.log(&format!("Loaded secondary structure from: {}", structure_file))?;
            logger.log(&format!("Structure length: {} positions", secondary_structure.positions.len()))?;

            // Collect rates and positions for structure-based recursive analysis
            let mut mod_mut_rates = Vec::new();
            let mut unmod_mut_rates = Vec::new();
            let mut mod_stop_rates = Vec::new();
            let mut unmod_stop_rates = Vec::new();
            let mut mut_positions = Vec::new();
            let mut stop_positions = Vec::new();
            let mut mut_is_single = Vec::new();
            let mut stop_is_single = Vec::new();

            // Extract gene ID from structure file name to filter data
            // Structure file name format: ref/{gene_id}.dp (e.g., ref/EC_16S.dp)
            // This is the reference gene used for k-value calculation (Human_18S for Human, EC_16S for EC)
            // The calculated k value will be applied to all genes in the reactivity file
            let structure_gene_id = std::path::Path::new(structure_file)
                .file_stem()
                .and_then(|s| s.to_str())
                .unwrap_or("");
            
            logger.log(&format!("Using reference gene '{}' for k-value calculation. The calculated k value will be applied to all genes.", structure_gene_id))?;
            
            // Count total genes in reactivity data
            let mut all_gene_ids = std::collections::HashSet::new();
            for (mod_idx, _, _, _, _) in &all_pairs {
                let m = &mod_records[*mod_idx];
                let chr_id = &m[0];
                all_gene_ids.insert(chr_id.clone());
            }
            let gene_list: Vec<String> = all_gene_ids.iter().map(|s| s.clone()).collect();
            logger.log(&format!("Reactivity data contains {} gene(s): {}", 
                all_gene_ids.len(), 
                gene_list.join(", ")))?;
            
            for (mod_idx, m_mut, m_stop, u_mut, u_stop) in &all_pairs {
                let m = &mod_records[*mod_idx];
                let chr_id = &m[0]; // Gene ID from reactivity data
                let position: u32 = m[2].parse().unwrap_or(0);
                
                // Only match positions from the reference gene (used for k-value calculation)
                // Other genes will use the same k value calculated from the reference gene
                if chr_id != structure_gene_id {
                    continue;
                }
                
                // Find position in structure
                if let Some(struct_idx) = secondary_structure.positions.iter().position(|&p| p == position) {
                    let is_single = secondary_structure.is_single_stranded[struct_idx];
                    
                    if !m_mut.is_nan() && !u_mut.is_nan() {
                        mod_mut_rates.push(*m_mut);
                        unmod_mut_rates.push(*u_mut);
                        mut_positions.push(position);
                        mut_is_single.push(is_single);
                    }
                    if !m_stop.is_nan() && !u_stop.is_nan() {
                        mod_stop_rates.push(*m_stop);
                        unmod_stop_rates.push(*u_stop);
                        stop_positions.push(position);
                        stop_is_single.push(is_single);
                    }
                }
            }

            logger.log(&format!("Matched {} mutation positions, {} stop positions from reference gene '{}' for k-value calculation", 
                mut_positions.len(), stop_positions.len(), structure_gene_id))?;

            // Predict k_mut using structure-based recursive method
            if let Some(predicted_k_mut) = predict_k_by_recursive_with_structure(
                &mod_mut_rates,
                &unmod_mut_rates,
                &mut_positions,
                &mut_is_single,
                (0.0, 2.0),
                0.1,  // coarse_step
                0.01, // fine_step
                0.1,  // fine_range
            ) {
                k_mut = predicted_k_mut;
                // Calculate and log metrics for the selected k value
                if let Some(metrics) = calculate_metrics_for_k(k_mut, &mod_mut_rates, &unmod_mut_rates, &mut_is_single) {
                    logger.log(&format!("k_mut predicted using structure-based recursive method: {:.4}", k_mut))?;
                    logger.log(&format!("  Metrics: AUC={:.4}, Sensitivity={:.4}, Specificity={:.4}, NegativeRatio={:.4}", 
                        metrics.auc, metrics.sensitivity, metrics.specificity, metrics.negative_ratio))?;
                    logger.log(&format!("  Zero values: Double-stranded={}/{} ({:.2}%), Single-stranded={}/{} ({:.2}%)", 
                        metrics.zero_count_double_stranded,
                        mod_mut_rates.len() - mut_is_single.iter().filter(|&&x| x).count(),
                        metrics.zero_ratio_double_stranded * 100.0,
                        metrics.zero_count_single_stranded,
                        mut_is_single.iter().filter(|&&x| x).count(),
                        metrics.zero_ratio_single_stranded * 100.0))?;
                } else {
                    logger.log(&format!("k_mut predicted using structure-based recursive method: {:.4}", k_mut))?;
                }
            } else {
                // Fallback to background method if recursive prediction fails
                if !mut_mod_bg.is_empty() && !mut_unmod_bg.is_empty() {
                    let mod_mean = mut_mod_bg.iter().sum::<f64>() / mut_mod_bg.len() as f64;
                    let unmod_mean = mut_unmod_bg.iter().sum::<f64>() / mut_unmod_bg.len() as f64;
                    if unmod_mean > 0.0 {
                        k_mut = mod_mean / unmod_mean;
                        logger.log(&format!("k_mut structure-based recursive prediction failed, using background method: {:.4}", k_mut))?;
                    } else {
                        logger.log("k_mut structure-based recursive prediction failed, using background method (k_mut=1.0 default)")?;
                    }
                } else {
                    logger.log("k_mut structure-based recursive prediction failed, using background method (k_mut=1.0 default)")?;
                }
            }

            // Predict k_stop using structure-based recursive method
            if let Some(predicted_k_stop) = predict_k_by_recursive_with_structure(
                &mod_stop_rates,
                &unmod_stop_rates,
                &stop_positions,
                &stop_is_single,
                (0.0, 2.0),
                0.1,  // coarse_step
                0.01, // fine_step
                0.1,  // fine_range
            ) {
                k_stop = predicted_k_stop;
                // Calculate and log metrics for the selected k value
                if let Some(metrics) = calculate_metrics_for_k(k_stop, &mod_stop_rates, &unmod_stop_rates, &stop_is_single) {
                    logger.log(&format!("k_stop predicted using structure-based recursive method: {:.4}", k_stop))?;
                    logger.log(&format!("  Metrics: AUC={:.4}, Sensitivity={:.4}, Specificity={:.4}, NegativeRatio={:.4}", 
                        metrics.auc, metrics.sensitivity, metrics.specificity, metrics.negative_ratio))?;
                    logger.log(&format!("  Zero values: Double-stranded={}/{} ({:.2}%), Single-stranded={}/{} ({:.2}%)", 
                        metrics.zero_count_double_stranded,
                        mod_stop_rates.len() - stop_is_single.iter().filter(|&&x| x).count(),
                        metrics.zero_ratio_double_stranded * 100.0,
                        metrics.zero_count_single_stranded,
                        stop_is_single.iter().filter(|&&x| x).count(),
                        metrics.zero_ratio_single_stranded * 100.0))?;
                } else {
                    logger.log(&format!("k_stop predicted using structure-based recursive method: {:.4}", k_stop))?;
                }
            } else {
                // Fallback to background method if recursive prediction fails
                if !stop_mod_bg.is_empty() && !stop_unmod_bg.is_empty() {
                    let mod_mean = stop_mod_bg.iter().sum::<f64>() / stop_mod_bg.len() as f64;
                    let unmod_mean = stop_unmod_bg.iter().sum::<f64>() / stop_unmod_bg.len() as f64;
                    if unmod_mean > 0.0 {
                        k_stop = mod_mean / unmod_mean;
                        logger.log(&format!("k_stop structure-based recursive prediction failed, using background method: {:.4}", k_stop))?;
                    } else {
                        logger.log("k_stop structure-based recursive prediction failed, using background method (k_stop=1.0 default)")?;
                    }
                } else {
                    logger.log("k_stop structure-based recursive prediction failed, using background method (k_stop=1.0 default)")?;
                }
            }
        } else {
            // Fallback to distribution-based recursive method if no structure file provided
            logger.log("No structure file provided, using distribution-based recursive method")?;
            
            // Collect all mutation rates for recursive analysis
            let mut mod_mut_rates = Vec::new();
            let mut unmod_mut_rates = Vec::new();
            let mut mod_stop_rates = Vec::new();
            let mut unmod_stop_rates = Vec::new();

            for (_, m_mut, m_stop, u_mut, u_stop) in &all_pairs {
                if !m_mut.is_nan() && !u_mut.is_nan() {
                    mod_mut_rates.push(*m_mut);
                    unmod_mut_rates.push(*u_mut);
                }
                if !m_stop.is_nan() && !u_stop.is_nan() {
                    mod_stop_rates.push(*m_stop);
                    unmod_stop_rates.push(*u_stop);
                }
            }

            // Predict k_mut using distribution-based recursive method
            if let Some(predicted_k_mut) = predict_k_by_recursive(
                &mod_mut_rates,
                &unmod_mut_rates,
                (0.0, 2.0),
                0.1,  // coarse_step
                0.01, // fine_step
                0.1,  // fine_range
                thread_count,
            ) {
                k_mut = predicted_k_mut;
                logger.log(&format!("k_mut predicted using distribution-based recursive method: {:.4}", k_mut))?;
            } else {
                // Fallback to background method if recursive prediction fails
                if !mut_mod_bg.is_empty() && !mut_unmod_bg.is_empty() {
                    let mod_mean = mut_mod_bg.iter().sum::<f64>() / mut_mod_bg.len() as f64;
                    let unmod_mean = mut_unmod_bg.iter().sum::<f64>() / mut_unmod_bg.len() as f64;
                    if unmod_mean > 0.0 {
                        k_mut = mod_mean / unmod_mean;
                        logger.log(&format!("k_mut distribution-based recursive prediction failed, using background method: {:.4}", k_mut))?;
                    } else {
                        logger.log("k_mut distribution-based recursive prediction failed, using background method (k_mut=1.0 default)")?;
                    }
                } else {
                    logger.log("k_mut distribution-based recursive prediction failed, using background method (k_mut=1.0 default)")?;
                }
            }

            // Predict k_stop using distribution-based recursive method
            if let Some(predicted_k_stop) = predict_k_by_recursive(
                &mod_stop_rates,
                &unmod_stop_rates,
                (0.0, 2.0),
                0.1,  // coarse_step
                0.01, // fine_step
                0.1,  // fine_range
                thread_count,
            ) {
                k_stop = predicted_k_stop;
                logger.log(&format!("k_stop predicted using distribution-based recursive method: {:.4}", k_stop))?;
            } else {
                // Fallback to background method if recursive prediction fails
                if !stop_mod_bg.is_empty() && !stop_unmod_bg.is_empty() {
                    let mod_mean = stop_mod_bg.iter().sum::<f64>() / stop_mod_bg.len() as f64;
                    let unmod_mean = stop_unmod_bg.iter().sum::<f64>() / stop_unmod_bg.len() as f64;
                    if unmod_mean > 0.0 {
                        k_stop = mod_mean / unmod_mean;
                        logger.log(&format!("k_stop distribution-based recursive prediction failed, using background method: {:.4}", k_stop))?;
                    } else {
                        logger.log("k_stop distribution-based recursive prediction failed, using background method (k_stop=1.0 default)")?;
                    }
                } else {
                    logger.log("k_stop distribution-based recursive prediction failed, using background method (k_stop=1.0 default)")?;
                }
            }
        }
    } else {
        // Use traditional background region method (default)
        if !stop_mod_bg.is_empty() && !stop_unmod_bg.is_empty() {
            let mod_mean = stop_mod_bg.iter().sum::<f64>() / stop_mod_bg.len() as f64;
            let unmod_mean = stop_unmod_bg.iter().sum::<f64>() / stop_unmod_bg.len() as f64;

            if unmod_mean > 0.0 {
                k_stop = mod_mean / unmod_mean;
                logger.log(&format!("k_stop calculated using background region method: {:.4}", k_stop))?;
            }
        }

        if !mut_mod_bg.is_empty() && !mut_unmod_bg.is_empty() {
            let mod_mean = mut_mod_bg.iter().sum::<f64>() / mut_mod_bg.len() as f64;
            let unmod_mean = mut_unmod_bg.iter().sum::<f64>() / mut_unmod_bg.len() as f64;

            if unmod_mean > 0.0 {
                k_mut = mod_mean / unmod_mean;
                logger.log(&format!("k_mut calculated using background region method: {:.4}", k_mut))?;
            }
        }
    }

    println!("\r[k value]                                ");
    if use_distribution_method {
        println!("    Method: distribution-based prediction");
    } else if use_recursive_method {
        if args.structure_file.is_some() {
            println!("    Method: recursive (structure-based, parallel two-phase search)");
        } else {
            println!("    Method: recursive (distribution-based, parallel two-phase search)");
        }
    } else {
        if let Some(ref gene_id) = args.k_background_gene_id {
            println!("    Method: background region ratio (gene_id: {})", gene_id);
        } else {
            println!("    Method: background region ratio (default regions)");
        }
    }
    println!(
        "    stop signal mean: mod={:.4}, unmod={:.4}. k_stop = {:.3}",
        stop_mod_bg.iter().sum::<f64>() / stop_mod_bg.len().max(1) as f64,
        stop_unmod_bg.iter().sum::<f64>() / stop_unmod_bg.len().max(1) as f64,
        k_stop
    );
    println!(
        "    mutation signal mean: mod={:.4}, unmod={:.4}. k_mut = {:.3}",
        mut_mod_bg.iter().sum::<f64>() / mut_mod_bg.len().max(1) as f64,
        mut_unmod_bg.iter().sum::<f64>() / mut_unmod_bg.len().max(1) as f64,
        k_mut
    );
    
    // Log that k values will be applied to all genes
    if use_recursive_method && args.structure_file.is_some() {
        if let Some(ref structure_file) = args.structure_file {
            let structure_gene_id = std::path::Path::new(structure_file)
                .file_stem()
                .and_then(|s| s.to_str())
                .unwrap_or("");
            logger.log(&format!("k values (k_mut={:.4}, k_stop={:.4}) calculated from reference gene '{}' will be applied to all genes in the reactivity calculation", 
                k_mut, k_stop, structure_gene_id))?;
        }
    }

    // 2. Calculate reactivity in parallel
    let total_pairs = all_pairs.len();
    let mut calc_progress = crate::progress::SimpleProgress::new(total_pairs);

    print!("\r[Reactivity calculating] ");
    std::io::stdout().flush()?;

    // Calculate reactivity in parallel
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
                // If input value is NaN (SNP filtered position), directly output NaN
                if m_mut.is_nan() || m_stop.is_nan() {
                    local_results.push((f64::NAN, f64::NAN));
                    continue;
                }
                
                // Stop signal method selection
                let stop_reactivity = match stop_method {
                    StopReactivityMethod::Kfactor => {
                        calculate_k_corrected_reactivity(m_stop, u_stop, k_stop)
                    }
                    StopReactivityMethod::Ding => calculate_reactivity_ding_style(
                        m_stop,
                        u_stop,
                        k_stop,
                        args.pseudocount,
                        args.maxscore,
                    ),
                    StopReactivityMethod::Rouskin => m_stop, // Use only mod sample stop signal
                };

                // Mutation signal method selection
                let mut_reactivity = match mutation_method {
                    MutationReactivityMethod::Kfactor => {
                        calculate_k_corrected_reactivity(m_mut, u_mut, k_mut)
                    }
                    MutationReactivityMethod::Siegfried => {
                        // Siegfried method: S_i = (n_Ti / c_Ti) - (n_Ui / c_Ui)
                        // Directly calculate mutation rate difference, not using k-factor correction
                        // Allow negative values (as in Shapemapper) - negative values indicate
                        // background mutation rate is higher than modified sample
                        m_mut - u_mut
                    }
                    MutationReactivityMethod::Zubradt => m_mut, // Use only mod sample mutation signal
                };

                local_results.push((stop_reactivity, mut_reactivity));
            }

            local_results
        })
        .collect();

    calc_progress.finish()?;

    // 3. Generate output file in parallel
    // Ensure output directory exists
    if let Some(parent) = Path::new(&args.output).parent() {
        std::fs::create_dir_all(parent)?;
    }
    let mut output_file = File::create(&args.output)?;

    // Write header
    let header_line = "ChrID,Strand,Position,Base,reactivity_stop,reactivity_mutation".to_string();
    writeln!(output_file, "{}", header_line)?;

    // Prepare output data in parallel
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

                // Get parallel calculated reactivity results
                let chunk_index = i / reactivity_chunk_size;
                let local_index = i % reactivity_chunk_size;
                let (stop_reactivity, mut_reactivity) =
                    reactivity_results[chunk_index][local_index];

                // Simplified output format: only keep core columns
                let chr_id = &rec[0];
                let strand = &rec[1];
                let position = &rec[2];
                let base = &rec[3];

                // Format output, NaN values directly output "NaN"
                let stop_str = if stop_reactivity.is_nan() {
                    "NaN".to_string()
                } else {
                    format!("{:.6}", stop_reactivity)
                };
                let mut_str = if mut_reactivity.is_nan() {
                    "NaN".to_string()
                } else {
                    format!("{:.6}", mut_reactivity)
                };
                
                let line = format!(
                    "{},{},{},{},{},{}",
                    chr_id, strand, position, base, stop_str, mut_str
                );
                lines.push(line);
            }

            lines
        })
        .collect();

    // Write files in parallel
    let mut write_progress = crate::progress::SimpleProgress::new(total_pairs);
    print!("\r[Writing files] ");
    std::io::stdout().flush()?;

    // Parallel write strategy: split output into multiple temporary files, then merge
    // Use unique temp directory per output file to avoid conflicts in parallel execution
    let output_filename = Path::new(&args.output)
        .file_name()
        .and_then(|n| n.to_str())
        .unwrap_or("reactivity");
    let temp_dir = std::env::temp_dir()
        .join("moddetector_reactivity")
        .join(format!("{}_{}", output_filename, std::process::id()));
    std::fs::create_dir_all(&temp_dir)?;

    let write_results: Vec<Vec<String>> = prepared_outputs
        .into_par_iter()
        .enumerate()
        .map(|(chunk_id, lines)| {
            let mut chunk_output = Vec::new();

            // Create temporary file
            let temp_file = temp_dir.join(format!("reactivity_chunk_{}.tmp", chunk_id));

            // Write to temporary file
            let mut temp_out = File::create(&temp_file).unwrap();

            for line in lines.into_iter() {
                writeln!(temp_out, "{}", line).unwrap();
                chunk_output.push(line);
            }

            // Close temporary file
            drop(temp_out);

            chunk_output
        })
        .collect();

    // Merge all temporary files to final output
    let mut processed_count = 0;
    for (chunk_id, chunk_output) in write_results.into_iter().enumerate() {
        // Read temporary file and write to final output
        let temp_file = temp_dir.join(format!("reactivity_chunk_{}.tmp", chunk_id));

        if temp_file.exists() {
            let content = std::fs::read_to_string(&temp_file)?;

            output_file.write_all(content.as_bytes())?;

            // Delete temporary file
            std::fs::remove_file(temp_file).ok();

            processed_count += chunk_output.len();
            if processed_count % 10000 == 0 {
                write_progress.update(processed_count)?;
            }
        }
    }

    // Clean up temporary directory
    std::fs::remove_dir_all(temp_dir).ok();

    write_progress.finish()?;

    let elapsed = start_time.elapsed();

    // Display output information based on actual method used
    let stop_method_name = match stop_method {
        StopReactivityMethod::Kfactor => "Kfactor",
        StopReactivityMethod::Ding => "Ding",
        StopReactivityMethod::Rouskin => "Rouskin",
    };

    let mutation_method_name = match mutation_method {
        MutationReactivityMethod::Kfactor => "Kfactor",
        MutationReactivityMethod::Siegfried => "Siegfried",
        MutationReactivityMethod::Zubradt => "Zubradt",
    };

    println!(
        "\r[Reactivity outputs] stop={}, mutation={}",
        stop_method_name, mutation_method_name
    );
    println!("    output={}", args.output);
    println!("{}", crate::progress::format_time_used(elapsed));

    // Log completion
    logger.log(&format!("Reactivity calculation completed"))?;
    logger.log(&format!("Reactivity Output: {}", args.output))?;
    logger.log(&format!("Stop Method: {}", stop_method_name))?;
    logger.log(&format!("Mutation Method: {}", mutation_method_name))?;
    if !has_unmod {
        logger.log("Mode: mod-only (no unmodified sample provided)")?;
    }
    logger.log(&format!("Total time: {:.2}s", elapsed.as_secs_f64()))?;

    Ok(())
}

/// Load CSV data in chunks in parallel (original implementation, kept as backup)
fn load_data_in_chunks(
    mod_csv: &str,
    unmod_csv: &str,
    thread_count: usize,
) -> Result<
    (
        Vec<Vec<String>>,
        Vec<Vec<String>>,
        HashMap<String, usize>,
        HashMap<String, usize>,
        String,
    ),
    Box<dyn Error>,
> {
    use std::sync::atomic::{AtomicUsize, Ordering};
    use std::sync::Arc;

    println!("[Loading] Starting parallel data loading...");

    // Read mod file
    let mod_file = File::open(mod_csv)?;
    let mut mod_lines: Vec<String> = BufReader::new(mod_file)
        .lines()
        .collect::<Result<Vec<String>, _>>()?;
    let header = mod_lines.remove(0);

    // Read unmod file
    let unmod_file = File::open(unmod_csv)?;
    let mut unmod_lines: Vec<String> = BufReader::new(unmod_file)
        .lines()
        .collect::<Result<Vec<String>, _>>()?;
    let _ = unmod_lines.remove(0); // Skip unmod header

    println!(
        "[Loading] Mod records: {}, Unmod records: {}",
        mod_lines.len(),
        unmod_lines.len()
    );

    // Process mod data in chunks
    let mod_chunk_size = (mod_lines.len() + thread_count - 1) / thread_count;
    let mod_chunks: Vec<Vec<String>> = mod_lines
        .chunks(mod_chunk_size)
        .map(|chunk| chunk.to_vec())
        .collect();

    // Process unmod data in chunks
    let unmod_chunk_size = (unmod_lines.len() + thread_count - 1) / thread_count;
    let unmod_chunks: Vec<Vec<String>> = unmod_lines
        .chunks(unmod_chunk_size)
        .map(|chunk| chunk.to_vec())
        .collect();

    let num_chunks = mod_chunks.len() + unmod_chunks.len();
    let completed = Arc::new(AtomicUsize::new(0));

    print!(
        "\r[Loading] Processing {} chunks with {} threads",
        num_chunks, thread_count
    );
    std::io::stdout().flush().ok();

    // Process mod data in parallel
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

            // Update progress
            let completed_count = completed.fetch_add(1, Ordering::Relaxed) + 1;
            let percentage = (completed_count as f64 * 100.0) / num_chunks as f64;
            print!(
                "\r[Loading] Processing {}/{} chunks ({:.1}%)",
                completed_count, num_chunks, percentage
            );
            std::io::stdout().flush().ok();

            (records, key_map)
        })
        .collect();

    // Process unmod data in parallel
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

            // Update progress
            let completed_count = completed.fetch_add(1, Ordering::Relaxed) + 1;
            let percentage = (completed_count as f64 * 100.0) / num_chunks as f64;
            print!(
                "\r[Loading] Processing {}/{} chunks ({:.1}%)",
                completed_count, num_chunks, percentage
            );
            std::io::stdout().flush().ok();

            (records, key_map)
        })
        .collect();

    println!("\r[Loading] Data loading completed");

    // Merge results
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

    println!(
        "[Loading] Loaded {} mod records, {} unmod records",
        mod_records.len(),
        unmod_records.len()
    );

    Ok((
        mod_records,
        unmod_records,
        mod_key_map,
        unmod_key_map,
        header,
    ))
}

fn load_mod_data_only(
    mod_csv: &str,
    thread_count: usize,
) -> Result<Vec<Vec<String>>, Box<dyn Error>> {
    use std::sync::atomic::{AtomicUsize, Ordering};
    use std::sync::Arc;

    println!("[Loading] Starting parallel data loading (mod-only)...");

    let mod_file = File::open(mod_csv)?;
    let mut mod_lines: Vec<String> = BufReader::new(mod_file)
        .lines()
        .collect::<Result<Vec<String>, _>>()?;
    if mod_lines.is_empty() {
        println!("[Loading] No data rows found in mod file");
        return Ok(Vec::new());
    }
    let _header = mod_lines.remove(0);

    println!("[Loading] Mod records: {}", mod_lines.len());

    if mod_lines.is_empty() {
        println!("[Loading] No mod data rows after header");
        return Ok(Vec::new());
    }

    let mod_chunk_size = (mod_lines.len() + thread_count - 1) / thread_count;
    let mod_chunks: Vec<Vec<String>> = mod_lines
        .chunks(mod_chunk_size.max(1))
        .map(|chunk| chunk.to_vec())
        .collect();

    let num_chunks = mod_chunks.len();
    let completed = Arc::new(AtomicUsize::new(0));

    print!(
        "\r[Loading] Processing {} mod chunks with {} threads",
        num_chunks, thread_count
    );
    std::io::stdout().flush().ok();

    let mod_results: Vec<Vec<Vec<String>>> = mod_chunks
        .into_par_iter()
        .map(|chunk| {
            let mut records = Vec::with_capacity(chunk.len());

            for line in &chunk {
                let fields: Vec<String> = line.split(',').map(|s| s.to_string()).collect();
                if fields.len() >= 3 {
                    records.push(fields);
                }
            }

            let completed_count = completed.fetch_add(1, Ordering::Relaxed) + 1;
            let percentage = if num_chunks == 0 {
                100.0
            } else {
                (completed_count as f64 * 100.0) / num_chunks as f64
            };
            print!(
                "\r[Loading] Processing {}/{} mod chunks ({:.1}%)",
                completed_count, num_chunks, percentage
            );
            std::io::stdout().flush().ok();

            records
        })
        .collect();

    println!("\r[Loading] Mod-only data loading completed");

    let mut mod_records = Vec::with_capacity(mod_lines.len());
    for records in mod_results {
        mod_records.extend(records);
    }

    println!(
        "[Loading] Loaded {} mod records (mod-only)",
        mod_records.len()
    );

    Ok(mod_records)
}
