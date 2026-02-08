use rand::seq::SliceRandom;
use rand::thread_rng;
use std::collections::HashMap;
use std::error::Error;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::time::Instant;
// use chrono;

// New: Comparison mode enumeration
#[derive(Debug, Clone)]
pub enum CompareMode {
    /// Original mod and unmod comparison
    ModVsUnmod,
    /// Reactivity comparison between different groups
    ReactivityGroups,
    /// Multiple biological replicates statistical testing
    BiologicalReplicates,
}

// New: Statistical test type enumeration
#[derive(Debug, Clone)]
pub enum StatisticalTest {
    /// Chi-square test
    ChiSquare,
    /// Continuity correction test
    ContinuityCorrection,
    /// t-test
    TTest,
    /// Mann-Whitney U test
    MannWhitneyU,
    /// Wilcoxon signed-rank test
    WilcoxonSignedRank,
    /// DiffScan method (normalization + Wilcoxon + scanning statistics)
    DiffScan,
    /// deltaSHAPE method (Z-factor + Z-score + sliding window)
    DeltaSHAPE,
}

// New: Statistical result structure
#[derive(Debug, Clone)]
pub struct StatisticalResult {
    pub test_type: StatisticalTest,
    pub statistic: f64,
    pub p_value: f64,
    pub significant: bool,
    pub effect_size: Option<f64>,
}

// New: Replicate data structure
#[derive(Debug, Clone)]
pub struct ReplicateData {
    pub group_name: String,
    pub replicates: Vec<PositionData>,
}

// New: Comparison result structure
#[derive(Debug, Clone)]
pub struct ComparisonResult {
    pub chr_id: String,
    pub strand: char,
    pub position: u32,
    pub ref_base: char,
    pub group1_values: Vec<f64>,
    pub group2_values: Vec<f64>,
    pub fold_change: f64,
    pub statistical_result: Option<StatisticalResult>,
    pub modification_score: f64,
}

// New: DiffScan related structure
#[derive(Debug, Clone)]
pub struct DiffScanResult {
    pub chr_id: String,
    pub strand: char,
    pub start_pos: u32,
    pub end_pos: u32,
    pub normalized_reactivity1: f64,
    pub normalized_reactivity2: f64,
    pub wilcoxon_p_value: f64,
    pub scan_statistic: f64,
    pub monte_carlo_p_value: f64,
    pub significant: bool,
    pub region_size: usize,
}

// deltaSHAPE result structure
#[derive(Debug, Clone)]
pub struct DeltaSHAPEResult {
    pub chr_id: String,
    pub strand: char,
    pub position: u32,
    pub base: char,
    pub smoothed_reactivity1: f64,
    pub smoothed_reactivity2: f64,
    pub delta_reactivity: f64,
    pub z_factor: f64,
    pub z_score: f64,
    pub significant: bool,
}

// New: Normalization parameters
#[derive(Debug, Clone)]
pub struct NormalizationParams {
    pub method: String,              // "quantile", "robust", "median"
    pub reference_condition: String, // Reference condition
}

#[derive(Debug, Clone)]
struct PositionData {
    chr_id: String,
    strand: char,
    position: u32,
    ref_base: char,
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

// New: Load reactivity data
fn load_reactivity_data(file_path: &str) -> Result<HashMap<String, Vec<f64>>, Box<dyn Error>> {
    let mut data = HashMap::new();
    let file = File::open(file_path)?;
    let reader = BufReader::new(file);

    for (line_num, line) in reader.lines().enumerate() {
        if line_num == 0 {
            continue; // Skip header row
        }

        let line = line?;
        let fields: Vec<&str> = line.split(',').collect();

        if fields.len() >= 4 {
            let chr_id = fields[0].trim_matches('"').to_string();
            let strand = fields[1].chars().next().unwrap_or('+');
            let position = fields[2].parse::<u32>()?;

            // Support two formats:
            // 1. ChrID,Strand,Position,Reactivity (4 columns)
            // 2. ChrID,Strand,Position,Base,reactivity_stop,reactivity_mutation (6 columns)
            let reactivity = if fields.len() >= 6 {
                // 6-column format: use reactivity_stop column (column 5, index 4)
                fields[4].parse::<f64>()?
            } else {
                // 4-column format: use column 4 (index 3)
                fields[3].parse::<f64>()?
            };

            let key = format!("{}:{}{}", chr_id, strand, position);
            data.entry(key).or_insert_with(Vec::new).push(reactivity);
        }
    }

    Ok(data)
}

// New: Perform statistical test
fn perform_statistical_test(
    group1: &[f64],
    group2: &[f64],
    test_type: &StatisticalTest,
    alpha: f64,
) -> StatisticalResult {
    match test_type {
        StatisticalTest::TTest => perform_t_test(group1, group2, alpha),
        StatisticalTest::MannWhitneyU => perform_mann_whitney_u_test(group1, group2, alpha),
        StatisticalTest::WilcoxonSignedRank => perform_wilcoxon_test(group1, group2, alpha),
        StatisticalTest::ChiSquare => perform_chi_square_test(group1, group2, alpha),
        StatisticalTest::ContinuityCorrection => {
            perform_continuity_correction_test(group1, group2, alpha)
        }
        StatisticalTest::DiffScan => perform_diffscan_test(group1, group2, alpha),
        StatisticalTest::DeltaSHAPE => perform_deltashape_test(group1, group2, alpha),
    }
}

// New: DiffScan statistical test
fn perform_diffscan_test(group1: &[f64], group2: &[f64], alpha: f64) -> StatisticalResult {
    // For single value comparison, use simplified DiffScan method (no internal normalization/clipping, fully trust upstream normalized input)
    if group1.len() == 1 && group2.len() == 1 {
        let reactivity1 = group1[0];
        let reactivity2 = group2[0];

        // Directly use input values to calculate difference
        let normalized_diff = (reactivity1 - reactivity2).abs();

        // Stricter significance criteria to reduce false positives (threshold can be parameterized later)
        let is_significant = normalized_diff > 0.3
            && (reactivity1 > 0.2 || reactivity2 > 0.2)
            && (reactivity1 > 0.05 && reactivity2 > 0.05);

        let wilcoxon_p = if is_significant {
            if normalized_diff > 0.6 {
                0.001
            } else if normalized_diff > 0.4 {
                0.01
            } else {
                0.05
            }
        } else {
            0.5
        };

        StatisticalResult {
            test_type: StatisticalTest::DiffScan,
            statistic: normalized_diff,
            p_value: wilcoxon_p,
            significant: is_significant,
            effect_size: Some(normalized_diff),
        }
    } else {
        // Multi-value case uses standard Wilcoxon test
        perform_wilcoxon_test(group1, group2, alpha)
    }
}

// New: deltaSHAPE statistical test
// Note: deltaSHAPE uses DiffScan logic at single position level, only uses deltaSHAPE-specific region detection at region level
fn perform_deltashape_test(group1: &[f64], group2: &[f64], alpha: f64) -> StatisticalResult {
    // deltaSHAPE uses DiffScan logic at single position level (already reasonable)
    // Region-level deltaSHAPE detection is implemented separately in compare_reactivity_results function
    perform_diffscan_test(group1, group2, alpha)
}

// New: Calculate Wilcoxon position p-values based on sliding window
fn compute_window_pvalues(
    series1: &[(u32, f64)],
    series2: &[(u32, f64)],
    window_radius: usize,
) -> Vec<(u32, f64)> {
    // Assume both conditions share the same position set, input is sorted by position in ascending order
    let mut pos_to_idx1 = std::collections::HashMap::new();
    for (i, (pos, _)) in series1.iter().enumerate() {
        pos_to_idx1.insert(*pos, i);
    }
    let mut pos_to_idx2 = std::collections::HashMap::new();
    for (i, (pos, _)) in series2.iter().enumerate() {
        pos_to_idx2.insert(*pos, i);
    }

    let len = series1.len().min(series2.len());
    let mut results = Vec::with_capacity(len);
    for i in 0..len {
        let center_pos = series1[i].0;
        if !pos_to_idx2.contains_key(&center_pos) {
            continue;
        }

        // Take window range [i - r, i + r]
        let start1 = i.saturating_sub(window_radius);
        let end1 = (i + window_radius).min(series1.len() - 1);

        let j = *pos_to_idx2.get(&center_pos).unwrap();
        let start2 = j.saturating_sub(window_radius);
        let end2 = (j + window_radius).min(series2.len() - 1);

        let win1: Vec<f64> = series1[start1..=end1].iter().map(|(_, v)| *v).collect();
        let win2: Vec<f64> = series2[start2..=end2].iter().map(|(_, v)| *v).collect();

        let stat = perform_wilcoxon_test(&win1, &win2, 0.05);
        results.push((center_pos, stat.p_value));
    }
    results
}

// New: Calculate region scan statistic Q(R) = -sum(log(pj)) / sqrt(|R|)
fn compute_scan_Q(
    pos_pvals: &[(u32, f64)],
    min_len: usize,
    max_len: usize,
) -> Vec<(u32, u32, f64)> {
    let n = pos_pvals.len();
    if n == 0 {
        return vec![];
    }
    // Prefix sum: -log(p)
    let mut prefix: Vec<f64> = Vec::with_capacity(n + 1);
    prefix.push(0.0);
    for &(_, p) in pos_pvals {
        let v = if p > 0.0 { -(p.ln()) } else { 0.0 };
        prefix.push(prefix.last().unwrap() + v);
    }

    let mut regions = Vec::new();
    for len in min_len..=max_len {
        if len == 0 {
            continue;
        }
        if len > n {
            continue;
        }
        for i in 0..=n.saturating_sub(len) {
            let j = i + len;
            let s = prefix[j] - prefix[i];
            let q = s / (len as f64).sqrt();
            let start_pos = pos_pvals[i].0;
            let end_pos = pos_pvals[j - 1].0;
            regions.push((start_pos, end_pos, q));
        }
    }
    regions
}

// New: Monte Carlo evaluation threshold (FWER control). Returns Q critical value at given alpha
fn monte_carlo_q_threshold(
    pos_pvals: &[(u32, f64)],
    min_len: usize,
    max_len: usize,
    alpha: f64,
    iterations: usize,
) -> f64 {
    if pos_pvals.is_empty() {
        return f64::INFINITY;
    }
    // Extract p-values and permute to simulate null hypothesis
    let mut rng = thread_rng();
    let mut base: Vec<f64> = pos_pvals.iter().map(|(_, p)| *p).collect();
    let n = base.len();
    let mut maxima: Vec<f64> = Vec::with_capacity(iterations);
    for _ in 0..iterations {
        base.shuffle(&mut rng);
        let shuffled: Vec<(u32, f64)> = pos_pvals
            .iter()
            .enumerate()
            .map(|(i, (pos, _))| (*pos, base[i]))
            .collect();
        let regions = compute_scan_Q(&shuffled, min_len, max_len);
        // Take maximum Q value of regions in this permutation
        let m = regions
            .iter()
            .map(|r| r.2)
            .fold(f64::NEG_INFINITY, f64::max);
        if m.is_finite() {
            maxima.push(m);
        }
    }
    // Take (1-alpha) quantile as threshold
    maxima.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let idx = ((1.0 - alpha) * (maxima.len() as f64 - 1.0)).round() as usize;
    maxima.get(idx).cloned().unwrap_or(f64::INFINITY)
}

// New: deltaSHAPE region detection function
// Returns (start_pos, end_pos, z_factor, z_score, significant)
fn compute_deltashape_regions(
    series1: &[(u32, f64)],
    series2: &[(u32, f64)],
    window_radius: usize,
    min_len: usize,
    max_len: usize,
) -> Vec<(u32, u32, f64, f64, bool)> {
    let mut regions = Vec::new();

    // Calculate Z-factor and Z-score for each position
    let mut pos_stats = Vec::new();
    for i in 0..series1.len() {
        let (pos, v1, v2) = (series1[i].0, series1[i].1, series2[i].1);

        // Calculate difference
        let delta_reactivity = v1 - v2;

        // Simplified Z-factor calculation (assume error is 10% of reactivity)
        let error1 = v1 * 0.1;
        let error2 = v2 * 0.1;
        let z_factor = if delta_reactivity.abs() == 0.0 {
            0.0
        } else {
            1.0 - (1.96 * (error1 + error2) / delta_reactivity.abs())
        };

        // Simplified Z-score calculation
        let z_score = if delta_reactivity.abs() > 0.0 {
            delta_reactivity.abs() / 0.1 // Assume standard deviation is 0.1
        } else {
            0.0
        };

        pos_stats.push((pos, z_factor, z_score));
    }

    // Find continuous regions using sliding window method
    let mut i = 0;
    while i < pos_stats.len() {
        let mut region_start = i;
        let mut region_end = i;
        let mut max_z_factor = pos_stats[i].1;
        let mut max_z_score = pos_stats[i].2;
        let mut significant_count = 0;

        // Expand region, find continuous significant sites
        let mut j = i;
        while j < pos_stats.len() && (j - i) < max_len {
            let (_, z_factor, z_score) = pos_stats[j];

            // deltaSHAPE significance judgment: Z-factor > 0 and |Z-score| >= 1
            let is_significant = z_factor > 0.0 && z_score.abs() >= 1.0;

            if is_significant {
                significant_count += 1;
                region_end = j;
                max_z_factor = max_z_factor.max(z_factor);
                max_z_score = max_z_score.max(z_score.abs());
            } else if significant_count > 0 {
                // If there are already significant sites but current site is not significant, check if region should end
                break;
            }

            j += 1;
        }

        // If region length meets requirements and has enough significant sites, add region
        let region_len = region_end - region_start + 1;
        if region_len >= min_len && significant_count >= (region_len / 2) {
            regions.push((
                pos_stats[region_start].0,
                pos_stats[region_end].0,
                max_z_factor,
                max_z_score,
                true,
            ));
        }

        i = region_end + 1;
    }

    regions
}

// New: Normalization function
fn normalize_reactivity_data(
    data1: &HashMap<String, Vec<f64>>,
    data2: &HashMap<String, Vec<f64>>,
    method: &str,
) -> (HashMap<String, Vec<f64>>, HashMap<String, Vec<f64>>) {
    match method {
        "quantile" => quantile_normalization(data1, data2),
        "robust" => robust_normalization(data1, data2),
        "median" => median_normalization(data1, data2),
        _ => (data1.clone(), data2.clone()),
    }
}

// New: Quantile normalization
fn quantile_normalization(
    data1: &HashMap<String, Vec<f64>>,
    data2: &HashMap<String, Vec<f64>>,
) -> (HashMap<String, Vec<f64>>, HashMap<String, Vec<f64>>) {
    // Collect all values
    let mut all_values1: Vec<f64> = data1.values().flatten().cloned().collect();
    let mut all_values2: Vec<f64> = data2.values().flatten().cloned().collect();

    all_values1.sort_by(|a, b| a.partial_cmp(b).unwrap());
    all_values2.sort_by(|a, b| a.partial_cmp(b).unwrap());

    // Calculate quantiles
    let mut normalized_data1 = HashMap::new();
    let mut normalized_data2 = HashMap::new();

    for (key, values1) in data1 {
        if let Some(values2) = data2.get(key) {
            let mut norm_values1 = Vec::new();
            let mut norm_values2 = Vec::new();

            for &val1 in values1 {
                let rank1 = all_values1
                    .iter()
                    .position(|&x| x >= val1)
                    .unwrap_or(all_values1.len());
                let percentile1 = rank1 as f64 / all_values1.len() as f64;
                let norm_val1 = all_values2[(percentile1 * all_values2.len() as f64) as usize];
                norm_values1.push(norm_val1);
            }

            for &val2 in values2 {
                let rank2 = all_values2
                    .iter()
                    .position(|&x| x >= val2)
                    .unwrap_or(all_values2.len());
                let percentile2 = rank2 as f64 / all_values2.len() as f64;
                let norm_val2 = all_values1[(percentile2 * all_values1.len() as f64) as usize];
                norm_values2.push(norm_val2);
            }

            normalized_data1.insert(key.clone(), norm_values1);
            normalized_data2.insert(key.clone(), norm_values2);
        }
    }

    (normalized_data1, normalized_data2)
}

// New: Robust normalization
fn robust_normalization(
    data1: &HashMap<String, Vec<f64>>,
    data2: &HashMap<String, Vec<f64>>,
) -> (HashMap<String, Vec<f64>>, HashMap<String, Vec<f64>>) {
    // Calculate median and MAD (median absolute deviation)
    let mut all_values1: Vec<f64> = data1.values().flatten().cloned().collect();
    let mut all_values2: Vec<f64> = data2.values().flatten().cloned().collect();

    all_values1.sort_by(|a, b| a.partial_cmp(b).unwrap());
    all_values2.sort_by(|a, b| a.partial_cmp(b).unwrap());

    let median1 = all_values1[all_values1.len() / 2];
    let median2 = all_values2[all_values2.len() / 2];

    let mad1 = calculate_mad(&all_values1, median1);
    let mad2 = calculate_mad(&all_values2, median2);

    // Normalize
    let mut normalized_data1 = HashMap::new();
    let mut normalized_data2 = HashMap::new();

    for (key, values1) in data1 {
        if let Some(values2) = data2.get(key) {
            let norm_values1: Vec<f64> = values1.iter().map(|&x| (x - median1) / mad1).collect();
            let norm_values2: Vec<f64> = values2.iter().map(|&x| (x - median2) / mad2).collect();

            normalized_data1.insert(key.clone(), norm_values1);
            normalized_data2.insert(key.clone(), norm_values2);
        }
    }

    (normalized_data1, normalized_data2)
}

// New: Median normalization
fn median_normalization(
    data1: &HashMap<String, Vec<f64>>,
    data2: &HashMap<String, Vec<f64>>,
) -> (HashMap<String, Vec<f64>>, HashMap<String, Vec<f64>>) {
    // Calculate median
    let mut all_values1: Vec<f64> = data1.values().flatten().cloned().collect();
    let mut all_values2: Vec<f64> = data2.values().flatten().cloned().collect();

    all_values1.sort_by(|a, b| a.partial_cmp(b).unwrap());
    all_values2.sort_by(|a, b| a.partial_cmp(b).unwrap());

    let median1 = all_values1[all_values1.len() / 2];
    let median2 = all_values2[all_values2.len() / 2];

    // Normalize to same median
    let target_median = (median1 + median2) / 2.0;

    let mut normalized_data1 = HashMap::new();
    let mut normalized_data2 = HashMap::new();

    for (key, values1) in data1 {
        if let Some(values2) = data2.get(key) {
            let norm_values1: Vec<f64> = values1
                .iter()
                .map(|&x| x * target_median / median1)
                .collect();
            let norm_values2: Vec<f64> = values2
                .iter()
                .map(|&x| x * target_median / median2)
                .collect();

            normalized_data1.insert(key.clone(), norm_values1);
            normalized_data2.insert(key.clone(), norm_values2);
        }
    }

    (normalized_data1, normalized_data2)
}

// New: Calculate MAD
fn calculate_mad(values: &[f64], median: f64) -> f64 {
    let mut deviations: Vec<f64> = values.iter().map(|&x| (x - median).abs()).collect();
    deviations.sort_by(|a, b| a.partial_cmp(b).unwrap());

    if deviations.is_empty() {
        1.0
    } else {
        deviations[deviations.len() / 2] * 1.4826 // Convert to standard deviation estimate
    }
}

// New: t-test implementation
fn perform_t_test(group1: &[f64], group2: &[f64], alpha: f64) -> StatisticalResult {
    if group1.len() < 2 || group2.len() < 2 {
        return StatisticalResult {
            test_type: StatisticalTest::TTest,
            statistic: 0.0,
            p_value: 1.0,
            significant: false,
            effect_size: None,
        };
    }

    let n1 = group1.len() as f64;
    let n2 = group2.len() as f64;

    let mean1 = group1.iter().sum::<f64>() / n1;
    let mean2 = group2.iter().sum::<f64>() / n2;

    let var1 = group1.iter().map(|x| (x - mean1).powi(2)).sum::<f64>() / (n1 - 1.0);
    let var2 = group2.iter().map(|x| (x - mean2).powi(2)).sum::<f64>() / (n2 - 1.0);

    let pooled_var = ((n1 - 1.0) * var1 + (n2 - 1.0) * var2) / (n1 + n2 - 2.0);
    let se = (pooled_var * (1.0 / n1 + 1.0 / n2)).sqrt();

    let t_stat = (mean1 - mean2) / se;
    let _df = n1 + n2 - 2.0;

    // Simplified p-value calculation (using normal distribution approximation)
    let p_value = 2.0 * (1.0 - normal_cdf(t_stat.abs()));

    let effect_size = (mean1 - mean2) / pooled_var.sqrt();

    StatisticalResult {
        test_type: StatisticalTest::TTest,
        statistic: t_stat,
        p_value,
        significant: p_value < alpha,
        effect_size: Some(effect_size),
    }
}

// New: Mann-Whitney U test implementation
fn perform_mann_whitney_u_test(group1: &[f64], group2: &[f64], alpha: f64) -> StatisticalResult {
    if group1.is_empty() || group2.is_empty() {
        return StatisticalResult {
            test_type: StatisticalTest::MannWhitneyU,
            statistic: 0.0,
            p_value: 1.0,
            significant: false,
            effect_size: None,
        };
    }

    // Merge data and sort
    let mut all_data = Vec::new();
    for (i, &val) in group1.iter().enumerate() {
        all_data.push((val, 1, i));
    }
    for (i, &val) in group2.iter().enumerate() {
        all_data.push((val, 2, i));
    }
    all_data.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

    // Calculate ranks
    let mut ranks = vec![0.0; all_data.len()];
    let mut current_rank = 1.0;
    let mut i = 0;

    while i < all_data.len() {
        let mut j = i;
        while j < all_data.len() && all_data[j].0 == all_data[i].0 {
            j += 1;
        }

        let avg_rank = (current_rank + (j - 1) as f64) / 2.0;
        for k in i..j {
            ranks[k] = avg_rank;
        }

        current_rank = j as f64 + 1.0;
        i = j;
    }

    // Calculate U statistic
    let mut u1 = 0.0;
    let mut u2 = 0.0;

    for (i, (_, group, _)) in all_data.iter().enumerate() {
        if *group == 1 {
            u1 += ranks[i];
        } else {
            u2 += ranks[i];
        }
    }

    let n1 = group1.len() as f64;
    let n2 = group2.len() as f64;

    u1 = n1 * n2 + n1 * (n1 + 1.0) / 2.0 - u1;
    u2 = n1 * n2 + n2 * (n2 + 1.0) / 2.0 - u2;

    let u_stat = u1.min(u2);

    // Calculate expected value and standard deviation
    let expected_u = n1 * n2 / 2.0;
    let var_u = n1 * n2 * (n1 + n2 + 1.0) / 12.0;
    let se_u = var_u.sqrt();

    // Calculate z statistic
    let z_stat = (u_stat - expected_u) / se_u;
    let p_value = 2.0 * (1.0 - normal_cdf(z_stat.abs()));

    let effect_size = (u_stat - expected_u) / (n1 * n2);

    StatisticalResult {
        test_type: StatisticalTest::MannWhitneyU,
        statistic: u_stat,
        p_value,
        significant: p_value < alpha,
        effect_size: Some(effect_size),
    }
}

// New: Wilcoxon signed-rank test implementation
fn perform_wilcoxon_test(group1: &[f64], group2: &[f64], alpha: f64) -> StatisticalResult {
    if group1.len() != group2.len() || group1.is_empty() {
        return StatisticalResult {
            test_type: StatisticalTest::WilcoxonSignedRank,
            statistic: 0.0,
            p_value: 1.0,
            significant: false,
            effect_size: None,
        };
    }

    // Calculate differences
    let mut differences: Vec<f64> = group1
        .iter()
        .zip(group2.iter())
        .map(|(a, b)| a - b)
        .collect();

    // Remove zero differences
    differences.retain(|&x| x != 0.0);

    if differences.is_empty() {
        return StatisticalResult {
            test_type: StatisticalTest::WilcoxonSignedRank,
            statistic: 0.0,
            p_value: 1.0,
            significant: false,
            effect_size: None,
        };
    }

    // Calculate ranks of absolute differences
    let mut abs_diffs: Vec<(f64, usize)> = differences
        .iter()
        .map(|&x| x.abs())
        .enumerate()
        .map(|(i, x)| (x, i))
        .collect();
    abs_diffs.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

    let mut ranks = vec![0.0; abs_diffs.len()];
    let mut current_rank = 1.0;
    let mut i = 0;

    while i < abs_diffs.len() {
        let mut j = i;
        while j < abs_diffs.len() && abs_diffs[j].0 == abs_diffs[i].0 {
            j += 1;
        }

        let avg_rank = (current_rank + (j - 1) as f64) / 2.0;
        for k in i..j {
            ranks[k] = avg_rank;
        }

        current_rank = j as f64 + 1.0;
        i = j;
    }

    // Calculate W statistic
    let mut w_plus = 0.0;
    let mut w_minus = 0.0;

    for (i, &diff) in differences.iter().enumerate() {
        if diff > 0.0 {
            w_plus += ranks[i];
        } else {
            w_minus += ranks[i];
        }
    }

    let w_stat = w_plus.min(w_minus);
    let n = differences.len() as f64;

    // Calculate expected value and standard deviation
    let expected_w = n * (n + 1.0) / 4.0;
    let var_w = n * (n + 1.0) * (2.0 * n + 1.0) / 24.0;
    let se_w = var_w.sqrt();

    // Calculate z statistic
    let z_stat = (w_stat - expected_w) / se_w;
    let p_value = 2.0 * (1.0 - normal_cdf(z_stat.abs()));

    let effect_size = (w_stat - expected_w) / (n * (n + 1.0) / 2.0);

    StatisticalResult {
        test_type: StatisticalTest::WilcoxonSignedRank,
        statistic: w_stat,
        p_value,
        significant: p_value < alpha,
        effect_size: Some(effect_size),
    }
}

// New: Chi-square test implementation
fn perform_chi_square_test(group1: &[f64], group2: &[f64], alpha: f64) -> StatisticalResult {
    // Simplified chi-square test, discretize continuous data
    let threshold = (group1.iter().sum::<f64>() + group2.iter().sum::<f64>())
        / (group1.len() + group2.len()) as f64;

    let mut observed = [[0; 2]; 2]; // [group][above_threshold]

    for &val in group1 {
        if val > threshold {
            observed[0][1] += 1;
        } else {
            observed[0][0] += 1;
        }
    }

    for &val in group2 {
        if val > threshold {
            observed[1][1] += 1;
        } else {
            observed[1][0] += 1;
        }
    }

    let total = observed
        .iter()
        .map(|row| row.iter().sum::<i32>())
        .sum::<i32>();
    if total == 0 {
        return StatisticalResult {
            test_type: StatisticalTest::ChiSquare,
            statistic: 0.0,
            p_value: 1.0,
            significant: false,
            effect_size: None,
        };
    }

    let mut chi_square = 0.0;
    for i in 0..2 {
        for j in 0..2 {
            let row_sum: i32 = observed[i].iter().sum();
            let col_sum: i32 = observed.iter().map(|row| row[j]).sum();
            let expected = (row_sum * col_sum) as f64 / total as f64;

            if expected > 0.0 {
                chi_square += (observed[i][j] as f64 - expected).powi(2) / expected;
            }
        }
    }

    // Simplified p-value calculation (using chi-square distribution approximation)
    let p_value = 1.0 - chi_square_cdf(chi_square, 1);

    StatisticalResult {
        test_type: StatisticalTest::ChiSquare,
        statistic: chi_square,
        p_value,
        significant: p_value < alpha,
        effect_size: None,
    }
}

// New: Continuity correction test implementation
fn perform_continuity_correction_test(
    group1: &[f64],
    group2: &[f64],
    alpha: f64,
) -> StatisticalResult {
    // Continuity correction based on chi-square test
    let result = perform_chi_square_test(group1, group2, alpha);

    // Apply Yates continuity correction
    let corrected_statistic = (result.statistic - 0.5).max(0.0);
    let corrected_p_value = 1.0 - chi_square_cdf(corrected_statistic, 1);

    StatisticalResult {
        test_type: StatisticalTest::ContinuityCorrection,
        statistic: corrected_statistic,
        p_value: corrected_p_value,
        significant: corrected_p_value < alpha,
        effect_size: result.effect_size,
    }
}

// New: Normal distribution cumulative distribution function approximation
fn normal_cdf(x: f64) -> f64 {
    0.5 * (1.0 + erf(x / 2.0_f64.sqrt()))
}

// New: Error function approximation
fn erf(x: f64) -> f64 {
    let a1 = 0.254829592;
    let a2 = -0.284496736;
    let a3 = 1.421413741;
    let a4 = -1.453152027;
    let a5 = 1.061405429;
    let p = 0.3275911;

    let sign = if x < 0.0 { -1.0 } else { 1.0 };
    let x = x.abs();

    let t = 1.0 / (1.0 + p * x);
    let y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * (-x * x).exp();

    sign * y
}

// New: Chi-square distribution cumulative distribution function approximation
fn chi_square_cdf(x: f64, df: i32) -> f64 {
    // Simplified chi-square distribution CDF approximation
    if x <= 0.0 {
        return 0.0;
    }

    // Use normal distribution approximation
    let mu = df as f64;
    let sigma = (2.0 * df as f64).sqrt();
    let z = (x - mu) / sigma;

    normal_cdf(z)
}

fn load_csv_data(file_path: &str) -> Result<HashMap<String, PositionData>, Box<dyn Error>> {
    let mut data = HashMap::new();
    let file = File::open(file_path)?;
    let reader = BufReader::new(file);

    for (line_num, line) in reader.lines().enumerate() {
        if line_num == 0 {
            continue; // Skip header row
        }

        let line = line?;
        let fields: Vec<&str> = line.split(',').collect();

        if fields.len() >= 13 {
            let chr_id = fields[0].trim_matches('"').to_string();
            let strand = fields[1].chars().next().unwrap_or('+');
            let position = fields[2].parse::<u32>()?;
            let ref_base = fields[3].chars().next().unwrap_or('N');
            let rfc = fields[4].parse::<usize>()?;
            let pipc = fields[5].parse::<usize>()?;
            let depth = fields[6].parse::<usize>()?;
            let ins = fields[7].parse::<usize>()?;
            let del = fields[8].parse::<usize>()?;
            let a = fields[9].parse::<usize>()?;
            let c = fields[10].parse::<usize>()?;
            let g = fields[11].parse::<usize>()?;
            // Fix: handle floating point format in last column
            let t = if fields[12].contains('.') {
                // If floating point format, parse as f64 first then convert to u32
                fields[12].parse::<f64>().unwrap_or(0.0) as u32
            } else {
                // If integer format, parse directly
                fields[12].parse::<u32>()?
            };

            // Use chr_id+strand+position as key, ensure positive and negative strands are processed separately
            let key = format!("{}:{}{}", chr_id, strand, position);
            data.insert(
                key,
                PositionData {
                    chr_id,
                    strand,
                    position,
                    ref_base,
                    rfc,
                    pipc,
                    depth,
                    ins,
                    del,
                    a: a as u32,
                    c: c as u32,
                    g: g as u32,
                    t,
                },
            );
        }
    }

    Ok(data)
}

// New: Multiple replicate comparison function
pub fn compare_biological_replicates(
    group1_files: &[String],
    group2_files: &[String],
    output_file: &str,
    test_type: StatisticalTest,
    alpha: f64,
    min_depth: usize,
    logger: &mut crate::Logger,
) -> Result<(), Box<dyn Error>> {
    let start_time = Instant::now();

    logger.log("=== ModDetector Biological Replicates Comparison ===")?;
    logger.log(&format!("Group 1 files: {:?}", group1_files))?;
    logger.log(&format!("Group 2 files: {:?}", group2_files))?;
    logger.log(&format!("Statistical test: {:?}", test_type))?;
    logger.log(&format!("Alpha: {}", alpha))?;

    println!("[Loading data]");
    println!("    Group 1: {} files", group1_files.len());
    println!("    Group 2: {} files", group2_files.len());
    println!();

    // Load all replicate data
    let mut group1_data = HashMap::new();
    let mut group2_data = HashMap::new();

    for file in group1_files {
        let data = load_csv_data(file)?;
        for (key, pos_data) in data {
            group1_data
                .entry(key)
                .or_insert_with(Vec::new)
                .push(pos_data);
        }
    }

    for file in group2_files {
        let data = load_csv_data(file)?;
        for (key, pos_data) in data {
            group2_data
                .entry(key)
                .or_insert_with(Vec::new)
                .push(pos_data);
        }
    }

    // Find shared positions
    let shared_positions: Vec<String> = group1_data
        .keys()
        .filter(|k| group2_data.contains_key(*k))
        .cloned()
        .collect();

    println!("[Data info]");
    println!("    Group 1 positions: {}", group1_data.len());
    println!("    Group 2 positions: {}", group2_data.len());
    println!("    Shared positions: {}", shared_positions.len());
    println!();

    // Create output file
    let mut output = File::create(output_file)?;
    writeln!(output, "ChrID,Strand,Position,RefBase,Group1_Mean_RFC,Group1_Mean_PIPC,Group1_Mean_Depth,Group2_Mean_RFC,Group2_Mean_PIPC,Group2_Mean_Depth,RFC_FoldChange,PIPC_FoldChange,Test_Statistic,P_Value,Significant,Effect_Size,Modification_Score")?;

    let mut significant_sites = 0;
    let mut processed_count = 0;

    for key in &shared_positions {
        processed_count += 1;
        if processed_count % 10000 == 0 {
            println!(
                "\r[Progressing] Processing {}/{} positions...",
                processed_count,
                shared_positions.len()
            );
            std::io::stdout().flush()?;
        }

        let group1_replicates = &group1_data[key];
        let group2_replicates = &group2_data[key];

        // Calculate average depth
        let group1_avg_depth: f64 = group1_replicates
            .iter()
            .map(|x| x.depth as f64)
            .sum::<f64>()
            / group1_replicates.len() as f64;
        let group2_avg_depth: f64 = group2_replicates
            .iter()
            .map(|x| x.depth as f64)
            .sum::<f64>()
            / group2_replicates.len() as f64;

        // Filter low depth positions
        if group1_avg_depth < min_depth as f64 || group2_avg_depth < min_depth as f64 {
            continue;
        }

        // Extract RFC and PIPC values
        let group1_rfc: Vec<f64> = group1_replicates.iter().map(|x| x.rfc as f64).collect();
        let group1_pipc: Vec<f64> = group1_replicates.iter().map(|x| x.pipc as f64).collect();
        let group2_rfc: Vec<f64> = group2_replicates.iter().map(|x| x.rfc as f64).collect();
        let group2_pipc: Vec<f64> = group2_replicates.iter().map(|x| x.pipc as f64).collect();

        // Perform statistical test
        let rfc_test = perform_statistical_test(&group1_rfc, &group2_rfc, &test_type, alpha);
        let pipc_test = perform_statistical_test(&group1_pipc, &group2_pipc, &test_type, alpha);

        // Calculate fold change
        let group1_mean_rfc = group1_rfc.iter().sum::<f64>() / group1_rfc.len() as f64;
        let group2_mean_rfc = group2_rfc.iter().sum::<f64>() / group2_rfc.len() as f64;
        let group1_mean_pipc = group1_pipc.iter().sum::<f64>() / group1_pipc.len() as f64;
        let group2_mean_pipc = group2_pipc.iter().sum::<f64>() / group2_pipc.len() as f64;

        let rfc_fold_change = if group2_mean_rfc > 0.0 {
            group1_mean_rfc / group2_mean_rfc
        } else if group1_mean_rfc > 0.0 {
            f64::INFINITY
        } else {
            1.0
        };

        let pipc_fold_change = if group2_mean_pipc > 0.0 {
            group1_mean_pipc / group2_mean_pipc
        } else if group1_mean_pipc > 0.0 {
            f64::INFINITY
        } else {
            1.0
        };

        // Calculate modification score
        let modification_score = (rfc_fold_change + pipc_fold_change) / 2.0;

        // Check significance
        let is_significant = rfc_test.significant || pipc_test.significant;
        if is_significant {
            significant_sites += 1;
        }

        // Parse position information
        let parts: Vec<&str> = key.split(':').collect();
        if parts.len() >= 2 {
            let chr_id = parts[0];
            let strand_pos = parts[1];
            let strand = strand_pos.chars().next().unwrap_or('+');
            let position = strand_pos[1..].parse::<u32>().unwrap_or(0);
            let ref_base = group1_replicates[0].ref_base;

            writeln!(output, "{},{},{},{},{:.3},{:.3},{:.3},{:.3},{:.3},{:.3},{:.3},{:.3},{:.6},{:.6},{},{:.3},{:.3}",
                chr_id,
                strand,
                position,
                ref_base,
                group1_mean_rfc,
                group1_mean_pipc,
                group1_avg_depth,
                group2_mean_rfc,
                group2_mean_pipc,
                group2_avg_depth,
                rfc_fold_change,
                pipc_fold_change,
                rfc_test.statistic,
                rfc_test.p_value,
                if is_significant { "TRUE" } else { "FALSE" },
                rfc_test.effect_size.unwrap_or(0.0),
                modification_score
            )?;
        }
    }

    let elapsed = start_time.elapsed();

    println!("\r[Output]");
    println!("    Output: {}", output_file);
    println!("    Significant positions: {}", significant_sites);
    println!("    Total processed: {}", processed_count);
    println!("{}", crate::progress::format_time_used(elapsed));

    logger.log(&format!("Biological replicates comparison completed"))?;
    logger.log(&format!("Significant sites: {}", significant_sites))?;
    logger.log(&format!("Total time: {:.2}s", elapsed.as_secs_f64()))?;

    Ok(())
}

// New: Compare two reactivity result files function
pub fn compare_reactivity_results(
    reactivity1_file: &str,
    reactivity2_file: &str,
    output_file: &str,
    test_type: StatisticalTest,
    alpha: f64,
    logger: &mut crate::Logger,
) -> Result<(), Box<dyn Error>> {
    let start_time = Instant::now();

    logger.log("=== ModDetector Reactivity Results Comparison ===")?;
    logger.log(&format!("Reactivity 1 file: {}", reactivity1_file))?;
    logger.log(&format!("Reactivity 2 file: {}", reactivity2_file))?;
    logger.log(&format!("Statistical test: {:?}", test_type))?;
    logger.log(&format!("Alpha: {}", alpha))?;

    println!("[Loading data]");
    println!("    Reactivity 1: {}", reactivity1_file);
    println!("    Reactivity 2: {}", reactivity2_file);
    println!();

    // Load reactivity data
    let reactivity1_data = load_reactivity_data(reactivity1_file)?;
    let reactivity2_data = load_reactivity_data(reactivity2_file)?;

    // Find shared positions
    let shared_positions: Vec<String> = reactivity1_data
        .keys()
        .filter(|k| reactivity2_data.contains_key(*k))
        .cloned()
        .collect();

    println!("[Data info]");
    println!("    Reactivity 1 positions: {}", reactivity1_data.len());
    println!("    Reactivity 2 positions: {}", reactivity2_data.len());
    println!("    Shared positions: {}", shared_positions.len());
    println!();

    // Create output file (point-by-point results)
    let mut output = File::create(output_file)?;
    writeln!(output, "ChrID,Strand,Position,Reactivity1_Mean,Reactivity1_Std,Reactivity2_Mean,Reactivity2_Std,FoldChange,Test_Statistic,P_Value,Significant,Effect_Size")?;
    // Also create region result file (based on Q statistic SVR)
    let region_output_file = format!("{}{}", output_file, ".regions.csv");
    let mut region_out = File::create(&region_output_file)?;
    writeln!(
        region_out,
        "ChrID,Strand,Start,End,Length,MaxQ,MonteCarlo_Threshold,Significant"
    )?;

    let mut significant_sites = 0;
    let mut processed_count = 0;

    for key in &shared_positions {
        processed_count += 1;
        if processed_count % 10000 == 0 {
            println!(
                "\r[Progressing] Processing {}/{} positions...",
                processed_count,
                shared_positions.len()
            );
            std::io::stdout().flush()?;
        }

        let reactivity1_values = &reactivity1_data[key];
        let reactivity2_values = &reactivity2_data[key];

        // Calculate statistics
        let reactivity1_mean =
            reactivity1_values.iter().sum::<f64>() / reactivity1_values.len() as f64;
        let reactivity2_mean =
            reactivity2_values.iter().sum::<f64>() / reactivity2_values.len() as f64;

        let reactivity1_std = if reactivity1_values.len() > 1 {
            let reactivity1_var = reactivity1_values
                .iter()
                .map(|x| (x - reactivity1_mean).powi(2))
                .sum::<f64>()
                / (reactivity1_values.len() - 1) as f64;
            reactivity1_var.sqrt()
        } else {
            0.0 // Single value case, standard deviation is 0
        };

        let reactivity2_std = if reactivity2_values.len() > 1 {
            let reactivity2_var = reactivity2_values
                .iter()
                .map(|x| (x - reactivity2_mean).powi(2))
                .sum::<f64>()
                / (reactivity2_values.len() - 1) as f64;
            reactivity2_var.sqrt()
        } else {
            0.0 // Single value case, standard deviation is 0
        };

        let fold_change = if reactivity2_mean > 0.0 {
            reactivity1_mean / reactivity2_mean
        } else if reactivity1_mean > 0.0 {
            f64::INFINITY
        } else {
            1.0
        };

        // Perform statistical test
        let test_result =
            perform_statistical_test(reactivity1_values, reactivity2_values, &test_type, alpha);

        if test_result.significant {
            significant_sites += 1;
        }

        // Parse position information
        let parts: Vec<&str> = key.split(':').collect();
        if parts.len() >= 2 {
            let chr_id = parts[0];
            let strand_pos = parts[1];
            let strand = strand_pos.chars().next().unwrap_or('+');
            let position = strand_pos[1..].parse::<u32>().unwrap_or(0);

            writeln!(
                output,
                "{},{},{},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{},{}",
                chr_id,
                strand,
                position,
                reactivity1_mean,
                reactivity1_std,
                reactivity2_mean,
                reactivity2_std,
                fold_change,
                test_result.statistic,
                test_result.p_value,
                if test_result.significant {
                    "TRUE"
                } else {
                    "FALSE"
                },
                test_result.effect_size.unwrap_or(0.0)
            )?;
        }
    }

    // Perform window p-value and scan statistics separately for each transcript
    // Group shared positions by transcript (chr_id+strand) and sort by position
    let mut by_tx: HashMap<(String, char), Vec<(u32, f64, f64)>> = HashMap::new();
    for key in &shared_positions {
        let parts: Vec<&str> = key.split(':').collect();
        if parts.len() >= 2 {
            let chr_id = parts[0].to_string();
            let strand_pos = parts[1];
            let strand = strand_pos.chars().next().unwrap_or('+');
            let position = strand_pos[1..].parse::<u32>().unwrap_or(0);
            let v1 = reactivity1_data[key][0];
            let v2 = reactivity2_data[key][0];
            by_tx
                .entry((chr_id, strand))
                .or_insert_with(Vec::new)
                .push((position, v1, v2));
        }
    }
    let window_radius = 2usize; // Window radius, can be parameterized
    let min_len = 5usize; // Minimum region length
    let max_len = 50usize; // Maximum region length
    let mc_iter = 200usize; // Monte Carlo iteration count (can be parameterized)

    for ((chr_id, strand), mut rows) in by_tx {
        rows.sort_by_key(|r| r.0);
        let series1: Vec<(u32, f64)> = rows.iter().map(|(p, v1, _)| (*p, *v1)).collect();
        let series2: Vec<(u32, f64)> = rows.iter().map(|(p, _, v2)| (*p, *v2)).collect();

        match test_type {
            StatisticalTest::DeltaSHAPE => {
                // deltaSHAPE region detection: use sliding window method with Z-factor and Z-score
                let deltashape_regions =
                    compute_deltashape_regions(&series1, &series2, window_radius, min_len, max_len);
                for (start, end, z_factor, z_score, significant) in deltashape_regions {
                    if significant {
                        let len = if end >= start {
                            (end - start + 1) as usize
                        } else {
                            0
                        };
                        writeln!(
                            region_out,
                            "{},{},{},{},{},{:.6},{:.6},{}",
                            chr_id, strand, start, end, len, z_factor, z_score, "TRUE"
                        )?;
                    }
                }
            }
            _ => {
                // Other methods use DiffScan region detection logic
                let pos_pvals = compute_window_pvalues(&series1, &series2, window_radius);
                if pos_pvals.is_empty() {
                    continue;
                }
                let regions = compute_scan_Q(&pos_pvals, min_len, max_len);
                if regions.is_empty() {
                    continue;
                }
                let q_thresh =
                    monte_carlo_q_threshold(&pos_pvals, min_len, max_len, alpha, mc_iter);
                // Output regions exceeding threshold
                for (start, end, q) in regions.into_iter() {
                    let len = if end >= start {
                        (end - start + 1) as usize
                    } else {
                        0
                    };
                    let sig = q >= q_thresh;
                    if sig {
                        writeln!(
                            region_out,
                            "{},{},{},{},{},{:.6},{:.6},{}",
                            chr_id,
                            strand,
                            start,
                            end,
                            len,
                            q,
                            q_thresh,
                            if sig { "TRUE" } else { "FALSE" }
                        )?;
                    }
                }
            }
        }
    }

    let elapsed = start_time.elapsed();

    println!("\r[Output]");
    println!("    Output: {}", output_file);
    println!("    Significant positions: {}", significant_sites);
    println!("    Total processed: {}", processed_count);
    println!("{}", crate::progress::format_time_used(elapsed));

    logger.log(&format!("Reactivity results comparison completed"))?;
    logger.log(&format!("Significant sites: {}", significant_sites))?;
    logger.log(&format!("Region output: {}", region_output_file))?;
    logger.log(&format!("Total time: {:.2}s", elapsed.as_secs_f64()))?;

    Ok(())
}

// New: Reactivity group comparison function
pub fn compare_reactivity_groups(
    group1_files: &[String],
    group2_files: &[String],
    output_file: &str,
    test_type: StatisticalTest,
    alpha: f64,
    logger: &mut crate::Logger,
) -> Result<(), Box<dyn Error>> {
    let start_time = Instant::now();

    logger.log("=== ModDetector Reactivity Groups Comparison ===")?;
    logger.log(&format!("Group 1 files: {:?}", group1_files))?;
    logger.log(&format!("Group 2 files: {:?}", group2_files))?;
    logger.log(&format!("Statistical test: {:?}", test_type))?;
    logger.log(&format!("Alpha: {}", alpha))?;

    println!("[Loading data]");
    println!("    Group 1: {} files", group1_files.len());
    println!("    Group 2: {} files", group2_files.len());
    println!();

    // Load reactivity data
    let mut group1_data = HashMap::new();
    let mut group2_data = HashMap::new();

    for file in group1_files {
        let data = load_reactivity_data(file)?;
        for (key, values) in data {
            group1_data
                .entry(key)
                .or_insert_with(Vec::new)
                .extend(values);
        }
    }

    for file in group2_files {
        let data = load_reactivity_data(file)?;
        for (key, values) in data {
            group2_data
                .entry(key)
                .or_insert_with(Vec::new)
                .extend(values);
        }
    }

    // Find shared positions
    let shared_positions: Vec<String> = group1_data
        .keys()
        .filter(|k| group2_data.contains_key(*k))
        .cloned()
        .collect();

    println!("[Data info]");
    println!("    Group 1 positions: {}", group1_data.len());
    println!("    Group 2 positions: {}", group2_data.len());
    println!("    Shared positions: {}", shared_positions.len());
    println!();

    // Create output file
    let mut output = File::create(output_file)?;
    writeln!(output, "ChrID,Strand,Position,Group1_Mean_Reactivity,Group1_Std_Reactivity,Group2_Mean_Reactivity,Group2_Std_Reactivity,FoldChange,Test_Statistic,P_Value,Significant,Effect_Size")?;

    let mut significant_sites = 0;
    let mut processed_count = 0;

    for key in &shared_positions {
        processed_count += 1;
        if processed_count % 10000 == 0 {
            println!(
                "\r[Progressing] Processing {}/{} positions...",
                processed_count,
                shared_positions.len()
            );
            std::io::stdout().flush()?;
        }

        let group1_values = &group1_data[key];
        let group2_values = &group2_data[key];

        // Perform statistical test
        let test_result = perform_statistical_test(group1_values, group2_values, &test_type, alpha);

        // Calculate statistics
        let group1_mean = group1_values.iter().sum::<f64>() / group1_values.len() as f64;
        let group2_mean = group2_values.iter().sum::<f64>() / group2_values.len() as f64;

        let group1_var = group1_values
            .iter()
            .map(|x| (x - group1_mean).powi(2))
            .sum::<f64>()
            / (group1_values.len() - 1) as f64;
        let group2_var = group2_values
            .iter()
            .map(|x| (x - group2_mean).powi(2))
            .sum::<f64>()
            / (group2_values.len() - 1) as f64;

        let group1_std = group1_var.sqrt();
        let group2_std = group2_var.sqrt();

        let fold_change = if group2_mean > 0.0 {
            group1_mean / group2_mean
        } else if group1_mean > 0.0 {
            f64::INFINITY
        } else {
            1.0
        };

        if test_result.significant {
            significant_sites += 1;
        }

        // Parse position information
        let parts: Vec<&str> = key.split(':').collect();
        if parts.len() >= 2 {
            let chr_id = parts[0];
            let strand_pos = parts[1];
            let strand = strand_pos.chars().next().unwrap_or('+');
            let position = strand_pos[1..].parse::<u32>().unwrap_or(0);

            writeln!(
                output,
                "{},{},{},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{},{}",
                chr_id,
                strand,
                position,
                group1_mean,
                group1_std,
                group2_mean,
                group2_std,
                fold_change,
                test_result.statistic,
                test_result.p_value,
                if test_result.significant {
                    "TRUE"
                } else {
                    "FALSE"
                },
                test_result.effect_size.unwrap_or(0.0)
            )?;
        }
    }

    let elapsed = start_time.elapsed();

    println!("\r[Output]");
    println!("    Output: {}", output_file);
    println!("    Significant positions: {}", significant_sites);
    println!("    Total processed: {}", processed_count);
    println!("{}", crate::progress::format_time_used(elapsed));

    logger.log(&format!("Reactivity groups comparison completed"))?;
    logger.log(&format!("Significant sites: {}", significant_sites))?;
    logger.log(&format!("Total time: {:.2}s", elapsed.as_secs_f64()))?;

    Ok(())
}

    // Modify original compare_samples function, add statistical test support
pub fn compare_samples_with_stats(
    mod_file: &str,
    unmod_file: &str,
    output_file: &str,
    min_depth: usize,
    min_fold_change: f64,
    test_type: Option<StatisticalTest>,
    alpha: f64,
    logger: &mut crate::Logger,
) -> Result<(), Box<dyn Error>> {
    let start_time = Instant::now();

    logger.log("=== ModDetector Enhanced Compare Function ===")?;
    logger.log(&format!("Modified Sample File: {}", mod_file))?;
    logger.log(&format!("Unmodified Sample File: {}", unmod_file))?;
    logger.log(&format!("Output File: {}", output_file))?;
    logger.log(&format!("Minimum Depth: {}", min_depth))?;
    logger.log(&format!("Minimum Fold Change: {}", min_fold_change))?;
    if let Some(test) = &test_type {
        logger.log(&format!("Statistical Test: {:?}", test))?;
        logger.log(&format!("Alpha: {}", alpha))?;
    }

    println!("[Loading data]");
    println!("    Mod: {}", mod_file);
    println!("    Unmod: {}", unmod_file);
    println!();

    println!("[Params]");
    println!(
        "    Min depth: {}. Min fold change: {}.",
        min_depth, min_fold_change
    );
    if let Some(test) = &test_type {
        println!("    Statistical test: {:?}. Alpha: {}.", test, alpha);
    }
    println!();

    let mod_data = load_csv_data(mod_file)?;
    let unmod_data = load_csv_data(unmod_file)?;

    let shared_positions = mod_data
        .keys()
        .filter(|k| unmod_data.contains_key(*k))
        .count();

    println!("[Data info]");
    println!(
        "    mod position: {}. unmod position: {}.",
        mod_data.len(),
        unmod_data.len()
    );
    println!("    shared position: {}.", shared_positions);
    println!();

    let mut output = File::create(output_file)?;

    // Choose output format based on whether statistical test is present
    if let Some(_) = test_type {
        writeln!(output, "ChrID,Strand,Position,RefBase,Mod_RFC,Mod_PIPC,Mod_Depth,Unmod_RFC,Unmod_PIPC,Unmod_Depth,RFC_FoldChange,PIPC_FoldChange,Test_Statistic,P_Value,Significant,Effect_Size,Modification_Score")?;
    } else {
        writeln!(output, "ChrID,Strand,Position,RefBase,Mod_RFC,Mod_PIPC,Mod_Depth,Unmod_RFC,Unmod_PIPC,Unmod_Depth,RFC_FoldChange,PIPC_FoldChange,Modification_Score")?;
    }

    let mut detected_sites = 0;
    let mut processed_count = 0;

    for (key, mod_pos) in &mod_data {
        if let Some(unmod_pos) = unmod_data.get(key) {
            processed_count += 1;
            if processed_count % 10000 == 0 {
                println!(
                    "\r[Progressing] Comparing {}/{} shared positions...",
                    processed_count, shared_positions
                );
                std::io::stdout().flush()?;
            }

            if mod_pos.depth < min_depth || unmod_pos.depth < min_depth {
                continue;
            }

            let rfc_fold_change = if unmod_pos.rfc > 0 {
                mod_pos.rfc as f64 / unmod_pos.rfc as f64
            } else if mod_pos.rfc > 0 {
                f64::INFINITY
            } else {
                1.0
            };

            let pipc_fold_change = if unmod_pos.pipc > 0 {
                mod_pos.pipc as f64 / unmod_pos.pipc as f64
            } else if mod_pos.pipc > 0 {
                f64::INFINITY
            } else {
                1.0
            };

            let modification_score = (rfc_fold_change + pipc_fold_change) / 2.0;

            let mut should_output =
                rfc_fold_change >= min_fold_change || pipc_fold_change >= min_fold_change;
            let mut test_result = None;

            // Perform statistical test
            if let Some(test_type) = &test_type {
                let group1_rfc = vec![mod_pos.rfc as f64];
                let group2_rfc = vec![unmod_pos.rfc as f64];
                let group1_pipc = vec![mod_pos.pipc as f64];
                let group2_pipc = vec![unmod_pos.pipc as f64];

                let rfc_test = perform_statistical_test(&group1_rfc, &group2_rfc, test_type, alpha);
                let pipc_test =
                    perform_statistical_test(&group1_pipc, &group2_pipc, test_type, alpha);

                // Use stricter significance criteria
                should_output = should_output && (rfc_test.significant || pipc_test.significant);
                test_result = Some(rfc_test);
            }

            if should_output {
                if let Some(test_result) = test_result {
                    writeln!(
                        output,
                        "{},{},{},{},{},{},{},{},{},{},{:.3},{:.3},{:.6},{:.6},{},{:.3},{:.3}",
                        mod_pos.chr_id,
                        mod_pos.strand,
                        mod_pos.position,
                        mod_pos.ref_base,
                        mod_pos.rfc,
                        mod_pos.pipc,
                        mod_pos.depth,
                        unmod_pos.rfc,
                        unmod_pos.pipc,
                        unmod_pos.depth,
                        rfc_fold_change,
                        pipc_fold_change,
                        test_result.statistic,
                        test_result.p_value,
                        if test_result.significant {
                            "TRUE"
                        } else {
                            "FALSE"
                        },
                        test_result.effect_size.unwrap_or(0.0),
                        modification_score
                    )?;
                } else {
                    writeln!(
                        output,
                        "{},{},{},{},{},{},{},{},{},{},{:.3},{:.3},{:.3}",
                        mod_pos.chr_id,
                        mod_pos.strand,
                        mod_pos.position,
                        mod_pos.ref_base,
                        mod_pos.rfc,
                        mod_pos.pipc,
                        mod_pos.depth,
                        unmod_pos.rfc,
                        unmod_pos.pipc,
                        unmod_pos.depth,
                        rfc_fold_change,
                        pipc_fold_change,
                        modification_score
                    )?;
                }
                detected_sites += 1;
            }
        }
    }

    let elapsed = start_time.elapsed();

    println!("\r[Output]");
    println!("    Shared: {}", output_file);
    println!("    Significant modified position: {}.", detected_sites);
    println!("{}", crate::progress::format_time_used(elapsed));

    logger.log(&format!("Enhanced sample comparison completed"))?;
    logger.log(&format!("Detected modification sites: {}", detected_sites))?;
    logger.log(&format!("Total time: {:.2}s", elapsed.as_secs_f64()))?;

    Ok(())
}

// Keep original compare_samples function for backward compatibility
pub fn compare_samples(
    mod_file: &str,
    unmod_file: &str,
    output_file: &str,
    min_depth: usize,
    min_fold_change: f64,
    logger: &mut crate::Logger,
) -> Result<(), Box<dyn Error>> {
    compare_samples_with_stats(
        mod_file,
        unmod_file,
        output_file,
        min_depth,
        min_fold_change,
        None,
        0.05,
        logger,
    )
}
