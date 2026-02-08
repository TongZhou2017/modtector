// PCR bias correction module
// Applies Chi-Square distribution-based depth correction to pileup CSV files

use std::collections::HashMap;
use std::error::Error;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use clap::Args;

/// Command-line arguments for PCR bias correction
#[derive(Args, Debug)]
pub struct CorrectArgs {
    /// Input pileup CSV file
    #[arg(short = 'i', long = "input")]
    pub input: String,
    /// Output pileup CSV file (with corrected effective_depth)
    #[arg(short = 'o', long = "output")]
    pub output: String,
    /// Weight for increasing depth correction (default: 1.0)
    #[arg(long = "weight-increase", default_value_t = 1.0)]
    pub weight_increase: f64,
    /// Weight for decreasing depth correction (default: 0.5)
    #[arg(long = "weight-decrease", default_value_t = 0.5)]
    pub weight_decrease: f64,
    /// Species name for filtering (optional)
    #[arg(long = "species")]
    pub species: Option<String>,
}

/// Chi-Square distribution model (df=2)
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
#[inline]
fn chi2_pdf_df2(x_over_scale: f64) -> f64 {
    (-x_over_scale / 2.0).exp() * 0.5
}

/// Given scale, compute closed-form optimal a and SSE
fn sse_with_closed_form_a(
    binned: &[(f64, f64)],
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

/// Brent's method for 1D minimization
fn brent_minimize(
    f: &dyn Fn(f64) -> f64,
    mut a: f64,
    mut b: f64,
    mut c: f64,
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
fn minimize_scale_logspace(
    binned: &[(f64, f64)],
    scale_lb: f64,
    scale_ub: f64,
    a_lb: f64,
    a_ub: f64,
) -> (f64, f64, f64) {
    let log_lb = scale_lb.ln();
    let log_ub = scale_ub.ln();

    // 1) Coarse grid search
    let grid_n = 120;
    let mut best_i = 0usize;
    let mut best_val = f64::INFINITY;
    let mut best_scale_grid = 0.0;

    for i in 0..grid_n {
        let t = i as f64 / (grid_n as f64 - 1.0);
        let s = log_lb + t * (log_ub - log_lb);
        let scale = s.exp();
        let (sse, _a) = sse_with_closed_form_a(binned, scale, a_lb, a_ub);
        if sse < best_val {
            best_val = sse;
            best_i = i;
            best_scale_grid = scale;
        }
        if i % 20 == 0 || i == grid_n - 1 {
            eprintln!("[DEBUG] Grid i={}: scale={:.2}, sse={:.6e}", i, scale, sse);
        }
    }
    
    eprintln!("[DEBUG] Grid search: best_i={}, best_scale={:.2}, best_sse={:.6e}", best_i, best_scale_grid, best_val);

    // 2) Construct bracket
    let i0 = best_i.saturating_sub(1);
    let i1 = best_i;
    let i2 = (best_i + 1).min(grid_n - 1);

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
    
    let (sse_a, _) = sse_with_closed_form_a(binned, s_a.exp(), a_lb, a_ub);
    let (sse_b, _) = sse_with_closed_form_a(binned, s_b.exp(), a_lb, a_ub);
    let (sse_c, _) = sse_with_closed_form_a(binned, s_c.exp(), a_lb, a_ub);
    eprintln!("[DEBUG] Brent bracket SSE: f(s_a)={:.6e}, f(s_b)={:.6e}, f(s_c)={:.6e}", sse_a, sse_b, sse_c);
    
    let (s_best, sse_best) = brent_minimize(&f, s_a, s_b, s_c, 1e-10, 200);
    let scale_best = s_best.exp();
    let (sse_final, a_final) = sse_with_closed_form_a(binned, scale_best, a_lb, a_ub);
    
    eprintln!("[DEBUG] Brent result: s_best={:.6e}, scale_best={:.2}, sse_best={:.6e}, sse_final={:.6e}, a_final={:.6e}", 
              s_best, scale_best, sse_best, sse_final, a_final);
    
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

/// Fit Chi-Square distribution to data
fn fit_chi_square_to_data(
    depths: &[f64],
    mutation_rates: &[f64],
) -> Option<(f64, f64, f64)> {
    if depths.len() < 10 {
        return None;
    }
    
    // Bin data by depth (matching Python's pd.cut)
    let mut data: Vec<(f64, f64)> = depths.iter().zip(mutation_rates.iter())
        .map(|(d, mr)| (*d, *mr))
        .collect();
    data.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
    
    if data.is_empty() {
        return None;
    }
    
    let min_depth = data[0].0;
    let max_depth = data[data.len() - 1].0;
    let n_bins = 50;
    
    if (max_depth - min_depth) < 1e-6 {
        return None;
    }
    
    // Create value-based bins (matching pd.cut)
    let bin_width = (max_depth - min_depth) / n_bins as f64;
    let mut binned: Vec<(f64, f64)> = Vec::new();
    
    for i in 0..n_bins {
        let bin_left = if i == 0 {
            min_depth - 1e-6
        } else {
            min_depth + (i as f64 * bin_width)
        };
        let bin_right = if i == n_bins - 1 {
            max_depth + 1e-6
        } else {
            min_depth + ((i + 1) as f64 * bin_width)
        };
        
        // Collect points in this bin
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
    let median_depth = {
        let mut depths: Vec<f64> = binned.iter().map(|(d, _)| *d).collect();
        depths.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        let n = depths.len();
        if n == 0 {
            0.0
        } else if n % 2 == 0 {
            (depths[n / 2 - 1] + depths[n / 2]) / 2.0
        } else {
            depths[n / 2]
        }
    };
    
    let bounds_a = (1e-6, 100.0);
    let bounds_scale = (1e-3, median_depth * 0.2);
    
    eprintln!("[DEBUG] minimize_scale_logspace: bounds_scale=[{:.6e}, {:.2}], bounds_a=[{:.6e}, {:.1}], binned.len()={}", 
              bounds_scale.0, bounds_scale.1, bounds_a.0, bounds_a.1, binned.len());
    
    eprintln!("[DEBUG] First 5 binned points: {:?}", &binned[..binned.len().min(5)]);
    
    let (scale, a, sse) = minimize_scale_logspace(
        &binned,
        bounds_scale.0,
        bounds_scale.1,
        bounds_a.0,
        bounds_a.1,
    );
    
    eprintln!("[DEBUG] minimize_scale_logspace result: scale={:.2}, a={:.6e}, sse={:.6e}", scale, a, sse);
    
    // Calculate R²
    let mut y_obs: Vec<f64> = Vec::new();
    let mut y_pred: Vec<f64> = Vec::new();
    for &(d, mr_obs) in &binned {
        y_obs.push(mr_obs);
        y_pred.push(chi_square_distribution_simple(d, a, scale));
    }
    let mean_y = y_obs.iter().sum::<f64>() / (y_obs.len() as f64);
    let ss_res = sse;
    let ss_tot = y_obs.iter().map(|yy| (yy - mean_y) * (yy - mean_y)).sum::<f64>();
    let r2 = if ss_tot > 0.0 { 1.0 - ss_res / ss_tot } else { f64::NEG_INFINITY };
    
    Some((a, scale, r2))
}

/// Calculate correction factor for a single position
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
    
    // Calculate expected mutation rate at this depth
    let y_at_depth = chi_square_distribution_simple(depth, a, scale);
    
    // Calculate correction factor
    if y_at_depth > target_y {
        // Depth too low, need to increase
        let distance = (y_at_depth - target_y) / target_y;
        let sqrt_distance = distance.sqrt();
        1.0 + sqrt_distance * weight_increase
    } else {
        // Depth too high, need to decrease
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

/// Apply PCR bias correction to a pileup CSV file
pub fn apply_pcr_bias_correction(
    args: &CorrectArgs,
    logger: &mut crate::Logger,
) -> Result<(), Box<dyn Error>> {
    logger.log("=== ModDetector PCR Bias Correction ===")?;
    logger.log(&format!("Input file: {}", args.input))?;
    logger.log(&format!("Output file: {}", args.output))?;
    logger.log(&format!("Weight increase: {}", args.weight_increase))?;
    logger.log(&format!("Weight decrease: {}", args.weight_decrease))?;
    
    // Read CSV file
    let file = File::open(&args.input)?;
    let reader = BufReader::new(file);
    let mut lines = reader.lines();
    
    // Read header
    let header = lines.next()
        .ok_or("Empty file")??;
    let header_fields: Vec<&str> = header.split(',').collect();
    
    // Find column indices
    let depth_idx = header_fields.iter().position(|&s| s == "depth")
        .ok_or("Missing 'depth' column")?;
    let mutation_count_idx = header_fields.iter().position(|&s| s == "rf_mutation_Count")
        .ok_or("Missing 'rf_mutation_Count' column")?;
    let species_idx = header_fields.iter().position(|&s| s == "species");
    
    // Collect data
    let mut depths = Vec::new();
    let mut mutation_rates = Vec::new();
    let mut records: Vec<Vec<String>> = Vec::new();
    
    for line in lines {
        let line = line?;
        let fields: Vec<String> = line.split(',').map(|s| s.to_string()).collect();
        
        // Filter by species if specified
        if let (Some(species_filter), Some(idx)) = (&args.species, species_idx) {
            if fields.get(idx) != Some(species_filter) {
                continue;
            }
        }
        
        // Parse depth and mutation count
        if let (Some(depth_str), Some(mut_count_str)) = (fields.get(depth_idx), fields.get(mutation_count_idx)) {
            if let (Ok(depth), Ok(mut_count)) = (depth_str.parse::<f64>(), mut_count_str.parse::<f64>()) {
                if depth > 0.0 {
                    let mutation_rate = mut_count / depth;
                    if mutation_rate > 0.0 && mutation_rate <= 0.1 {
                        depths.push(depth);
                        mutation_rates.push(mutation_rate);
                    }
                }
            }
        }
        
        records.push(fields);
    }
    
    logger.log(&format!("Loaded {} records", records.len()))?;
    logger.log(&format!("Valid data points: {}", depths.len()))?;
    
    if depths.len() < 10 {
        return Err("Not enough valid data points (need at least 10)".into());
    }
    
    // Fit Chi-Square distribution
    logger.log("Fitting Chi-Square distribution...")?;
    let chi2_params = fit_chi_square_to_data(&depths, &mutation_rates);
    
    if chi2_params.is_none() {
        return Err("Failed to fit Chi-Square distribution".into());
    }
    
    let (a, scale, r2) = chi2_params.unwrap();
    logger.log(&format!("Chi-Square fit: a={:.6e}, scale={:.2}, R²={:.4}", a, scale, r2))?;
    
    // Calculate target_y (minimum area line)
    let max_depth = depths.iter().fold(0.0_f64, |a, &b| a.max(b));
    let x_fit: Vec<f64> = (0..200)
        .map(|i| (i as f64 / 199.0) * max_depth * 1.2)
        .collect();
    let y_fit: Vec<f64> = x_fit.iter().map(|&x| chi_square_distribution_simple(x, a, scale)).collect();
    let mut target_y = calculate_min_area_line(&y_fit);
    
    if target_y <= 1e-6 {
        let mut rates_sorted = mutation_rates.clone();
        rates_sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        target_y = rates_sorted[rates_sorted.len() / 2];
        logger.log(&format!("Warning: Target line too small, using data median: y={:.6}", target_y))?;
    } else {
        logger.log(&format!("Target line (min area): y={:.6}", target_y))?;
    }
    
    // Apply correction to all records
    logger.log("Applying correction factors...")?;
    let mut output_file = BufWriter::new(File::create(&args.output)?);
    
    // Write header with effective_depth column
    let header_fields_vec: Vec<&str> = header.split(',').collect();
    let mut output_header = header.clone();
    let has_effective_depth = header_fields_vec.iter().any(|&s| s == "effective_depth");
    
    if !has_effective_depth {
        // Insert effective_depth after depth
        let mut header_fields: Vec<&str> = header_fields_vec.clone();
        if let Some(depth_pos) = header_fields.iter().position(|&s| s == "depth") {
            header_fields.insert(depth_pos + 1, "effective_depth");
        } else {
            header_fields.push("effective_depth");
        }
        output_header = header_fields.join(",");
    }
    writeln!(output_file, "{}", output_header)?;
    
    let mut correction_factors = Vec::new();
    
    for fields in records {
        let depth = if let Some(depth_str) = fields.get(depth_idx) {
            depth_str.parse::<f64>().unwrap_or(0.0)
        } else {
            0.0
        };
        
        let mutation_rate = if depth > 0.0 {
            if let Some(mut_count_str) = fields.get(mutation_count_idx) {
                mut_count_str.parse::<f64>().unwrap_or(0.0) / depth
            } else {
                0.0
            }
        } else {
            0.0
        };
        
        let effective_depth = if mutation_rate > 0.0 && mutation_rate <= 0.1 {
            let correction_factor = calculate_correction_factor(
                depth,
                mutation_rate,
                Some((a, scale)),
                target_y,
                args.weight_increase,
                args.weight_decrease,
            );
            
            let final_correction = correction_factor.max(0.5_f64).min(10.0_f64);
            let eff_depth = ((depth * final_correction).round() as usize).max(1);
            correction_factors.push(correction_factor);
            eff_depth
        } else {
            depth.max(1.0) as usize
        };
        
        // Write record with effective_depth
        let mut output_fields = fields.clone();
        if let Some(depth_pos) = header_fields_vec.iter().position(|&s| s == "depth") {
            // Check if effective_depth column already exists in header
            if let Some(eff_depth_pos) = header_fields_vec.iter().position(|&s| s == "effective_depth") {
                // effective_depth column exists, update it
                if output_fields.len() > eff_depth_pos {
                    output_fields[eff_depth_pos] = effective_depth.to_string();
                } else {
                    // Pad if needed
                    while output_fields.len() <= eff_depth_pos {
                        output_fields.push(String::new());
                    }
                    output_fields[eff_depth_pos] = effective_depth.to_string();
                }
            } else {
                // Insert effective_depth after depth
                if output_fields.len() > depth_pos {
                    output_fields.insert(depth_pos + 1, effective_depth.to_string());
                } else {
                    output_fields.push(effective_depth.to_string());
                }
            }
        } else {
            output_fields.push(effective_depth.to_string());
        }
        
        writeln!(output_file, "{}", output_fields.join(","))?;
    }
    
    // Statistics
    if !correction_factors.is_empty() {
        let mean_cf = correction_factors.iter().sum::<f64>() / correction_factors.len() as f64;
        let mut sorted_cf = correction_factors.clone();
        sorted_cf.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        let median_cf = sorted_cf[sorted_cf.len() / 2];
        let min_cf = correction_factors.iter().fold(f64::INFINITY, |a, &b| a.min(b));
        let max_cf = correction_factors.iter().fold(0.0_f64, |a, &b| a.max(b));
        
        logger.log(&format!("Correction statistics:"))?;
        logger.log(&format!("  Mean correction factor: {:.4}", mean_cf))?;
        logger.log(&format!("  Median correction factor: {:.4}", median_cf))?;
        logger.log(&format!("  Min correction factor: {:.4}", min_cf))?;
        logger.log(&format!("  Max correction factor: {:.4}", max_cf))?;
    }
    
    logger.log(&format!("Corrected pileup saved to: {}", args.output))?;
    
    Ok(())
}

