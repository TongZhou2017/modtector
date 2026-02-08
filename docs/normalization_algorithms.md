# Normalization Algorithms for RNA Structure Probing Data

## Abstract

This document provides a comprehensive description of the normalization algorithms implemented in the ModDetector software for processing RNA structure probing reactivity data. The normalization methods are designed to eliminate experimental batch effects and systematic biases, enabling quantitative comparison of reactivity scores across different samples and experimental conditions. Three primary normalization algorithms are implemented: Percentile28 (2-8% method), Winsor90 (90% Winsorizing), and Boxplot normalization. Additionally, an optional piecewise linear mapping transformation (Zarringhalam method) can be applied post-normalization to optimize reactivity score distributions.

## 1. Introduction

RNA structure probing experiments generate reactivity scores that reflect the accessibility and structural characteristics of RNA molecules. However, raw reactivity values are subject to experimental variations, batch effects, and systematic biases that must be corrected through normalization. The normalization algorithms described herein operate on reactivity value vectors extracted from genomic regions, where each region is defined by a unique chromosome identifier and strand orientation.

## 2. Input Data Specification

### 2.1 Data Format

The input data is provided in CSV (Comma-Separated Values) format with the following structure:

- **Column 1**: Chromosome/Reference sequence identifier (ChrID) - `string`
- **Column 2**: Strand orientation (Strand) - `string` (typically "+" or "-")
- **Column 3**: Position information - `numeric` or `string`
- **Column 4**: Stop signal reactivity score (stop_reactivity) - `numeric` (non-negative float)
- **Column 5**: Mutation signal reactivity score (mutation_reactivity) - `numeric` (non-negative float)
- **Column 6+**: Additional metadata (optional)

### 2.2 Input Parameters

For each normalization algorithm, the following parameters are required:

- **Input vector**: `V = {v₁, v₂, ..., vₙ}`, where `vᵢ ≥ 0` for all `i ∈ [1, n]` and `n` is the number of data points
- **Minimum sample size**: Algorithm-specific threshold (see Section 3 for details)

### 2.3 Data Preprocessing

Prior to normalization, the input reactivity values are:
1. Extracted from the appropriate column (stop_reactivity or mutation_reactivity)
2. Grouped by region identifier (ChrID|Strand)
3. Formed into numerical vectors for each region
4. Validated for non-negative values and missing data handling

## 3. Normalization Algorithms

### 3.1 Percentile28 Normalization (2-8% Method)

#### 3.1.1 Algorithm Description

The Percentile28 normalization method, also referred to as the 2-8% method, normalizes reactivity scores by dividing each value by the mean of reactivity values within the 2nd to 10th percentile range. This approach uses the lower percentile interval as a normalization baseline, effectively reducing the influence of extreme values while preserving the main distribution characteristics of the data.

#### 3.1.2 Mathematical Formulation

Given an input vector `V = {v₁, v₂, ..., vₙ}` of `n` reactivity values:

**Step 1: Sorting**
Sort the input vector in ascending order to obtain the sorted vector:
```
V_sorted = {v_(1), v_(2), ..., v_(n)}
```
where `v_(1) ≤ v_(2) ≤ ... ≤ v_(n)`.

**Step 2: Minimum Sample Size Check**
If `n < 10`, return the original vector without normalization:
```
V_norm = V
```

**Step 3: Percentile Index Calculation**
Calculate the position indices for the 2nd and 10th percentiles:
```
p₂ = ⌈n × 0.02⌉
p₁₀ = ⌈n × 0.10⌉
```
where `⌈·⌉` denotes the ceiling function. To ensure array bounds are not exceeded:
```
p₂ = min(p₂, n - 1)
p₁₀ = min(p₁₀, n - 1)
```

**Step 4: Normalization Factor Calculation**
Calculate the mean value within the percentile interval [p₂, p₁₀):
```
avg = { (1/(p₁₀ - p₂)) × Σᵢ₌ₚ₂^(p₁₀-1) v_(i)   if p₁₀ > p₂ and Σ > 0
      { 1.0                                        otherwise
```

**Step 5: Normalization**
Apply normalization to each value in the original (unsorted) vector:
```
v_norm,i = vᵢ / avg    for all i ∈ [1, n]
```

#### 3.1.3 Input-Output Specification

- **Input**: 
  - Vector `V` of length `n` with non-negative reactivity values
  - Minimum requirement: `n ≥ 1` (algorithm handles `n < 10` by returning original values)

- **Output**: 
  - Normalized vector `V_norm` of length `n` with normalized reactivity scores
  - Output range: `[0, +∞)` (typically `[0, 10]` for well-behaved data)

- **Parameters**:
  - No user-configurable parameters
  - Fixed percentile thresholds: 2% (lower) and 10% (upper)

#### 3.1.4 Algorithm Properties

- **Robustness**: Moderate - uses lower percentiles to reduce outlier influence
- **Sensitivity**: Low sensitivity to extreme high values
- **Computational complexity**: O(n log n) due to sorting step
- **Edge cases**: Returns original values if `n < 10` or if all values in percentile interval are zero

### 3.2 Winsor90 Normalization (90% Winsorizing)

#### 3.2.1 Algorithm Description

The Winsor90 normalization method combines Winsorization (tail trimming) with percentile-based normalization. The algorithm first truncates extreme values at the 5th and 95th percentiles, then normalizes the Winsorized values by dividing by the 95th percentile value. This approach effectively handles outliers while maintaining the distribution characteristics of the central data.

#### 3.2.2 Mathematical Formulation

Given an input vector `V = {v₁, v₂, ..., vₙ}` of `n` reactivity values:

**Step 1: Sorting**
Sort the input vector in ascending order:
```
V_sorted = {v_(1), v_(2), ..., v_(n)}
```

**Step 2: Minimum Sample Size Check**
If `n < 2`, return the original vector:
```
V_norm = V
```

**Step 3: Percentile Index Calculation**
Calculate the position indices for the 5th and 95th percentiles:
```
p₅ = ⌊n × 0.05⌋
p₉₅ = ⌈n × 0.95⌉
```
where `⌊·⌋` and `⌈·⌉` denote the floor and ceiling functions, respectively. Constrain indices:
```
p₅ = min(p₅, n - 1)
p₉₅ = min(p₉₅, n - 1)
```

**Step 4: Percentile Value Extraction**
Extract the percentile threshold values:
```
v₅ = V_sorted[p₅]
v₉₅ = V_sorted[p₉₅]
```

**Step 5: Winsorization**
Apply Winsorization to each value in the original vector:
```
v_winsor,i = { v₅      if vᵢ < v₅
             { v₉₅     if vᵢ > v₉₅
             { vᵢ       otherwise
```

**Step 6: Normalization Factor Calculation**
Calculate the normalization factor:
```
norm_val = max(v₉₅, 10⁻⁶)
```
The minimum value of `10⁻⁶` prevents division by zero.

**Step 7: Normalization**
Apply normalization:
```
v_norm,i = v_winsor,i / norm_val    for all i ∈ [1, n]
```

#### 3.2.3 Input-Output Specification

- **Input**: 
  - Vector `V` of length `n` with non-negative reactivity values
  - Minimum requirement: `n ≥ 1` (algorithm handles `n < 2` by returning original values)

- **Output**: 
  - Normalized vector `V_norm` of length `n` with normalized reactivity scores
  - Output range: `[0, 1]` (values above 95th percentile are truncated and normalized to 1.0)

- **Parameters**:
  - No user-configurable parameters
  - Fixed percentile thresholds: 5% (lower) and 95% (upper)
  - Minimum normalization factor: `ε = 10⁻⁶`

#### 3.2.4 Algorithm Properties

- **Robustness**: High - truncates extreme values before normalization
- **Sensitivity**: Low sensitivity to outliers due to Winsorization
- **Computational complexity**: O(n log n) due to sorting step
- **Edge cases**: Returns original values if `n < 2`; handles zero/negative values through minimum normalization factor

### 3.3 Boxplot Normalization

#### 3.3.1 Algorithm Description

The Boxplot normalization method employs boxplot statistical principles to identify and filter extreme values before normalization. The algorithm calculates the interquartile range (IQR) and identifies outliers using the standard boxplot criterion (values exceeding Q₃ + 1.5×IQR). After filtering outliers, the mean of the top 10% of remaining values is used as the normalization factor. This approach emphasizes high-reactivity signals while robustly handling outliers.

#### 3.3.2 Mathematical Formulation

Given an input vector `V = {v₁, v₂, ..., vₙ}` of `n` reactivity values:

**Step 1: Sorting**
Sort the input vector in ascending order:
```
V_sorted = {v_(1), v_(2), ..., v_(n)}
```

**Step 2: Minimum Sample Size Check**
If `n < 4`, return the original vector:
```
V_norm = V
```

**Step 3: Quartile Index Calculation**
Calculate the position indices for the first quartile (Q₁) and third quartile (Q₃):
```
q₁ = ⌈n × 0.25⌉
q₃ = ⌈n × 0.75⌉
```
Constrain indices:
```
q₁ = min(q₁, n - 1)
q₃ = min(q₃, n - 1)
```

**Step 4: Quartile Value Extraction**
Extract quartile values:
```
Q₁ = V_sorted[q₁]
Q₃ = V_sorted[q₃]
```

**Step 5: Interquartile Range and Upper Bound Calculation**
Calculate the interquartile range and upper bound for outlier detection:
```
IQR = Q₃ - Q₁
upper = Q₃ + 1.5 × IQR
```

**Step 6: Outlier Filtering**
Filter out values exceeding the upper bound:
```
V_filtered = {v_(i) ∈ V_sorted : v_(i) ≤ upper}
```
Let `m = |V_filtered|` denote the size of the filtered dataset.

**Step 7: Top 10% Mean Calculation**
Calculate the number of top 10% values:
```
top10 = ⌈m × 0.10⌉
```

Calculate the mean of the top 10% values in the filtered dataset:
```
avg = { (1/top10) × Σᵢ₌ₘ₋ₜₒₚ₁₀^(m-1) V_filtered[i]   if top10 > 0 and Σ > 0
      { 1.0                                              otherwise
```

**Step 8: Normalization**
Apply normalization using the original (unfiltered) values:
```
v_norm,i = vᵢ / avg    for all i ∈ [1, n]
```

Note: The normalization factor is calculated from filtered data, but normalization is applied to all original values.

#### 3.3.3 Input-Output Specification

- **Input**: 
  - Vector `V` of length `n` with non-negative reactivity values
  - Minimum requirement: `n ≥ 1` (algorithm handles `n < 4` by returning original values)

- **Output**: 
  - Normalized vector `V_norm` of length `n` with normalized reactivity scores
  - Output range: `[0, +∞)` (typically `[0, 5]` for well-behaved data after outlier filtering)

- **Parameters**:
  - No user-configurable parameters
  - Fixed quartile thresholds: 25% (Q₁) and 75% (Q₃)
  - Outlier detection multiplier: `k = 1.5` (standard boxplot criterion)
  - Top percentile for normalization factor: 10%

#### 3.3.4 Algorithm Properties

- **Robustness**: Very high - filters outliers using statistical principles
- **Sensitivity**: Low sensitivity to extreme outliers due to filtering
- **Computational complexity**: O(n log n) due to sorting and filtering steps
- **Edge cases**: Returns original values if `n < 4`; handles empty filtered dataset by using default normalization factor

## 4. Post-Normalization Transformation

### 4.1 Zarringhalam Piecewise Linear Mapping

#### 4.1.1 Algorithm Description

The Zarringhalam piecewise linear mapping is an optional post-normalization transformation that remaps normalized reactivity scores from the range [0, 1] to a new distribution through four piecewise linear functions. This transformation adjusts the score distribution to be more uniform across low, medium, and high reactivity regions, facilitating subsequent structure prediction and analysis.

#### 4.1.2 Mathematical Formulation

Given a normalized reactivity score `x ∈ [0, +∞)`, the piecewise linear mapping function `f(x)` is defined as:

```
f(x) = { (x / 0.25) × 0.35                    if x < 0.25
       { 0.35 + ((x - 0.25) / 0.05) × 0.20    if 0.25 ≤ x < 0.3
       { 0.55 + ((x - 0.3) / 0.4) × 0.30      if 0.3 ≤ x < 0.7
       { 0.85 + ((x - 0.7) / 0.3) × 0.15      if x ≥ 0.7
```

This can be expressed more compactly as:

```
f(x) = { (x / 0.25) × 0.35                    if x ∈ [0, 0.25)
       { 0.35 + (x - 0.25) × 4.0              if x ∈ [0.25, 0.3)
       { 0.55 + (x - 0.3) × 0.75              if x ∈ [0.3, 0.7)
       { 0.85 + (x - 0.7) × 0.5               if x ∈ [0.7, +∞)
```

#### 4.1.3 Mapping Intervals

The piecewise linear mapping divides the input range into four intervals with different slopes:

1. **Interval 1** (`x ∈ [0, 0.25)`): Maps to `[0, 0.35)` with slope `1.4`
2. **Interval 2** (`x ∈ [0.25, 0.3)`): Maps to `[0.35, 0.55)` with slope `4.0`
3. **Interval 3** (`x ∈ [0.3, 0.7)`): Maps to `[0.55, 0.85)` with slope `0.75`
4. **Interval 4** (`x ∈ [0.7, +∞)`): Maps to `[0.85, +∞)` with slope `0.5`

#### 4.1.4 Input-Output Specification

- **Input**: 
  - Normalized reactivity score vector `V_norm` (output from any normalization algorithm)
  - Values typically in range `[0, +∞)`, though designed for `[0, 1]` range

- **Output**: 
  - Transformed reactivity score vector `V_transformed` of same length
  - Output range: `[0, +∞)` with enhanced separation in different reactivity regions

- **Parameters**:
  - Fixed breakpoints: `{0.25, 0.3, 0.7}`
  - Fixed target values: `{0.35, 0.55, 0.85, 1.0}`

#### 4.1.5 Algorithm Properties

- **Purpose**: Redistribute reactivity scores to improve separation between reactivity levels
- **Linearity**: Piecewise linear transformation preserves relative ordering
- **Computational complexity**: O(n) - single pass through data
- **Continuity**: Function is continuous at breakpoints

## 5. Implementation Details

### 5.1 Regional Grouping Strategy

Normalization is performed independently for each genomic region, where a region is defined by the combination of chromosome identifier (ChrID) and strand orientation. This ensures that:

1. Data from different chromosomes are normalized independently
2. Forward and reverse strand data are processed separately
3. Regional variations in signal intensity do not affect normalization of other regions

The region identifier format is: `"ChrID|Strand"`

### 5.2 Independent Processing of Signal Types

Stop signals and mutation signals are normalized independently:

- **Stop signals**: Extracted from column 4 (stop_reactivity), normalized, and written back to column 4
- **Mutation signals**: Extracted from column 5 (mutation_reactivity), normalized, and written back to column 5

This independent processing allows each signal type to be optimized according to its own distribution characteristics.

### 5.3 Window-Based Normalization (Optional)

The implementation supports optional window-based normalization modes:

- **Fixed window mode**: Normalization within fixed-size sliding windows
- **Dynamic window mode**: Adaptive window sizing based on reactive base density

These modes are available for advanced use cases but are not the default behavior.

## 6. Output Data Specification

### 6.1 Output Format

The output maintains the same CSV format and structure as the input, with the following modifications:

- **Column 4** (stop_reactivity): Contains normalized stop signal reactivity scores
- **Column 5** (mutation_reactivity): Contains normalized mutation signal reactivity scores
- **Other columns**: Unchanged from input

### 6.2 Output Precision

Normalized reactivity scores are output with 6 decimal places precision, formatted as:
```
"0.000000"
```

### 6.3 Data Validation

The implementation handles edge cases:

- **Missing values**: Replaced with `"0.000000"`
- **Invalid values** (NaN, empty strings): Replaced with `"0.000000"`
- **Insufficient data**: Returns original values (algorithm-specific minimum sample size requirements)

## 7. Parameter Configuration

### 7.1 Normalization Method Selection

The normalization method is selected via the `--method` parameter with the following options:

- `percentile28`, `2-8`, or `2_8`: Percentile28 normalization
- `winsor90` or `winsorizing`: Winsor90 normalization
- `boxplot`: Boxplot normalization

### 7.2 Additional Parameters

- **`--cutoff`**: SNP threshold for mutation signal filtering (range: [0.0, 1.0], default: 0.25)
- **`--bases`**: Reactive base types to consider (default: "ACGT", options: any combination of A, C, G, T)
- **`--linear`**: Enable Zarringhalam piecewise linear mapping (boolean, default: false)
- **`--window`**: Fixed window size for window-based normalization (default: 0, disabled)
- **`--window-offset`**: Window offset for fixed window mode (default: 0)
- **`--dynamic`**: Enable dynamic window mode (boolean, default: false)
- **`--window-size-dynamic`**: Minimum reactive bases per dynamic window (range: [1, 1000], default: 50)
- **`--log`**: Log file path (optional, default: "norm.log")

## 8. Algorithm Comparison

| Property | Percentile28 | Winsor90 | Boxplot |
|----------|--------------|----------|---------|
| **Robustness to outliers** | Moderate | High | Very High |
| **Normalization baseline** | Lower percentiles (2-10%) | Upper percentile (95%) | Top 10% of filtered data |
| **Outlier handling** | None | Truncation | Filtering |
| **Output range** | [0, +∞) | [0, 1] | [0, +∞) |
| **Minimum sample size** | 10 | 2 | 4 |
| **Computational complexity** | O(n log n) | O(n log n) | O(n log n) |
| **Best for** | General purpose, moderate outliers | High outlier presence | Extreme outliers, high-reactivity emphasis |

## 9. References and Notes

### 9.1 Algorithm Origins

- **Percentile28**: Based on the 2-8% normalization method commonly used in SHAPE-MaP and related RNA structure probing analyses
- **Winsor90**: Implements standard Winsorization with 90% central range (5th to 95th percentiles)
- **Boxplot**: Based on standard boxplot outlier detection (Tukey's method) with IQR-based filtering
- **Zarringhalam mapping**: Piecewise linear transformation method for reactivity score redistribution

### 9.2 Implementation Notes

- All algorithms handle edge cases (small sample sizes, zero values, etc.) gracefully
- Sorting is performed using stable partial comparison to handle floating-point precision issues
- Array bounds are carefully checked to prevent index out-of-bounds errors
- Default normalization factors (typically 1.0) are used when calculation is not possible

## 10. Conclusion

The normalization algorithms implemented in ModDetector provide robust methods for processing RNA structure probing reactivity data. Each algorithm offers different trade-offs in terms of robustness, sensitivity, and computational requirements, allowing users to select the most appropriate method for their specific data characteristics and analysis goals. The optional piecewise linear mapping further enhances the utility of normalized scores for downstream structure prediction and analysis applications.






