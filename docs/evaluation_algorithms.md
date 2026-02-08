# Accuracy Evaluation Algorithms for RNA Structure Probing Data

## Abstract

This document provides a comprehensive description of the accuracy evaluation algorithms implemented in the ModDetector software for assessing the performance of RNA structure probing reactivity data in predicting RNA secondary structure. The evaluation framework employs binary classification metrics to compare predicted reactivity scores against experimentally determined or computationally predicted secondary structure annotations. The algorithms feature sophisticated sequence alignment mechanisms, including automatic shift correction for handling positional discrepancies between reactivity data and reference secondary structures. Multiple evaluation modes are supported, including overall assessment, base-specific analysis, and signal-type-specific evaluation (stop vs. mutation signals).

## 1. Introduction

RNA structure probing experiments generate reactivity scores that reflect the accessibility of individual nucleotides, with higher reactivity typically associated with single-stranded regions and lower reactivity with double-stranded (base-paired) regions. The accuracy evaluation algorithms described herein assess the predictive power of reactivity scores by treating secondary structure prediction as a binary classification problem: single-stranded (positive class) versus double-stranded (negative class). The evaluation framework implements comprehensive metrics including receiver operating characteristic (ROC) analysis, precision-recall curves, and various classification performance measures.

A key challenge in evaluating reactivity data against reference secondary structures is the potential misalignment between reactivity measurement positions and reference structure positions. This misalignment can arise from sequence differences, annotation offsets, or experimental artifacts. The algorithms address this challenge through automatic shift correction algorithms that optimize sequence alignment based on base matching.

## 2. Data Structures and Input Specification

### 2.1 Secondary Structure Representation

The reference secondary structure is represented using the `SecondaryStructure` data structure:

```rust
pub struct SecondaryStructure {
    pub sequence: String,              // RNA sequence (ACGU characters)
    pub structure: String,              // Dot-bracket notation (e.g., "..((..))..")
    pub positions: Vec<u32>,            // 1-based positions of valid bases
    pub is_single_stranded: Vec<bool>,  // Boolean vector: true = single-stranded
}
```

**Structure Notation**:
- `.` (dot): Single-stranded (unpaired) nucleotide
- `(`, `)`, `[`, `]`, `{`, `}`, `<`, `>`: Paired nucleotides (double-stranded)
- Other characters (e.g., spaces, `a`, `A`): Ignored during parsing

**Position Mapping**:
- Only positions corresponding to standard RNA bases (A, C, G, U) are included
- Positions are 1-based (first base = position 1)
- The `is_single_stranded` vector has the same length as `positions`, with `true` indicating single-stranded and `false` indicating double-stranded

### 2.2 Reactivity Data Format

Reactivity data is provided in CSV format with the following structure:

- **Column 1**: Chromosome/Reference sequence identifier (ChrID) - `string`
- **Column 2**: Strand orientation (Strand) - `string` (typically "+" or "-")
- **Column 3**: Position - `u32` (1-based genomic position)
- **Column 4**: Base - `char` (A, C, G, T, or U)
- **Column 5**: Stop signal reactivity (reactivity_stop) - `f64` (non-negative float)
- **Column 6**: Mutation signal reactivity (reactivity_mutation) - `f64` (non-negative float)

### 2.3 Input Parameters

For evaluation algorithms, the following parameters are required:

- **Reactivity file path**: Path to CSV file containing reactivity data
- **Secondary structure file path**: Path to file containing reference secondary structure
- **Gene ID**: Identifier for the target gene/transcript
- **Strand**: Strand orientation ("+" or "-")
- **Signal type**: Either "stop" or "mutation" to specify which reactivity column to use
- **Output directory**: Directory for saving evaluation results and plots

### 2.4 Evaluation Result Structure

The evaluation results are stored in the `EvaluationResult` structure:

```rust
pub struct EvaluationResult {
    pub auc: f64,                    // Area Under ROC Curve
    pub f1_score: f64,              // F1-score (harmonic mean of precision and recall)
    pub sensitivity: f64,            // True Positive Rate (TPR)
    pub specificity: f64,            // True Negative Rate (TNR)
    pub ppv: f64,                    // Positive Predictive Value (Precision)
    pub npv: f64,                    // Negative Predictive Value
    pub accuracy: f64,               // Overall classification accuracy
    pub thresholds: Vec<f64>,         // Threshold values used for ROC/PR curves
    pub tpr: Vec<f64>,               // True Positive Rate at each threshold
    pub fpr: Vec<f64>,               // False Positive Rate at each threshold
    pub precision: Vec<f64>,         // Precision at each threshold
    pub recall: Vec<f64>,            // Recall at each threshold
}
```

## 3. Secondary Structure Alignment Algorithms

### 3.1 Automatic Shift Correction Algorithm

#### 3.1.1 Algorithm Description

The automatic shift correction algorithm addresses positional misalignment between reactivity data and reference secondary structures. This misalignment can occur due to sequence differences, annotation offsets, or experimental artifacts. The algorithm performs an exhaustive search over a range of shift values to find the optimal alignment that maximizes base matching between the reactivity sequence and the reference structure sequence.

#### 3.1.2 Mathematical Formulation

Given:
- Reactivity data: `R = {(pos₁, base₁, react₁), (pos₂, base₂, react₂), ..., (posₘ, baseₘ, reactₘ)}`
- Secondary structure: `S = {(pos₁', base₁', struct₁'), (pos₂', base₂', struct₂'), ..., (posₙ', baseₙ', structₙ')}`

where `structᵢ' ∈ {single, double}` indicates the structural state.

**Step 1: Structure Base Extraction**

Extract valid bases from the secondary structure sequence:
```
S_bases = {(pos, base, is_single) | base ∈ {A, C, G, U}}
```

where `is_single = true` if the corresponding structure character is `.`, and `is_single = false` otherwise.

**Step 2: Shift Range Calculation**

Calculate the maximum shift range based on length difference:
```
length_diff = |m - n|
max_shift = min(|length_diff| + 5, 20)
```

The shift range is: `shift ∈ [-max_shift, max_shift]`

**Step 3: Base Matching with T/U Equivalence**

For each reactivity position `posᵢ` and base `baseᵢ`, and for each shift value `s ∈ [-max_shift, max_shift]`:

Calculate the shifted position:
```
shifted_pos = posᵢ + s
```

If `shifted_pos ∈ [1, n]`, check base matching:
```
base_match = (baseᵢ == baseⱼ') ∨ ((baseᵢ == 'T' ∧ baseⱼ' == 'U') ∨ (baseᵢ == 'U' ∧ baseⱼ' == 'T'))
```

where `j` is the index in `S_bases` corresponding to `shifted_pos`.

**Step 4: Optimal Shift Selection**

For each shift value `s`, calculate:
```
match_count(s) = Σᵢ I(shifted_posᵢ ∈ [1, n] ∧ base_match(shifted_posᵢ))
total_count(s) = Σᵢ I(shifted_posᵢ ∈ [1, n])
match_rate(s) = match_count(s) / total_count(s)
```

where `I(condition)` is the indicator function (1 if condition is true, 0 otherwise).

Select the optimal shift:
```
best_shift = argmax_{s ∈ [-max_shift, max_shift]} match_count(s)
```

**Step 5: Data Matching**

Using the optimal shift `best_shift`, create matched data pairs:
```
matched_data = {(reactᵢ, is_singleⱼ) | shifted_posᵢ ∈ [1, n] ∧ base_match(shifted_posᵢ)}
```

where `j` corresponds to the structure position matching `shifted_posᵢ`.

#### 3.3.3 Input-Output Specification

- **Input**:
  - Reactivity data vector: `R = {(pos, base, reactivity)}` of length `m`
  - Secondary structure object: `SecondaryStructure` with `n` valid positions
  - Signal type: `"stop"` or `"mutation"` (used for logging purposes)

- **Output**:
  - Matched data vector: `matched_data = {(reactivity, is_single)}` of length `k ≤ m`
  - Only positions with valid matches are included
  - Positions outside the structure range are excluded

- **Parameters**:
  - `max_shift`: Automatically calculated as `min(|m - n| + 5, 20)`
  - Base equivalence: T ↔ U treated as equivalent

#### 3.3.4 Algorithm Characteristics

1. **Base Equivalence Handling**: The algorithm treats T and U as equivalent bases, accommodating differences between DNA and RNA sequence representations.

2. **Robustness to Length Differences**: The shift range adapts to the length difference between reactivity and structure data, with a maximum limit of 20 positions.

3. **Match Rate Optimization**: The algorithm selects the shift that maximizes the number of matching bases, ensuring optimal alignment.

4. **Position Filtering**: Only positions within the valid structure range are included in the matched data, preventing out-of-bounds errors.

## 4. Accuracy Evaluation Metrics

### 4.1 Binary Classification Framework

The evaluation treats secondary structure prediction as a binary classification problem:

- **Positive class (P)**: Single-stranded nucleotides (`is_single = true`)
- **Negative class (N)**: Double-stranded nucleotides (`is_single = false`)
- **Prediction**: Based on reactivity threshold `τ`
  - If `reactivity ≥ τ`: Predicted as single-stranded (positive)
  - If `reactivity < τ`: Predicted as double-stranded (negative)

### 4.2 Confusion Matrix

For a given threshold `τ`, the classification results form a confusion matrix:

|                    | Predicted Positive | Predicted Negative |
|--------------------|--------------------|--------------------|
| **Actual Positive**| True Positive (TP) | False Negative (FN)|
| **Actual Negative**| False Positive (FP)| True Negative (TN) |

### 4.3 Classification Metrics

#### 4.3.1 Sensitivity (True Positive Rate, Recall)

Sensitivity measures the proportion of actual single-stranded nucleotides correctly identified:

```
Sensitivity = TPR = Recall = TP / (TP + FN)
```

**Range**: [0, 1]
- 1.0: All single-stranded nucleotides correctly identified
- 0.0: No single-stranded nucleotides identified

#### 4.3.2 Specificity (True Negative Rate)

Specificity measures the proportion of actual double-stranded nucleotides correctly identified:

```
Specificity = TNR = TN / (TN + FP)
```

**Range**: [0, 1]
- 1.0: All double-stranded nucleotides correctly identified
- 0.0: No double-stranded nucleotides identified

#### 4.3.3 Positive Predictive Value (Precision)

PPV measures the proportion of predicted single-stranded nucleotides that are actually single-stranded:

```
PPV = Precision = TP / (TP + FP)
```

**Range**: [0, 1]
- 1.0: All predicted single-stranded nucleotides are correct
- 0.0: No predicted single-stranded nucleotides are correct

#### 4.3.4 Negative Predictive Value

NPV measures the proportion of predicted double-stranded nucleotides that are actually double-stranded:

```
NPV = TN / (TN + FN)
```

**Range**: [0, 1]

#### 4.3.5 Accuracy

Overall classification accuracy:

```
Accuracy = (TP + TN) / (TP + TN + FP + FN)
```

**Range**: [0, 1]

#### 4.3.6 F1-Score

Harmonic mean of precision and recall:

```
F1 = 2 × (Precision × Recall) / (Precision + Recall)
     = 2 × (PPV × Sensitivity) / (PPV + Sensitivity)
```

**Range**: [0, 1]
- Provides a balanced measure when classes are imbalanced

### 4.4 ROC Curve and AUC

#### 4.4.1 ROC Curve Generation

The Receiver Operating Characteristic (ROC) curve plots True Positive Rate (TPR) against False Positive Rate (FPR) across all possible thresholds.

**Step 1: Threshold Selection**

For matched data `matched_data = {(reactᵢ, is_singleᵢ)}` of length `n`:

Extract unique reactivity values:
```
unique_values = {reactᵢ | i ∈ [1, n]} (sorted, deduplicated)
```

**Step 2: Metric Calculation for Each Threshold**

For each threshold `τ ∈ unique_values`:

```
TP(τ) = Σᵢ I(reactᵢ ≥ τ ∧ is_singleᵢ = true)
FP(τ) = Σᵢ I(reactᵢ ≥ τ ∧ is_singleᵢ = false)
TN(τ) = Σᵢ I(reactᵢ < τ ∧ is_singleᵢ = false)
FN(τ) = Σᵢ I(reactᵢ < τ ∧ is_singleᵢ = true)
```

Then calculate:
```
TPR(τ) = TP(τ) / (TP(τ) + FN(τ))
FPR(τ) = FP(τ) / (FP(τ) + TN(τ))
```

**Step 3: ROC Point Generation**

Generate ROC points:
```
ROC_points = {(0.0, 0.0)} ∪ {(FPR(τ), TPR(τ)) | τ ∈ unique_values} ∪ {(1.0, 1.0)}
```

**Step 4: Monotonicity Enforcement**

Sort ROC points by FPR and ensure monotonicity:
- If multiple points have the same FPR, keep the one with maximum TPR
- This ensures the ROC curve is monotonically non-decreasing

#### 4.4.2 Area Under ROC Curve (AUC)

AUC is calculated using the trapezoidal rule:

```
AUC = Σᵢ₌₁ⁿ⁻¹ (FPRᵢ₊₁ - FPRᵢ) × (TPRᵢ + TPRᵢ₊₁) / 2
```

where `n` is the number of ROC points, and points are sorted by FPR.

**Range**: [0, 1]
- 1.0: Perfect classifier
- 0.5: Random classifier (diagonal line)
- < 0.5: Worse than random (rare, indicates systematic misclassification)

### 4.5 Precision-Recall Curve

The Precision-Recall (PR) curve plots Precision (PPV) against Recall (Sensitivity) across all thresholds.

**PR Point Generation**:

For each threshold `τ ∈ unique_values`:
```
Precision(τ) = TP(τ) / (TP(τ) + FP(τ))
Recall(τ) = TPR(τ) = TP(τ) / (TP(τ) + FN(τ))
```

**PR Curve Points**:
```
PR_points = {(Recall(τ), Precision(τ)) | τ ∈ unique_values}
```

### 4.6 Optimal Threshold Selection

The optimal threshold is selected based on F1-score maximization:

```
best_threshold = argmax_{τ} F1(τ)
```

where:
```
F1(τ) = 2 × Precision(τ) × Recall(τ) / (Precision(τ) + Recall(τ))
```

The final classification metrics (sensitivity, specificity, PPV, NPV, accuracy, F1-score) are calculated using this optimal threshold.

## 5. Evaluation Algorithms

### 5.1 Overall Evaluation Algorithm

#### 5.1.1 Algorithm Description

The overall evaluation algorithm assesses the predictive performance of reactivity scores across all nucleotide positions, regardless of base type. This provides a global measure of how well reactivity data predicts secondary structure.

#### 5.1.2 Algorithm Workflow

**Step 1: Data Extraction and Alignment**

1. Extract reactivity data for the target gene and strand:
   ```
   reactivity_data = extract_reactivity_with_bases(file, gene_id, strand, signal_type)
   ```

2. Parse secondary structure file:
   ```
   secondary_structure = parse_secondary_structure(structure_file)
   ```

3. Perform automatic shift correction and alignment:
   ```
   matched_data = match_reactivity_with_structure_auto_shift(
       reactivity_data, secondary_structure, signal_type
   )
   ```

**Step 2: Data Sorting**

Sort matched data by reactivity value (ascending):
```
sorted_data = sort(matched_data, key=reactivity)
```

**Step 3: ROC and PR Curve Generation**

Generate ROC and PR curves as described in Sections 4.4 and 4.5.

**Step 4: Metric Calculation**

Calculate all classification metrics using the optimal threshold (Section 4.6).

#### 5.1.3 Input-Output Specification

- **Input**:
  - Reactivity CSV file path
  - Secondary structure file path
  - Gene ID, strand, signal type
  - Output directory path

- **Output**:
  - `EvaluationResult` structure with all metrics
  - ROC curve plot (SVG format)
  - Precision-Recall curve plot (SVG format)
  - Text file with detailed results

- **Parameters**:
  - Signal type: `"stop"` or `"mutation"`
  - Automatic shift correction: Enabled by default

### 5.2 Base-Specific Evaluation Algorithm

#### 5.2.1 Algorithm Description

The base-specific evaluation algorithm assesses predictive performance separately for each nucleotide type (A, C, G, T/U). This analysis reveals whether reactivity scores have different predictive power for different bases, which can inform understanding of base-specific reactivity patterns and potential biases.

#### 5.2.2 Algorithm Workflow

**Step 1: Data Extraction and Alignment**

Same as overall evaluation (Section 5.1.2, Step 1).

**Step 2: Base Grouping**

Group matched data by base type:
```
base_groups = {
    'A': [(reactᵢ, is_singleᵢ) | baseᵢ = 'A'],
    'C': [(reactᵢ, is_singleᵢ) | baseᵢ = 'C'],
    'G': [(reactᵢ, is_singleᵢ) | baseᵢ = 'G'],
    'T': [(reactᵢ, is_singleᵢ) | baseᵢ ∈ {'T', 'U'}],
    'U': [(reactᵢ, is_singleᵢ) | baseᵢ ∈ {'T', 'U'}]
}
```

**Step 3: Per-Base Evaluation**

For each base type `b ∈ {A, C, G, T, U}`:

1. Check minimum sample size: `|base_groups[b]| ≥ 10`
2. If sufficient, calculate evaluation metrics:
   ```
   base_result[b] = evaluate_reactivity_accuracy_from_matched_data(
       base_groups[b], signal_type
   )
   ```
3. Record base-specific statistics:
   ```
   single_stranded_count[b] = Σᵢ I(is_singleᵢ = true)
   double_stranded_count[b] = Σᵢ I(is_singleᵢ = false)
   ```

**Step 4: Overall Evaluation**

Calculate overall metrics using all matched data (same as Section 5.1.2).

**Step 5: Comprehensive Result Compilation**

Create `ComprehensiveEvaluationResult`:
```rust
pub struct ComprehensiveEvaluationResult {
    pub overall_result: EvaluationResult,
    pub base_specific_results: Vec<BaseSpecificEvaluationResult>,
    pub total_positions: usize,
}
```

where `BaseSpecificEvaluationResult` includes:
- Base type
- Evaluation metrics
- Data count
- Single-stranded and double-stranded counts

#### 5.2.3 Input-Output Specification

- **Input**: Same as overall evaluation (Section 5.1.3)

- **Output**:
  - Overall evaluation results (same as Section 5.1.3)
  - Per-base evaluation results for each base type with sufficient data (≥10 positions)
  - Combined ROC curve plot showing overall + all base-specific curves
  - Individual ROC and PR plots for each base type
  - Comprehensive text report with all metrics

- **Parameters**:
  - Minimum sample size per base: 10 positions (default)
  - Base types evaluated: A, C, G, T, U

### 5.3 Signal-Type-Specific Evaluation

#### 5.3.1 Algorithm Description

The evaluation framework supports two types of reactivity signals:
- **Stop signal**: Reactivity derived from reverse transcription stop events
- **Mutation signal**: Reactivity derived from mutation events

The algorithm can evaluate both signals separately or simultaneously, enabling comparison of their predictive performance.

#### 5.3.2 Algorithm Workflow

**Step 1: Dual Signal Extraction**

For each signal type `s ∈ {"stop", "mutation"}`:

1. Extract reactivity data using the appropriate column:
   ```
   reactivity_data[s] = extract_reactivity_with_bases(
       file, gene_id, strand, s
   )
   ```

2. Perform alignment:
   ```
   matched_data[s] = match_reactivity_with_structure_auto_shift(
       reactivity_data[s], secondary_structure, s
   )
   ```

**Step 2: Independent Evaluation**

For each signal type:
```
result[s] = evaluate_reactivity_accuracy_from_matched_data(
    matched_data[s], s
)
```

**Step 3: Performance Comparison**

Compare metrics:
```
if result["stop"].auc > result["mutation"].auc:
    better_signal = "stop"
elif result["mutation"].auc > result["stop"].auc:
    better_signal = "mutation"
else:
    better_signal = "equal"
```

**Step 4: Base-Specific Analysis**

Optionally perform base-specific evaluation for each signal type (Section 5.2).

#### 5.3.3 Input-Output Specification

- **Input**: Same as overall evaluation, but signal type can be specified as "both"

- **Output**:
  - Separate evaluation results for stop and mutation signals
  - Performance comparison report
  - Individual plots and reports for each signal type
  - Combined report summarizing both signals

## 6. Implementation Details

### 6.1 Secondary Structure Parsing

The secondary structure parser (`parse_secondary_structure`) handles various file formats:

1. **Sequence Lines**: Lines containing only ACGU characters (case-insensitive)
2. **Structure Lines**: Lines containing dot-bracket notation (`.`, `(`, `)`, `[`, `]`, `{`, `}`, `<`, `>`)
3. **Comment Lines**: Lines starting with `#` are ignored
4. **Empty Lines**: Ignored

The parser:
- Concatenates multiple sequence lines
- Concatenates multiple structure lines
- Removes spaces from structure strings
- Validates that sequence and structure lengths match
- Filters to include only positions with valid RNA bases (A, C, G, U)

### 6.2 Reactivity Data Extraction

The reactivity data extractor (`extract_reactivity_with_bases`) processes CSV files:

1. Reads header row (skipped)
2. For each data row:
   - Validates chromosome ID and strand match
   - Extracts position, base, and reactivity value
   - Selects reactivity column based on signal type:
     - `signal_type == "mutation"`: Column 6 (reactivity_mutation)
     - `signal_type == "stop"`: Column 5 (reactivity_stop)
3. Sorts data by position (ascending)

### 6.3 Threshold Selection Strategy

The algorithm uses all unique reactivity values as thresholds, ensuring:
- Each threshold corresponds to an actual data point
- No interpolation is needed
- ROC and PR curves are based on empirical data

For optimal threshold selection, the algorithm:
1. Calculates F1-score for each threshold
2. Selects the threshold with maximum F1-score
3. Uses this threshold for final metric calculation

### 6.4 Numerical Stability

The algorithms include safeguards for edge cases:

1. **Division by Zero**: All division operations check for zero denominators:
   ```
   metric = if denominator > 0 then numerator / denominator else 0.0
   ```

2. **Empty Data**: Functions return errors if no matched data is found

3. **Insufficient Samples**: Base-specific evaluation requires minimum 10 positions per base

4. **Monotonic ROC Curves**: ROC points are sorted and deduplicated to ensure monotonicity

## 7. Output Formats

### 7.1 Text Results File

The evaluation results are saved in a tab-separated text file with the following structure:

```
# Reactivity Accuracy Evaluation Results
# Signal Type: {signal_type}
#
# Overall Metrics:
AUC	{value}
F1-score	{value}
Sensitivity	{value}
Specificity	{value}
PPV	{value}
NPV	{value}
Accuracy	{value}
#
# ROC Curve Data:
Threshold	TPR	FPR	Precision	Recall
{threshold₁}	{tpr₁}	{fpr₁}	{precision₁}	{recall₁}
{threshold₂}	{tpr₂}	{fpr₂}	{precision₂}	{recall₂}
...
```

### 7.2 Comprehensive Results File

For base-specific evaluations, a comprehensive results file includes:

```
=== Base-Specific Evaluation Results ===
Gene ID: {gene_id}
Strand: {strand}
Signal Type: {signal_type}
Total Positions: {total_positions}

=== Overall Evaluation Results ===
AUC: {value}
F1-score: {value}
Sensitivity: {value}
Specificity: {value}
PPV: {value}
NPV: {value}
Accuracy: {value}

=== Per-Base Evaluation Results ===
--- Base A ---
Position Count: {count} (Single-stranded: {single_count}, Double-stranded: {double_count})
AUC: {value}
F1-score: {value}
...
```

### 7.3 Visualization Outputs

1. **ROC Curve Plot** (SVG format):
   - X-axis: False Positive Rate (0 to 1)
   - Y-axis: True Positive Rate (0 to 1)
   - Diagonal reference line (random classifier)
   - AUC value annotation

2. **Precision-Recall Curve Plot** (SVG format):
   - X-axis: Recall (0 to 1)
   - Y-axis: Precision (0 to 1)
   - F1-score annotation

3. **Combined ROC Curve Plot** (for base-specific evaluation):
   - Multiple curves: Overall + A, C, G, T/U bases
   - Color-coded by base type
   - Legend with AUC values

## 8. Algorithm Complexity

### 8.1 Time Complexity

- **Shift Correction**: O(m × n × s), where:
  - `m`: Number of reactivity positions
  - `n`: Number of structure positions
  - `s`: Shift range (typically ≤ 41: -20 to +20)
  
- **ROC Curve Generation**: O(n × log(n)), where:
  - `n`: Number of matched positions
  - Dominated by sorting unique reactivity values

- **Metric Calculation**: O(n × t), where:
  - `n`: Number of matched positions
  - `t`: Number of unique thresholds (typically ≤ n)

- **Overall Complexity**: O(m × n × s + n × log(n) + n × t) ≈ O(m × n × s) for typical cases

### 8.2 Space Complexity

- **Data Storage**: O(m + n) for reactivity and structure data
- **Matched Data**: O(k), where k ≤ min(m, n)
- **ROC/PR Curves**: O(t), where t is the number of unique thresholds
- **Overall**: O(m + n + t)

## 9. Limitations and Considerations

### 9.1 Alignment Limitations

1. **Linear Shift Only**: The algorithm assumes a constant offset (shift) between reactivity and structure positions. Non-linear misalignments (e.g., insertions/deletions) are not handled.

2. **Shift Range**: Maximum shift is limited to 20 positions. Larger misalignments may not be corrected.

3. **Base Matching Requirement**: Positions are only matched if bases match (with T/U equivalence). Mismatched bases are excluded from evaluation.

### 9.2 Evaluation Assumptions

1. **Binary Classification**: The algorithm treats structure as binary (single- vs. double-stranded). More complex structural states (e.g., tertiary interactions) are not considered.

2. **Threshold Independence**: Metrics are calculated independently for each threshold. The optimal threshold is selected post-hoc based on F1-score.

3. **Class Balance**: The algorithm does not account for class imbalance in metric calculation, though F1-score provides some balance.

### 9.3 Data Quality Requirements

1. **Minimum Sample Size**: Base-specific evaluation requires ≥10 positions per base type.

2. **Complete Structure Annotation**: The reference structure must cover all positions in the reactivity data (after alignment).

3. **Non-negative Reactivity**: The algorithm assumes non-negative reactivity values. Negative values may indicate experimental artifacts.

## 10. References and Related Work

The evaluation algorithms implement standard binary classification metrics commonly used in machine learning and bioinformatics:

- **ROC Analysis**: Fawcett, T. (2006). "An introduction to ROC analysis". Pattern Recognition Letters, 27(8), 861-874.

- **Precision-Recall Curves**: Davis, J., & Goadrich, M. (2006). "The relationship between Precision-Recall and ROC curves". ICML '06.

- **RNA Structure Probing**: Deigan, K. E., et al. (2009). "Accurate SHAPE-directed RNA structure determination". PNAS, 106(1), 97-102.

- **Secondary Structure Alignment**: Similar to sequence alignment algorithms but adapted for structure annotation alignment.

---

**Document Version**: 1.0  
**Last Updated**: 2024  
**Software Version**: ModDetector v1.0






