    # Statistical Algorithms for Differential Site Detection

    ## Abstract

    This document provides a comprehensive description of the statistical algorithms implemented in the ModDetector software for detecting differential modification sites in RNA structure probing data. The algorithms operate at two distinct levels: single-nucleotide resolution and regional resolution. At the single-nucleotide level, seven statistical tests are available: Chi-square test, Continuity correction test, Student's t-test, Mann-Whitney U test, Wilcoxon signed-rank test, DiffScan method, and DeltaSHAPE method. At the regional level, two specialized methods are implemented: DiffScan (with scanning statistics) and DeltaSHAPE (with Z-factor and Z-score analysis). These methods enable robust detection of differentially modified nucleotides and structural regions across different experimental conditions, biological replicates, and sample groups.

    ## 1. Introduction

    RNA structure probing experiments generate quantitative reactivity data that reflect the accessibility of nucleotides and their modification states. Statistical comparison of reactivity values between experimental conditions is essential for identifying differentially modified sites. The statistical algorithms described herein are designed to compare reactivity data at both single-nucleotide and regional resolutions, accounting for experimental variability, biological replicates, and multiple comparison corrections.

    ## 2. Data Structure and Input Specification

    ### 2.1 Single-Nucleotide Level Data

    For single-nucleotide comparisons, the input consists of:

    - **Group 1 values**: `G₁ = {g₁₁, g₁₂, ..., g₁ₘ}` where `m ≥ 1` is the number of replicates or observations
    - **Group 2 values**: `G₂ = {g₂₁, g₂₂, ..., g₂ₙ}` where `n ≥ 1` is the number of replicates or observations
    - Each value `gᵢⱼ` represents a reactivity score, RFC (read fraction coverage), or PIPC (proportion of insertions and point mutations) value

    ### 2.2 Regional Level Data

    For regional comparisons, the input consists of:

    - **Series 1**: `S₁ = {(p₁₁, v₁₁), (p₁₂, v₁₂), ..., (p₁ₖ, v₁ₖ)}` where `p₁ᵢ` is position and `v₁ᵢ` is reactivity value
    - **Series 2**: `S₂ = {(p₂₁, v₂₁), (p₂₂, v₂₂), ..., (p₂ₖ, v₂ₖ)}` where positions correspond to Series 1
    - Both series are sorted by position and represent contiguous genomic regions (e.g., transcript regions)

    ### 2.3 Common Parameters

    - **Alpha (α)**: Significance level for hypothesis testing (default: 0.05)
    - **Min depth**: Minimum sequencing depth threshold for filtering low-coverage sites (optional)
    - **Min fold change**: Minimum fold change threshold for filtering (optional)

    ## 3. Single-Nucleotide Level Statistical Tests

    ### 3.1 Student's t-Test

    #### 3.1.1 Method Description

    The Student's t-test is a parametric statistical test used to compare the means of two independent groups, assuming that the data are approximately normally distributed. This test is appropriate when both groups have sufficient sample sizes (typically `n ≥ 2` for each group) and the data meet normality assumptions.

    #### 3.1.2 Mathematical Formulation

    **Null hypothesis**: `H₀: μ₁ = μ₂` (no difference between group means)

    **Alternative hypothesis**: `H₁: μ₁ ≠ μ₂` (two-tailed test)

    **Step 1: Sample Statistics**

    Calculate group means:
    ```
    μ̂₁ = (1/m) × Σᵢ₌₁ᵐ g₁ᵢ
    μ̂₂ = (1/n) × Σᵢ₌₁ⁿ g₂ᵢ
    ```

    Calculate group variances:
    ```
    σ̂₁² = (1/(m-1)) × Σᵢ₌₁ᵐ (g₁ᵢ - μ̂₁)²
    σ̂₂² = (1/(n-1)) × Σᵢ₌₁ⁿ (g₂ᵢ - μ̂₂)²
    ```

    **Step 2: Pooled Variance**

    For independent samples with equal variances (pooled t-test):
    ```
    σₚ² = ((m-1) × σ̂₁² + (n-1) × σ̂₂²) / (m + n - 2)
    ```

    **Step 3: Standard Error**

    ```
    SE = √(σₚ² × (1/m + 1/n))
    ```

    **Step 4: t-Statistic**

    ```
    t = (μ̂₁ - μ̂₂) / SE
    ```

    **Step 5: Degrees of Freedom**

    ```
    df = m + n - 2
    ```

    **Step 6: P-Value Calculation**

    The p-value is approximated using the standard normal distribution:
    ```
    p = 2 × (1 - Φ(|t|))
    ```
    where `Φ(x)` is the cumulative distribution function of the standard normal distribution.

    **Step 7: Effect Size (Cohen's d)**

    ```
    d = (μ̂₁ - μ̂₂) / √σₚ²
    ```

    #### 3.1.3 Input-Output Specification

    - **Input**:
    - Group 1: Vector of length `m ≥ 2` (returns non-significant result if `m < 2` or `n < 2`)
    - Group 2: Vector of length `n ≥ 2`
    - Alpha: Significance level (default: 0.05)

    - **Output**:
    - `statistic`: t-statistic value
    - `p_value`: Two-tailed p-value
    - `significant`: Boolean indicating if `p_value < alpha`
    - `effect_size`: Cohen's d (standardized mean difference)

    - **Parameters**:
    - No additional parameters required

    ### 3.2 Mann-Whitney U Test

    #### 3.2.1 Method Description

    The Mann-Whitney U test (also known as Wilcoxon rank-sum test) is a non-parametric statistical test that compares the distributions of two independent groups without assuming normality. It is based on ranks rather than raw values and is appropriate when data are not normally distributed or sample sizes are small.

    #### 3.2.2 Mathematical Formulation

    **Null hypothesis**: `H₀: P(G₁ > G₂) = P(G₂ > G₁)` (equal distributions)

    **Alternative hypothesis**: `H₁: P(G₁ > G₂) ≠ P(G₂ > G₁)` (two-tailed test)

    **Step 1: Pool and Rank Data**

    Combine both groups: `C = {g₁₁, ..., g₁ₘ, g₂₁, ..., g₂ₙ}`

    Sort `C` in ascending order and assign ranks. For tied values, assign the average rank.

    **Step 2: Calculate Rank Sums**

    ```
    R₁ = Σᵢ₌₁ᵐ rank(g₁ᵢ)
    R₂ = Σᵢ₌₁ⁿ rank(g₂ᵢ)
    ```

    **Step 3: Calculate U Statistics**

    ```
    U₁ = m × n + m × (m + 1) / 2 - R₁
    U₂ = m × n + n × (n + 1) / 2 - R₂
    ```

    **Step 4: U Statistic**

    ```
    U = min(U₁, U₂)
    ```

    **Step 5: Expected Value and Variance**

    Under the null hypothesis:
    ```
    E[U] = m × n / 2
    Var[U] = m × n × (m + n + 1) / 12
    ```

    **Step 6: Z-Statistic (Normal Approximation)**

    ```
    z = (U - E[U]) / √Var[U]
    ```

    **Step 7: P-Value**

    ```
    p = 2 × (1 - Φ(|z|))
    ```

    **Step 8: Effect Size**

    ```
    r = (U - E[U]) / (m × n)
    ```
    where `r` ranges from -1 to 1.

    #### 3.2.3 Input-Output Specification

    - **Input**:
    - Group 1: Vector of length `m ≥ 1` (returns non-significant if empty)
    - Group 2: Vector of length `n ≥ 1`
    - Alpha: Significance level

    - **Output**:
    - `statistic`: U-statistic value
    - `p_value`: Two-tailed p-value
    - `significant`: Boolean indicating significance
    - `effect_size`: Effect size `r`

    - **Parameters**:
    - No additional parameters required

    ### 3.3 Wilcoxon Signed-Rank Test

    #### 3.3.1 Method Description

    The Wilcoxon signed-rank test is a non-parametric test for paired samples or matched observations. It tests whether the median difference between paired observations is zero. This test is appropriate when comparing the same samples under different conditions.

    #### 3.3.2 Mathematical Formulation

    **Null hypothesis**: `H₀: median(G₁ - G₂) = 0` (median difference is zero)

    **Alternative hypothesis**: `H₁: median(G₁ - G₂) ≠ 0` (two-tailed test)

    **Prerequisite**: `m = n` (paired data)

    **Step 1: Calculate Differences**

    ```
    D = {d₁, d₂, ..., dₘ} where dᵢ = g₁ᵢ - g₂ᵢ
    ```

    Remove zero differences: `D' = {dᵢ | dᵢ ≠ 0}`

    **Step 2: Rank Absolute Differences**

    Sort absolute differences `|d'ᵢ|` in ascending order and assign ranks. For tied values, assign average ranks.

    **Step 3: Calculate Rank Sums**

    ```
    W₊ = Σ rank(d'ᵢ) for all d'ᵢ > 0
    W₋ = Σ rank(d'ᵢ) for all d'ᵢ < 0
    ```

    **Step 4: W Statistic**

    ```
    W = min(W₊, W₋)
    ```

    **Step 5: Expected Value and Variance**

    Under the null hypothesis (n' = number of non-zero differences):
    ```
    E[W] = n' × (n' + 1) / 4
    Var[W] = n' × (n' + 1) × (2n' + 1) / 24
    ```

    **Step 6: Z-Statistic**

    ```
    z = (W - E[W]) / √Var[W]
    ```

    **Step 7: P-Value**

    ```
    p = 2 × (1 - Φ(|z|))
    ```

    **Step 8: Effect Size**

    ```
    r = (W - E[W]) / (n' × (n' + 1) / 2)
    ```

    #### 3.3.3 Input-Output Specification

    - **Input**:
    - Group 1: Vector of length `m`
    - Group 2: Vector of length `n` where `n = m` (paired data)
    - Alpha: Significance level

    - **Output**:
    - `statistic`: W-statistic value
    - `p_value`: Two-tailed p-value
    - `significant`: Boolean indicating significance
    - `effect_size`: Effect size `r`

    - **Parameters**:
    - Requires paired data (`m = n`)

    ### 3.4 Chi-Square Test

    #### 3.4.1 Method Description

    The Chi-square test is applied to discrete data. For continuous reactivity data, the values are discretized into binary categories based on a threshold, and then a 2×2 contingency table is constructed for testing independence between groups and reactivity categories.

    #### 3.4.2 Mathematical Formulation

    **Step 1: Discretization**

    Calculate threshold:
    ```
    θ = (Σᵢ₌₁ᵐ g₁ᵢ + Σⱼ₌₁ⁿ g₂ⱼ) / (m + n)
    ```

    **Step 2: Construct Contingency Table**

    |              | Below/Equal θ | Above θ | Total |
    |--------------|---------------|---------|-------|
    | Group 1      | O₁₁           | O₁₂     | m     |
    | Group 2      | O₂₁           | O₂₂     | n     |
    | Total        | C₁            | C₂      | N     |

    where `Oᵢⱼ` are observed counts.

    **Step 3: Expected Frequencies**

    ```
    Eᵢⱼ = (Row_i × Col_j) / N
    ```

    **Step 4: Chi-Square Statistic**

    ```
    χ² = Σᵢ₌₁² Σⱼ₌₁² (Oᵢⱼ - Eᵢⱼ)² / Eᵢⱼ
    ```

    **Step 5: Degrees of Freedom**

    ```
    df = (2 - 1) × (2 - 1) = 1
    ```

    **Step 6: P-Value**

    ```
    p = 1 - F_χ²(χ², df=1)
    ```
    where `F_χ²` is the chi-square cumulative distribution function.

    #### 3.4.3 Input-Output Specification

    - **Input**:
    - Group 1: Vector of length `m ≥ 1`
    - Group 2: Vector of length `n ≥ 1`
    - Alpha: Significance level

    - **Output**:
    - `statistic`: Chi-square statistic
    - `p_value`: P-value
    - `significant`: Boolean indicating significance
    - `effect_size`: Not computed for this test

    - **Parameters**:
    - No additional parameters

    ### 3.5 Continuity Correction Test

    #### 3.5.1 Method Description

    The Continuity correction test (Yates' correction) is a modification of the Chi-square test that applies a correction factor to improve the approximation to the chi-square distribution, especially for small sample sizes.

    #### 3.5.2 Mathematical Formulation

    This test uses the Chi-square test formulation (Section 3.4) but applies Yates' correction:

    **Step 1-4**: Same as Chi-square test (Section 3.4.2)

    **Step 5: Apply Continuity Correction**

    ```
    χ²_corrected = max(0, χ² - 0.5)
    ```

    **Step 6: P-Value**

    ```
    p = 1 - F_χ²(χ²_corrected, df=1)
    ```

    #### 3.5.3 Input-Output Specification

    - **Input**: Same as Chi-square test (Section 3.4.3)
    - **Output**: Same as Chi-square test, with corrected statistic and p-value
    - **Parameters**: No additional parameters

    ### 3.6 DiffScan Method (Single-Nucleotide Level)

    #### 3.6.1 Method Description

    The DiffScan method at the single-nucleotide level is designed for comparing normalized reactivity values between two conditions. It uses a simplified approach that directly compares reactivity differences with predefined thresholds, suitable for cases where data are already normalized and ready for comparison.

    #### 3.6.2 Mathematical Formulation

    **Case 1: Single Value Comparison (`m = 1, n = 1`)**

    ```
    normalized_diff = |g₁₁ - g₂₁|
    ```

    Significance criteria:
    ```
    significant = (normalized_diff > 0.3) AND 
                (g₁₁ > 0.2 OR g₂₁ > 0.2) AND
                (g₁₁ > 0.05 AND g₂₁ > 0.05)
    ```

    P-value assignment:
    ```
    p = { 0.001  if normalized_diff > 0.6
        { 0.01   if normalized_diff > 0.4
        { 0.05   if normalized_diff > 0.3 AND significant
        { 0.5    otherwise
    ```

    **Case 2: Multiple Values (`m > 1` or `n > 1`)**

    Falls back to Wilcoxon signed-rank test (Section 3.3) if paired, or Mann-Whitney U test (Section 3.2) if unpaired.

    #### 3.6.3 Input-Output Specification

    - **Input**:
    - Group 1: Vector of reactivity values (typically normalized)
    - Group 2: Vector of reactivity values (typically normalized)
    - Alpha: Significance level

    - **Output**:
    - `statistic`: Absolute difference (single value case) or Wilcoxon/Mann-Whitney statistic
    - `p_value`: Assigned or computed p-value
    - `significant`: Boolean based on threshold criteria or statistical test
    - `effect_size`: Absolute difference or effect size from underlying test

    - **Parameters**:
    - Hard-coded thresholds: 0.3 (difference), 0.2 (minimum reactivity), 0.05 (minimum both values)

    ### 3.7 DeltaSHAPE Method (Single-Nucleotide Level)

    #### 3.7.1 Method Description

    At the single-nucleotide level, the DeltaSHAPE method uses the same logic as the DiffScan method (Section 3.6). The DeltaSHAPE-specific Z-factor and Z-score calculations are applied only at the regional level (Section 4.2).

    #### 3.7.2 Mathematical Formulation

    Same as DiffScan method (Section 3.6.2).

    #### 3.7.3 Input-Output Specification

    Same as DiffScan method (Section 3.6.3).

    ## 4. Regional Level Statistical Methods

    ### 4.1 DiffScan Method (Regional Level)

    #### 4.1.1 Method Description

    The DiffScan method at the regional level combines sliding-window Wilcoxon tests with scanning statistics to identify differentially modified regions. It uses a Monte Carlo permutation approach to control the family-wise error rate (FWER) across multiple region tests.

    #### 4.1.2 Mathematical Formulation

    **Input**: 
    - Series 1: `S₁ = {(p₁, v₁₁), (p₂, v₁₂), ..., (pₖ, v₁ₖ)}`
    - Series 2: `S₂ = {(p₁, v₂₁), (p₂, v₂₂), ..., (pₖ, v₂ₖ)}`
    - Window radius: `r` (default: 2)
    - Min region length: `L_min` (default: 5)
    - Max region length: `L_max` (default: 50)
    - Monte Carlo iterations: `M` (default: 200)

    **Step 1: Compute Window P-Values**

    For each position `pᵢ`:
    - Extract window: `W₁ᵢ = {v₁ⱼ | |j - i| ≤ r}`
    - Extract window: `W₂ᵢ = {v₂ⱼ | |j - i| ≤ r}`
    - Perform Wilcoxon signed-rank test (Section 3.3) on `W₁ᵢ` vs `W₂ᵢ`
    - Obtain p-value: `pᵢ`

    Result: `P = {(p₁, p₁), (p₂, p₂), ..., (pₖ, pₖ)}`

    **Step 2: Compute Scanning Statistic Q(R)**

    For each potential region `R` of length `L` where `L_min ≤ L ≤ L_max`:

    Compute prefix sums:
    ```
    prefix[0] = 0
    prefix[i] = prefix[i-1] - ln(pᵢ)   for i = 1, 2, ..., k
    ```

    For region from index `i` to `i+L-1`:
    ```
    S = prefix[i+L] - prefix[i]
    Q(R) = S / √L
    ```

    **Step 3: Monte Carlo Threshold Calculation**

    For iteration `m = 1, 2, ..., M`:
    1. Randomly permute p-values: `P^(m) = shuffle(P)`
    2. Compute `Q_max^(m) = max(Q(R) for all regions R in permuted data)`

    Sort `{Q_max^(1), Q_max^(2), ..., Q_max^(M)}` and find threshold:
    ```
    Q_threshold = Q_max^((1-α) × M)
    ```

    **Step 4: Identify Significant Regions**

    Region `R` is significant if:
    ```
    Q(R) ≥ Q_threshold
    ```

    #### 4.1.3 Input-Output Specification

    - **Input**:
    - Series 1 and Series 2: Position-reactivity pairs for contiguous regions
    - Window radius: `r` (default: 2)
    - Min region length: `L_min` (default: 5)
    - Max region length: `L_max` (default: 50)
    - Monte Carlo iterations: `M` (default: 200)
    - Alpha: Significance level

    - **Output**:
    - Regions file with columns: ChrID, Strand, Start, End, Length, MaxQ, MonteCarlo_Threshold, Significant
    - Each significant region represents a contiguous differentially modified region

    - **Parameters**:
    - `window_radius` (default: 2): Half-width of sliding window
    - `min_len` (default: 5): Minimum region length in nucleotides
    - `max_len` (default: 50): Maximum region length in nucleotides
    - `mc_iterations` (default: 200): Number of Monte Carlo permutations

    ### 4.2 DeltaSHAPE Method (Regional Level)

    #### 4.2.1 Method Description

    The DeltaSHAPE method identifies differentially modified regions using Z-factor and Z-score calculations combined with a sliding window approach. It is designed to detect regions with consistent reactivity differences between conditions.

    #### 4.2.2 Mathematical Formulation

    **Input**: Same as DiffScan regional method (Section 4.1.2)

    **Step 1: Compute Position-Level Statistics**

    For each position `pᵢ`:
    - Delta reactivity: `Δᵢ = v₁ᵢ - v₂ᵢ`
    - Error estimates: 
    ```
    ε₁ᵢ = v₁ᵢ × 0.1  (assumed 10% error)
    ε₂ᵢ = v₂ᵢ × 0.1
    ```
    - Z-factor:
    ```
    Z_factorᵢ = { 1 - (1.96 × (ε₁ᵢ + ε₂ᵢ) / |Δᵢ|)   if |Δᵢ| > 0
                { 0                                    otherwise
    ```
    - Z-score:
    ```
    Z_scoreᵢ = { |Δᵢ| / 0.1   if |Δᵢ| > 0
                { 0             otherwise
    ```
    where 0.1 is an assumed standard deviation.

    Result: `Stats = {(p₁, Z_factor₁, Z_score₁), ..., (pₖ, Z_factorₖ, Z_scoreₖ)}`

    **Step 2: Identify Significant Positions**

    Position `pᵢ` is significant if:
    ```
    Z_factorᵢ > 0 AND |Z_scoreᵢ| ≥ 1.0
    ```

    **Step 3: Merge into Regions**

    Starting from position `i = 0`:
    - Find consecutive significant positions
    - Extend region while:
    - Region length < `L_max`
    - Current position is significant OR (region has significant positions AND current is near-significant)
    - Record region if:
    - Region length `≥ L_min`
    - At least half of positions in region are significant

    For each region, compute:
    ```
    Z_factor_max = max(Z_factorᵢ for i in region)
    Z_score_max = max(|Z_scoreᵢ| for i in region)
    ```

    #### 4.2.3 Input-Output Specification

    - **Input**:
    - Series 1 and Series 2: Position-reactivity pairs
    - Window radius: `r` (used for smoothing, default: 2)
    - Min region length: `L_min` (default: 5)
    - Max region length: `L_max` (default: 50)

    - **Output**:
    - Regions file with columns: ChrID, Strand, Start, End, Length, Z_factor, Z_score, Significant
    - Each region contains at least `L_min/2` significant positions

    - **Parameters**:
    - Error percentage: 10% (hard-coded)
    - Standard deviation: 0.1 (hard-coded)
    - Z-factor threshold: > 0
    - Z-score threshold: ≥ 1.0
    - Min region length: `L_min` (default: 5)
    - Max region length: `L_max` (default: 50)

    ## 5. Statistical Function Implementations

    ### 5.1 Normal Distribution Approximation

    The cumulative distribution function (CDF) of the standard normal distribution is approximated using the error function:

    ```
    Φ(x) = 0.5 × (1 + erf(x / √2))
    ```

    where the error function `erf(x)` is computed using a polynomial approximation:

    ```
    erf(x) ≈ sign(x) × y
    ```

    with `y = 1 - (((((a₅×t + a₄)×t + a₃)×t + a₂)×t + a₁)×t × exp(-x²))`

    where:
    - `t = 1 / (1 + p×|x|)`
    - `p = 0.3275911`
    - `a₁ = 0.254829592, a₂ = -0.284496736, a₃ = 1.421413741, a₄ = -1.453152027, a₅ = 1.061405429`

    ### 5.2 Chi-Square Distribution Approximation

    The chi-square distribution CDF is approximated using a normal approximation:

    ```
    F_χ²(x, df) ≈ Φ((x - df) / √(2×df))
    ```

    ## 6. Additional Metrics and Scores

    ### 6.1 Fold Change Calculation

    Fold change is calculated as the ratio of group means:

    ```
    FC = { μ̂₁ / μ̂₂              if μ̂₂ > 0
        { +∞                   if μ̂₁ > 0 AND μ̂₂ = 0
        { 1.0                  if μ̂₁ = 0 AND μ̂₂ = 0
    ```

    For RFC and PIPC values, separate fold changes are computed:
    - RFC fold change: `FC_RFC = RFC₁ / RFC₂`
    - PIPC fold change: `FC_PIPC = PIPC₁ / PIPC₂`

    ### 6.2 Modification Score

    The modification score is a composite metric that combines RFC and PIPC fold changes:

    ```
    Modification_Score = (FC_RFC + FC_PIPC) / 2.0
    ```

    This score provides a unified measure of modification intensity, ranging from 0 to +∞ (typically 0 to 10 for well-behaved data). A score > 1.0 indicates increased modification in group 1 compared to group 2.

    ## 7. Output Format Specifications

    ### 7.1 Single-Nucleotide Level Output

    The output CSV file contains the following columns:

    - **ChrID**: Chromosome/transcript identifier
    - **Strand**: Strand orientation (+ or -)
    - **Position**: Genomic position
    - **RefBase**: Reference nucleotide (if available)
    - **Group1_Mean_Reactivity**: Mean reactivity for group 1 (or Mod/Reactivity1)
    - **Group1_Std_Reactivity**: Standard deviation for group 1
    - **Group2_Mean_Reactivity**: Mean reactivity for group 2 (or Unmod/Reactivity2)
    - **Group2_Std_Reactivity**: Standard deviation for group 2
    - **FoldChange**: Fold change (Group1/Group2)
    - **Test_Statistic**: Calculated test statistic
    - **P_Value**: P-value from statistical test
    - **Significant**: TRUE/FALSE indicating significance
    - **Effect_Size**: Effect size (if applicable)

    ### 7.2 Regional Level Output

    The output CSV file (`.regions.csv`) contains:

    - **ChrID**: Chromosome/transcript identifier
    - **Strand**: Strand orientation
    - **Start**: Start position of region
    - **End**: End position of region
    - **Length**: Region length in nucleotides
    - **MaxQ** (DiffScan) or **Z_factor** (DeltaSHAPE): Region statistic
    - **MonteCarlo_Threshold** (DiffScan) or **Z_score** (DeltaSHAPE): Threshold or score
    - **Significant**: TRUE/FALSE indicating significance

    ## 8. Performance Considerations

    ### 8.1 Computational Complexity

    - **Single-nucleotide tests**: O(m log m + n log n) for ranking-based tests, O(m + n) for parametric tests
    - **Regional DiffScan**: O(k × (2r+1) × log(2r+1) + M × k × L_max) where k is number of positions
    - **Regional DeltaSHAPE**: O(k × L_max) for region merging

    ### 8.2 Memory Requirements

    - Single-nucleotide tests: O(m + n)
    - Regional methods: O(k + M) where M is Monte Carlo iterations

    ### 8.3 Recommended Usage

    - **Small datasets** (< 10,000 positions): All methods suitable
    - **Medium datasets** (10,000 - 100,000 positions): All methods suitable with appropriate filtering
    - **Large datasets** (> 100,000 positions): Consider filtering by depth/fold change before statistical testing

    ## 9. References and Implementation Notes

    ### 9.1 Statistical Tests

    - Student's t-test: Standard independent samples t-test
    - Mann-Whitney U test: Non-parametric alternative to t-test
    - Wilcoxon signed-rank test: Paired non-parametric test
    - Chi-square test: Contingency table test for independence
    - Continuity correction: Yates' correction for chi-square

    ### 9.2 Regional Methods

    - DiffScan: Based on scanning statistics with Monte Carlo FWER control
    - DeltaSHAPE: Based on Z-factor methodology adapted from high-throughput screening

    ### 9.3 Implementation Details

    - All statistical tests are implemented in Rust
    - P-values are computed using approximations suitable for large-scale genomic applications
    - Effect sizes are computed where applicable to provide biological interpretation
    - Regional methods are designed to handle transcript-level comparisons

