# Reactivity Score Normalization Method

Normalization of reactivity scores is a critical step in RNA structure probing analysis, serving to eliminate experimental batch effects and systematic biases, thereby enabling comparability of reactivity scores across different samples. This method implements independent normalization processing for two types of reactivity data: mutation signals and stop signals, supporting multiple normalization algorithms, and optionally applying piecewise linear mapping to further optimize the distribution characteristics of reactivity scores.

## Input Data Format

The input data is provided in CSV format, containing at least 6 columns. The first column contains the chromosome or reference sequence identifier (ChrID), which identifies the genomic region from which the data originates. The second column contains strand orientation information (Strand), which distinguishes between forward and reverse strand data, ensuring that data from different strands are processed independently. The fourth column contains the stop signal reactivity score (stop_reactivity), representing reactivity values detected through reverse transcription termination events, which reflect the secondary structure characteristics of RNA at that position. The fifth column contains the mutation signal reactivity score (mutation_reactivity), representing reactivity values detected through nucleotide mismatches, which reflect chemical modifications or structural accessibility of RNA at that position. Other columns may contain auxiliary data such as position information and base types, but the normalization process primarily utilizes the aforementioned four columns.

## Data Grouping Strategy

The normalization process employs a regional grouping strategy, grouping data according to chromosome identifiers and strand orientation to form independent processing units. Specifically, each processing unit is defined by a unique identifier in the format "ChrID|Strand", ensuring that all positions from the same chromosome and strand are grouped together. This grouping strategy ensures that the normalization process operates within biologically relevant regions, avoiding the influence of signal intensity differences between different regions on normalization results. For each group, the program separately extracts stop signal reactivity values and mutation signal reactivity values from all positions within that region, forming two independent data vectors that serve as input for subsequent normalization calculations.

## Mutation Signal Normalization

Mutation signal normalization processing is performed on the mutation reactivity value vector within each regional group. The program first extracts all mutation_reactivity values from positions within that region to form a numerical vector. Subsequently, according to the user-specified normalization method, normalization calculations are performed on this vector. The normalization process converts raw reactivity values into standardized reactivity scores, enabling comparability of data across different regions and samples. The normalized results are written back to the fifth column (mutation_reactivity column) of the original data structure, replacing the original raw reactivity values. This independent processing approach ensures that the normalization processes for mutation signals and stop signals do not interfere with each other, allowing optimization of their respective signal characteristics.

## Stop Signal Normalization

Stop signal normalization processing employs the same grouping strategy and normalization algorithm as mutation signals, but operates on the stop_reactivity value vector. The program extracts stop_reactivity values from all positions within each regional group, forming an independent numerical vector. The normalization calculation process is identical to mutation signal processing, but operates on a different data column. The normalized results are written back to the fourth column (stop_reactivity column) of the original data structure, replacing the original raw reactivity values. Through this independent normalization processing, stop signals and mutation signals can each obtain normalization results most suitable for their respective data distribution characteristics, thereby more accurately reflecting RNA structure information in subsequent analyses.

## Normalization Algorithms

### Percentile28 Method

The Percentile28 normalization method (also known as the 2-8% method) is based on percentile statistics, calculating the mean value of a specific percentile interval as the normalization factor. Specifically, for a given reactivity value vector, all values are first sorted in ascending order. If the total number of data points is less than 10, the original values are returned without normalization. Subsequently, the position indices for the 2nd percentile (P2) and 10th percentile (P10) are calculated, where the P2 position index is computed using the formula:

\[
p_2 = \lceil n \times 0.02 \rceil
\]

and the P10 position index is computed using:

\[
p_{10} = \lceil n \times 0.10 \rceil
\]

where \(n\) is the total number of data points. To ensure indices do not exceed array bounds, \(p_2\) and \(p_{10}\) are constrained to be no greater than \(n-1\). The normalization factor (avg) is obtained by calculating the mean of all values in the P2 to P10 interval:

\[
avg = \frac{1}{p_{10} - p_2} \sum_{i=p_2}^{p_{10}-1} sorted[i]
\]

where \(sorted[i]\) represents the value at the \(i\)-th position after sorting. If this mean value is less than or equal to 0, or if \(p_{10} \leq p_2\), a default value of 1.0 is used as the normalization factor. The final normalized reactivity score is calculated using the formula:

\[
R_{norm} = \frac{R_{raw}}{avg}
\]

where \(R_{raw}\) is the raw reactivity value and \(R_{norm}\) is the normalized reactivity score. This method effectively reduces the influence of extreme values on normalization results by using the mean of lower percentile intervals as the normalization baseline, while preserving the main distribution characteristics of the data.

### Winsor90 Method

The Winsor90 normalization method employs Winsorization (tail trimming) technology combined with percentile normalization. First, all reactivity values are sorted. If the total number of data points is less than 2, the original values are returned without normalization. The position indices for the 5th percentile (P5) and 95th percentile (P95) are then calculated, where the P5 position index is computed using:

\[
p_5 = \lfloor n \times 0.05 \rfloor
\]

and the P95 position index is computed using:

\[
p_{95} = \lceil n \times 0.95 \rceil
\]

where \(n\) is the total number of data points. To ensure indices do not exceed array bounds, \(p_5\) and \(p_{95}\) are constrained to be no greater than \(n-1\). The corresponding percentile values are extracted: \(v_5 = sorted[\min(p_5, n-1)]\) and \(v_{95} = sorted[\min(p_{95}, n-1)]\). Subsequently, all raw reactivity values undergo Winsorization processing, truncating values less than \(v_5\) to \(v_5\) and values greater than \(v_{95}\) to \(v_{95}\):

\[
R_{winsor} = \begin{cases} 
v_5 & \text{if } R_{raw} < v_5 \\ 
v_{95} & \text{if } R_{raw} > v_{95} \\ 
R_{raw} & \text{otherwise}
\end{cases}
\]

The normalization factor uses the value of \(v_{95}\), but ensures it is at least \(10^{-6}\) to avoid division by zero:

\[
norm\_val = \max(v_{95}, 10^{-6})
\]

The final normalized reactivity score is calculated using:

\[
R_{norm} = \frac{R_{winsor}}{norm\_val}
\]

This method effectively handles outliers in the data by truncating extreme values and normalizing based on the upper percentile, while maintaining the main distribution characteristics of the data.

### Boxplot Method

The Boxplot normalization method is based on boxplot statistical principles, identifying and removing extreme values before normalization. First, all reactivity values are sorted. If the total number of data points is less than 4, the original values are returned without normalization. The position indices for the first quartile (Q1) and third quartile (Q3) are then calculated, where the Q1 position index is computed using:

\[
q_1 = \lceil n \times 0.25 \rceil
\]

and the Q3 position index is computed using:

\[
q_3 = \lceil n \times 0.75 \rceil
\]

where \(n\) is the total number of data points. To ensure indices do not exceed array bounds, \(q_1\) and \(q_3\) are constrained to be no greater than \(n-1\). The corresponding quartile values are extracted: \(Q_1 = sorted[\min(q_1, n-1)]\) and \(Q_3 = sorted[\min(q_3, n-1)]\). The interquartile range (IQR) is calculated as:

\[
IQR = Q_3 - Q_1
\]

and the upper bound is calculated as:

\[
upper = Q_3 + 1.5 \times IQR
\]

All extreme values greater than the upper bound are removed, retaining data points less than or equal to the upper bound to form a filtered dataset \(filtered\). In the filtered dataset, the mean of the top 10% of data points is calculated as the normalization factor:

\[
avg = \frac{1}{top10} \sum_{i=0}^{top10-1} filtered[i]
\]

where \(top10\) is the number of top 10% data points, calculated using:

\[
top10 = \lceil |filtered| \times 0.10 \rceil
\]

and \(|filtered|\) represents the size of the filtered dataset. If this mean value is less than or equal to 0, or if \(top10 = 0\), a default value of 1.0 is used as the normalization factor. The final normalized reactivity score is calculated using:

\[
R_{norm} = \frac{R_{raw}}{avg}
\]

where the normalization calculation uses the original reactivity values rather than the filtered values. This method effectively handles outliers in the data by identifying them through boxplot principles and using the mean of high reactivity regions as the normalization baseline, while highlighting high reactivity signals.

## Piecewise Linear Mapping

After completing basic normalization, users can optionally apply the Zarringhalam piecewise linear mapping method to further transform the normalized reactivity scores. The Zarringhalam remap function is implemented entirely within the normalization workflow to ensure consistent behavior whenever Zarringhalam is used (with or without prior normalization). The function treats the minimum value (e.g., Winsor90's v5/v95) as 0, then maps the range [0, max] using the piecewise linear rule. This preserves relative mod/unmod differences for Siegfried reactivity calculation instead of clamping negatives to 0 and losing information.

**Key Implementation Details** (v0.14.10+):
- The function uses the actual maximum value instead of a fixed 1.0 for high-value region mapping
- This matches rf-norm implementation and prevents compression of high values that exceed 1.0 after normalization, which was causing high false positive rates
- After winsor90 normalization, values may exceed 1.0. Instead of normalizing to [0, 1] with fixed 1.0 as max, the function uses the actual max_val for the high-value region mapping [0.7, max] → [0.85, 1.0]

The specific mapping rules are as follows: for a normalized score \(x\), if \(x < 0.25\), the mapped score is:

\[
f(x) = \frac{x}{0.25} \times 0.35
\]

if \(0.25 \leq x < 0.3\), the mapped score is:

\[
f(x) = 0.35 + \frac{x - 0.25}{0.05} \times (0.55 - 0.35)
\]

if \(0.3 \leq x < 0.7\), the mapped score is:

\[
f(x) = 0.55 + \frac{x - 0.3}{0.4} \times (0.85 - 0.55)
\]

and if \(x \geq 0.7\), the mapped score uses the actual maximum value:

\[
f(x) = 0.85 + \frac{x - 0.7}{max - 0.7} \times (1.0 - 0.85)
\]

where \(max\) is the actual maximum value in the data (not fixed at 1.0). If \(max \leq 0.7\), all values in this range map to 0.85.

This piecewise linear mapping adjusts the distribution of reactivity scores, making the score distribution more uniform across low reactivity regions, medium reactivity regions, and high reactivity regions, thereby facilitating subsequent structure prediction and analysis. The use of actual maximum values ensures that high reactivity signals are not compressed, maintaining the full dynamic range of the data.

## Output Data Format

The output data maintains the same CSV format and column structure as the input data, but the values in the fourth column (stop_reactivity) and fifth column (mutation_reactivity) have been replaced with normalized reactivity scores. Normalized reactivity scores are output in floating-point format with 6 decimal places precision, formatted as "0.000000". If certain positions in the original data have empty or invalid reactivity values (such as NaN), the program preserves NaN values in the output to maintain data integrity. The first row of the output file is the header row, consistent with the input file header, and each subsequent row corresponds to a genomic position, containing all original information for that position along with the normalized reactivity scores.

## Configurable Parameters

The normalization method is selected through the `--method` parameter, supporting three methods: `percentile28` (or `2-8`, `2_8`), `winsor90` (or `winsorizing`), and `boxplot`. Users can select the most appropriate normalization method based on data characteristics and analysis requirements. The `--cutoff` parameter is used to set the SNP threshold, with a value range from 0.0 to 1.0 and a default value of 0.25. This parameter is used for parameter validation and future feature extensions. The `--bases` parameter is used to specify reactive base types, with a default value of "ACGT" indicating that all four base types participate in normalization calculations. Users can specify particular base types according to experimental design, for example, "AC" indicates that only adenine and cytosine reactivity are considered. In the current implementation, normalization processing is primarily based on a regional grouping strategy, grouping data according to chromosome identifiers and strand orientation, with all positions within each group normalized as a whole. The `--window` and `--window-offset` parameters are used for fixed window mode settings, where `--window` specifies the window size and `--window-offset` specifies the window offset. These parameters are primarily used for parameter validation in the current implementation, with actual normalization processing based on regional grouping rather than sliding windows. The `--dynamic` parameter is used to enable dynamic sliding window mode. When set to true, the program dynamically adjusts window size according to the distribution of reactive bases. This functionality is defined in the current implementation but not activated in the main normalization workflow. The `--window-size-dynamic` parameter is used to specify the minimum number of reactive bases required per window in dynamic windows, with a default value of 50 and a value range from 1 to 1000. This parameter is used in conjunction with dynamic window mode. The `--linear` parameter controls whether to apply Zarringhalam piecewise linear mapping. When set to true, normalized reactivity scores undergo piecewise linear mapping transformation, which is applied to all normalized reactivity scores after completing basic normalization. The `--log` parameter is used to specify the log file path. If not specified, the program defaults to creating a log file named "norm.log", recording detailed information and parameter settings for normalization processing, including software version, runtime, input and output file paths, normalization method, SNP threshold, base types, window parameters, linear mapping settings, and other key information.