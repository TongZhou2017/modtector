# Reactivity Score Calculation Methods

The reactivity score calculation is a critical step in RNA structure probing analysis, enabling the quantitative assessment of RNA structural accessibility and chemical modification patterns by comparing signal intensities between modified and unmodified samples. This method implements independent reactivity calculation for both mutation signals and stop signals, supporting multiple calculation algorithms that account for background signal correction and experimental variations. The reactivity calculation process transforms raw signal counts into normalized reactivity scores that reflect the structural and chemical properties of RNA molecules at each genomic position.

## Input Data Format

The input data for reactivity calculation consists of CSV format files containing signal count information from the ModDetector count command. The modified sample CSV file (required) and unmodified sample CSV file (optional) share the same column structure, with each row representing a single genomic position and its associated signal measurements. The first column contains the chromosome or reference sequence identifier (ChrID), which uniquely identifies the genomic region from which the data originates. The second column contains the strand orientation information (Strand), distinguishing between forward strand (`+`) and reverse strand (`-`) data to ensure proper handling of strand-specific signals. The third column contains the genomic position (Position), representing the 1-based coordinate after signal correction in the reference coordinate system. The fourth column contains the reference base (Base) at the corresponding genomic position, indicating the nucleotide type (A, C, G, or T) present in the reference sequence.

The fifth column contains the mutation signal count (rf_mutation_Count), representing the number of nucleotide mismatches detected between the reference genome and sequencing reads at this position. This count reflects the frequency of base substitutions, insertions, and deletions that occur during reverse transcription, which are indicative of chemical modifications or structural accessibility. The sixth column contains the stop signal count (pipe_truncation_count), representing the number of reverse transcription truncation events detected at this position. This count reflects the frequency of RT stops, which occur when reverse transcriptase encounters structural barriers or modified nucleotides that prevent further extension. The seventh column contains the total read depth (depth) at this position, representing the total number of sequencing reads that cover this genomic coordinate. This depth value is essential for calculating normalized signal rates, as it provides the denominator for converting raw counts into proportions.

Additional columns may be present in the input files, including insertion counts (rf_mutation_ins), deletion counts (rf_mutation_del), and individual base counts (base_A, base_C, base_G, base_T), but these are not directly used in the reactivity calculation process. The program extracts the mutation signal rate by dividing the mutation count by the depth (mutation = rf_mutation_Count / depth), and the stop signal rate by dividing the stop count by the depth (stop = pipe_truncation_count / depth), ensuring that signal intensities are normalized by sequencing coverage.

## Background Region Selection and k-Factor Calculation

Prior to reactivity calculation, the method identifies background regions that are expected to exhibit minimal structural constraints and chemical modifications, typically corresponding to ribosomal RNA regions (rRNA, 16S, or 18S) or regions specified through file naming conventions. These background regions serve as internal controls for estimating systematic biases between modified and unmodified samples, enabling the calculation of correction factors that account for experimental variations, sequencing depth differences, and technical artifacts.

For each background region position, the method extracts both stop signal rates and mutation signal rates from both the modified and unmodified samples. The stop signal rate at each position is calculated by dividing the stop signal count (pipe_truncation_count) by the read depth (depth) from the input CSV file. Similarly, the mutation signal rate is calculated by dividing the mutation signal count (rf_mutation_Count) by the read depth (depth). These calculated signal rates from background positions are collected into four independent vectors: stop_mod_bg (a vector containing stop signal rates calculated as pipe_truncation_count / depth from background positions in the modified sample), stop_unmod_bg (a vector containing stop signal rates from background positions in the unmodified sample), mut_mod_bg (a vector containing mutation signal rates calculated as rf_mutation_Count / depth from background positions in the modified sample), and mut_unmod_bg (a vector containing mutation signal rates from background positions in the unmodified sample).

The k-factor for stop signals (k_stop) is calculated as the ratio of the mean stop signal rate in the modified sample to the mean stop signal rate in the unmodified sample across all background positions. Specifically:

\[
k_{\text{stop}} = \frac{\text{mean}(\text{stop\_mod\_bg})}{\text{mean}(\text{stop\_unmod\_bg})}
\]

where stop_mod_bg is the vector of stop signal rates (calculated as pipe_truncation_count / depth) from background positions in the modified sample, and stop_unmod_bg is the corresponding vector from the unmodified sample. Similarly, the k-factor for mutation signals (k_mut) is calculated as:

\[
k_{\text{mut}} = \frac{\text{mean}(\text{mut\_mod\_bg})}{\text{mean}(\text{mut\_unmod\_bg})}
\]

where mut_mod_bg is the vector of mutation signal rates (calculated as rf_mutation_Count / depth) from background positions in the modified sample, and mut_unmod_bg is the corresponding vector from the unmodified sample. Note that these vectors (stop_mod_bg, stop_unmod_bg, mut_mod_bg, mut_unmod_bg) are intermediate computational variables derived from the input file features, not direct columns in the input CSV files.

These k-factors represent systematic scaling factors that account for differences in overall signal intensity between the two samples, which may arise from variations in sequencing depth, library preparation efficiency, or experimental conditions. When the unmodified sample mean is zero or when no background regions are identified, the k-factors default to 1.0, indicating no correction is applied. The k-factor correction approach ensures that reactivity scores reflect true biological differences rather than technical variations between samples.

## Stop Signal Reactivity Calculation

The stop signal reactivity calculation quantifies the structural accessibility of RNA molecules by measuring the frequency of reverse transcription truncation events. Three distinct methods are implemented to accommodate different experimental designs and analytical requirements, each with specific assumptions and computational approaches.

### Current Method (k-Factor Correction)

The current method employs k-factor correction to account for systematic differences between modified and unmodified samples. For each genomic position i, we define the following variables: T_i represents the stop signal rate in the modified sample, calculated as T_i = pipe_truncation_count_mod / depth_mod, where pipe_truncation_count_mod is the stop signal count (pipe_truncation_count) from the modified sample CSV file and depth_mod is the read depth (depth) from the modified sample. Similarly, U_i represents the stop signal rate in the unmodified sample, calculated as U_i = pipe_truncation_count_unmod / depth_unmod, where pipe_truncation_count_unmod and depth_unmod are the corresponding values (pipe_truncation_count and depth) from the unmodified sample CSV file. The stop signal reactivity (R_i^stop) is then calculated using the formula:

\[
R_i^{\text{stop}} = \max(0, T_i - k_{\text{stop}} \times U_i)
\]

where k_stop is the background-corrected scaling factor calculated from background regions as described above. The max(0, ...) operation ensures that reactivity scores are non-negative, as negative values would be biologically meaningless and likely represent measurement noise or artifacts. This method assumes that background regions exhibit similar structural properties in both samples, allowing the k-factor to capture and correct for technical variations while preserving true biological signals.

### Ding Method

The Ding method, based on the approach described by Ding et al. (2014), employs logarithmic transformation to stabilize variance and handle the wide dynamic range of signal intensities. For each genomic position i, we define n_Ti as the raw stop signal count (pipe_truncation_count) from the modified sample CSV file, and n_Ui as the corresponding raw stop signal count (pipe_truncation_count) from the unmodified sample CSV file. The method first applies logarithmic transformation to both stop signal counts after adding a pseudocount parameter (p) to avoid undefined logarithmic values at zero counts. The transformed signals are calculated as:

\[
T_i' = \ln(n_{Ti} + p)
\]

\[
U_i' = \ln(n_{Ui} + p)
\]

where p is the pseudocount parameter (default: 1.0) specified by the `--pseudocount` parameter. The reactivity score is then calculated as:

\[
R_i^{\text{stop}} = \max(0, T_i' - U_i')
\]

representing the logarithmic difference between the two samples. This logarithmic approach helps normalize the data distribution and reduces the influence of extreme values, making it particularly suitable for datasets with high variance or when dealing with low-coverage regions. The pseudocount parameter prevents division by zero and logarithmic errors, but should be chosen carefully as it can influence the final reactivity scores, especially for positions with very low signal counts.

### Rouskin Method

The Rouskin method represents a simplified approach that uses only the modified sample stop signal, without subtraction of the unmodified sample. For each genomic position i, we define T_i as the stop signal rate in the modified sample, calculated as T_i = pipe_truncation_count_mod / depth_mod, where pipe_truncation_count_mod is the stop signal count (pipe_truncation_count) and depth_mod is the read depth (depth) from the modified sample CSV file. The reactivity score is calculated directly as:

\[
R_i^{\text{stop}} = T_i
\]

representing the stop signal rate in the modified sample without any background subtraction or correction. This method is appropriate when unmodified control samples are not available or when the experimental design assumes that background stop signals are negligible. The Rouskin method provides a straightforward reactivity measure but does not account for systematic biases or background noise that might be present in the modified sample alone. This approach is particularly useful in mod-only experimental designs where unmodified samples cannot be obtained or when the focus is on relative reactivity patterns within a single sample rather than absolute reactivity differences between conditions.

## Mutation Signal Reactivity Calculation

The mutation signal reactivity calculation quantifies the chemical modification patterns and structural accessibility by measuring the frequency of nucleotide mismatches during reverse transcription. Three distinct methods are implemented, each designed to handle different aspects of mutation signal processing and background correction.

### Current Method (k-Factor Correction)

The current method for mutation signals employs the same k-factor correction approach as used for stop signals, ensuring consistency in the reactivity calculation framework. For each genomic position i, we define T_i as the mutation signal rate in the modified sample, calculated as T_i = rf_mutation_Count_mod / depth_mod, where rf_mutation_Count_mod is the mutation signal count (rf_mutation_Count) from the modified sample CSV file and depth_mod is the read depth (depth) from the modified sample. Similarly, U_i represents the mutation signal rate in the unmodified sample, calculated as U_i = rf_mutation_Count_unmod / depth_unmod, where rf_mutation_Count_unmod and depth_unmod are the corresponding values (rf_mutation_Count and depth) from the unmodified sample CSV file. The mutation signal reactivity (R_i^mut) is then calculated using the formula:

\[
R_i^{\text{mut}} = \max(0, T_i - k_{\text{mut}} \times U_i)
\]

where k_mut is the background-corrected scaling factor specific to mutation signals, calculated independently from k_stop as described in the background region selection section. The k_mut factor is calculated independently from k_stop, as mutation and stop signals may exhibit different background characteristics and systematic biases. This independent calculation allows the method to account for signal-specific variations while maintaining the same correction principle. The non-negativity constraint ensures that reactivity scores reflect true biological signals rather than measurement artifacts.

### Siegfried Method

The Siegfried method calculates mutation reactivity by directly computing the difference between mutation rates in the modified and unmodified samples, without applying k-factor correction. For each genomic position i, we define T_i as the mutation signal rate in the modified sample, calculated as T_i = rf_mutation_Count_mod / depth_mod, where rf_mutation_Count_mod is the mutation signal count (rf_mutation_Count) and depth_mod is the read depth (depth) from the modified sample CSV file. Similarly, U_i represents the mutation signal rate in the unmodified sample, calculated as U_i = rf_mutation_Count_unmod / depth_unmod, where rf_mutation_Count_unmod and depth_unmod are the corresponding values (rf_mutation_Count and depth) from the unmodified sample CSV file. The reactivity score is then calculated as:

\[
R_i^{\text{mut}} = T_i - U_i
\]

representing the direct difference between mutation rates without k-factor correction. **Negative values are allowed** (as in Shapemapper), where negative reactivity indicates that the background mutation rate is higher than the modified sample rate at that position. This approach assumes that the mutation rates are already normalized by sequencing depth and that systematic biases between samples are minimal or can be ignored. The Siegfried method provides a straightforward reactivity measure that directly reflects the difference in mutation frequencies, making it particularly suitable when sample-to-sample variations are well-controlled or when the focus is on detecting strong modification signals that are robust to background variations. This method is computationally simpler than the k-factor correction approach and may be preferred when background regions are not reliably identified or when the k-factor calculation is unreliable.

### Zubradt Method

The Zubradt method represents a simplified approach that uses only the modified sample mutation signal, without subtraction of the unmodified sample. For each genomic position i, we define T_i as the mutation signal rate in the modified sample, calculated as T_i = rf_mutation_Count_mod / depth_mod, where rf_mutation_Count_mod is the mutation signal count (rf_mutation_Count) and depth_mod is the read depth (depth) from the modified sample CSV file. The reactivity score is calculated directly as:

\[
R_i^{\text{mut}} = T_i
\]

representing the mutation signal rate in the modified sample without any background subtraction or correction. This method is appropriate for mod-only experimental designs where unmodified control samples are not available, or when the experimental protocol ensures that background mutation signals are negligible. The Zubradt method provides a direct measure of mutation frequency but does not account for systematic biases, sequencing errors, or background noise that might be present in the modified sample. This approach is particularly useful when the focus is on identifying highly modified regions within a single sample or when unmodified controls cannot be obtained due to experimental constraints.

## Mod-Only Mode

The reactivity calculation supports a mod-only mode when unmodified sample files are not provided. In this mode, the program processes only the modified sample data, and reactivity scores are calculated using methods that do not require unmodified sample subtraction (Rouskin method for stop signals and Zubradt method for mutation signals, or user-specified methods that can handle single-sample input). The k-factor calculation is skipped in mod-only mode, as there are no unmodified sample signals available for comparison. This mode is particularly useful for exploratory analyses, quality control assessments, or when experimental designs do not include unmodified control samples. The mod-only mode maintains the same output format and data structure as the standard two-sample mode, ensuring compatibility with downstream analysis pipelines.

## Output Data Format

The output of the reactivity calculation is stored in a CSV format file containing six columns that provide comprehensive reactivity information for each genomic position. The first column contains the chromosome or reference sequence identifier (ChrID), maintaining consistency with the input file format and enabling proper genomic coordinate tracking. The second column contains the strand orientation (Strand), preserving the strand-specific information from the input data. The third column contains the genomic position (Position), representing the 1-based coordinate in the reference coordinate system. The fourth column contains the reference base (Base), indicating the nucleotide type at the corresponding position.

The fifth column contains the stop signal reactivity score (reactivity_stop), representing the calculated reactivity value based on reverse transcription truncation events. This value is computed using the selected stop signal method (Current, Ding, or Rouskin) and reflects the structural accessibility of the RNA molecule at this position. Higher reactivity scores indicate greater structural accessibility or more frequent RT stops, which may correspond to single-stranded regions, modified nucleotides, or structural disruptions. The sixth column contains the mutation signal reactivity score (reactivity_mutation), representing the calculated reactivity value based on nucleotide mismatch frequencies. This value is computed using the selected mutation signal method (Current, Siegfried, or Zubradt) and reflects the chemical modification patterns or structural accessibility that lead to base substitutions during reverse transcription.

Both reactivity scores are output as floating-point numbers with six decimal places of precision (format: "0.000000"), ensuring sufficient numerical resolution for downstream analysis and visualization. Reactivity scores are constrained to be non-negative, with negative values from the calculation formulas being set to zero to reflect the biological constraint that reactivity cannot be negative. The output file maintains the same row structure as the input file, with each row corresponding to a single genomic position, enabling easy integration with downstream analysis tools and visualization software.

## Configurable Parameters

The reactivity calculation provides several configurable parameters that allow users to customize the calculation process according to their experimental design and analytical requirements. The `--mod-csv` parameter (required) specifies the path to the modified sample CSV file, which must exist and have a `.csv` file extension. This file contains the signal count data from the modified sample, typically generated by the ModDetector count command. The `--unmod-csv` parameter (optional) specifies the path to the unmodified sample CSV file, which follows the same format requirements as the modified sample file. When this parameter is not provided, the program operates in mod-only mode, using methods that do not require unmodified sample subtraction.

The `--output` parameter (required) specifies the path to the output CSV file where reactivity scores will be written. The output file must have a `.csv` file extension and will be created or overwritten during the calculation process. The `--stop-method` parameter allows users to select the calculation method for stop signal reactivity, with three available options: `current` (default) for k-factor correction, `ding` for logarithmic transformation method, and `rouskin` for mod-only approach. The `--mutation-method` parameter allows users to select the calculation method for mutation signal reactivity, with three available options: `current` (default) for k-factor correction, `siegfried` for direct difference calculation, and `zubradt` for mod-only approach.

The `--pseudocount` parameter is specific to the Ding method and sets the pseudocount value added to signal counts before logarithmic transformation. The default value is 1.0, and the parameter must be greater than 0. This parameter prevents undefined logarithmic values at zero counts and influences the final reactivity scores, especially for positions with very low signal counts. Users should choose this value based on their data characteristics, with smaller values providing more sensitivity to low signals but potentially amplifying noise, and larger values providing more stability but potentially reducing sensitivity.

The `--maxscore` parameter is defined for the Ding method but is not currently used in the implementation, serving as a placeholder for future functionality that may limit the upper bound of reactivity values. The default value is 10.0, and the parameter must be between 0 and 100. The `--threads` parameter specifies the number of parallel threads to use during data processing and reactivity calculation. If not specified, the program defaults to 8 threads. The thread count can be set to any positive integer, with values between 64 and 128 threads recommended for optimal performance on high-performance computing systems. The program uses parallel processing for data loading, signal extraction, reactivity calculation, and output file writing, significantly reducing computation time for large datasets.

The `--log` parameter (optional) specifies the path to a log file where detailed processing information will be recorded. If not specified, the program may use default logging behavior or suppress detailed logging. The log file contains information about software version, runtime timestamp, input and output file paths, selected calculation methods, parameter values, processing progress, and completion status. This logging information is valuable for reproducibility, debugging, and tracking analysis parameters in large-scale processing pipelines.

These configurable parameters provide flexibility in adapting the reactivity calculation to different experimental designs, data characteristics, and computational environments, while maintaining consistency in the core calculation framework and output format.
