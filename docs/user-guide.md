# User Guide

This comprehensive guide covers all aspects of using modtector, from basic concepts to advanced features.

## Table of Contents

1. [Understanding RNA Modifications](#understanding-rna-modifications)
2. [modtector Workflow](#moddetector-workflow)
3. [Data Preparation](#data-preparation)
4. [Signal Types](#signal-types)
5. [Normalization Methods](#normalization-methods)
6. [Reactivity Calculation](#reactivity-calculation)
7. [Evaluation Metrics](#evaluation-metrics)
8. [Visualization](#visualization)
9. [Advanced Features](#advanced-features)
10. [Best Practices](#best-practices)

## Understanding RNA Modifications

### What are RNA Modifications?

RNA modifications are chemical alterations to RNA nucleotides that can affect:
- RNA structure and stability
- Protein-RNA interactions
- Translation efficiency
- RNA localization
- Gene expression regulation

### Common RNA Modifications

- **m6A**: N6-methyladenosine (most common)
- **m1A**: N1-methyladenosine
- **m5C**: 5-methylcytosine
- **Ψ**: Pseudouridine
- **D**: Dihydrouridine

### Detection Methods

modtector supports detection based on:
- **Stop signals**: Pipeline truncation during reverse transcription
- **Mutation signals**: Base mutations during reverse transcription

## modtector Workflow

### Overview

modtector follows a systematic workflow:

```
Raw BAM Files → Pileup Analysis → Normalization → Reactivity Calculation → Evaluation
     ↓              ↓                ↓                    ↓                ↓
  Alignment      Signal Counts    Filtered Data    Modification Scores   Accuracy
```

### Step-by-Step Process

1. **Data Input**: BAM files, reference sequences, structure files
2. **Pileup Analysis**: Count stop and mutation signals
3. **Normalization**: Filter noise and outliers
4. **Reactivity Calculation**: Compare modified vs unmodified samples
5. **Visualization**: Generate plots and charts
6. **Evaluation**: Assess accuracy using known structures

## Data Preparation

### Input Requirements

#### BAM Files
- **Format**: BAM format with proper alignment
- **Quality**: High-quality alignments with minimal mismatches
- **Coverage**: Sufficient read depth (recommended >50x)
- **Paired samples**: Modified and unmodified samples

#### Reference Sequences
- **Format**: FASTA format
- **Quality**: High-quality reference sequences
- **Matching**: Must match your BAM file alignments
- **Completeness**: Include all target regions

#### Secondary Structure Files
- **Format**: Dot-bracket notation (.dp files)
- **Source**: Experimentally determined or predicted structures
- **Accuracy**: High-confidence structures for evaluation
- **Coverage**: Cover all regions of interest

### Data Quality Checks

Before running modtector, verify:

1. **BAM File Quality**:
   ```bash
   samtools flagstat sample.bam
   samtools view -c sample.bam
   ```

2. **Reference Sequence**:
   ```bash
   samtools faidx reference.fa
   ```

3. **Alignment Quality**:
   ```bash
   samtools view sample.bam | head -100
   ```

4. **BAM Index Files** (Required for count command):
   ```bash
   # Check if index exists
   ls -lh sample.bam.bai
   
   # Create index if missing
   samtools index -b sample.bam -o sample.bam.bai -@ 8
   ```

### Batch Processing and Single-cell Mode

#### Batch Mode (`--batch`)

Batch mode allows you to process multiple BAM files sequentially, each file independently.

**Use Cases:**
- Processing multiple samples that need separate analysis
- Bulk RNA-seq data with multiple replicates
- When each file should be treated as an independent sample

**Example:**
```bash
modtector count --batch \
  -b "/path/to/data/*sort.bam" \
  -f reference.fa \
  -o output_dir/ \
  -t 8 \
  -w 10000
```

**Requirements:**
- BAM files must be sorted and indexed (`.bai` files must exist)
- Glob pattern must match at least one file
- Output path must be a directory

#### Single-cell Unified Mode (`--single-cell`)

Single-cell unified mode is optimized for single-cell RNA-seq data, providing 2-3x performance improvement through unified processing.

**Key Features:**
- Unified data distribution scanning (once instead of per-file)
- Automatic cell label extraction from filenames
- Cross-file parallel processing
- Reduced I/O overhead

**Cell Label Extraction:**
- Supports RHX pattern: `*RHX672.sort.bam` → `RHX672`
- Falls back to last underscore-separated part: `sample_cell123.bam` → `cell123`
- Final fallback: base filename without extension

**Example:**
```bash
modtector count --single-cell \
  -b "/path/to/single_cell/*sort.bam" \
  -f reference.fa \
  -o output_dir/ \
  -t 8 \
  -w 10000 \
  -l batch.log
```

**Output:**
- Each cell generates a separate CSV file: `RHX672.csv`, `RHX673.csv`, etc.
- Log files are also generated per cell

**Performance Comparison:**
- Batch mode: ~N × single_file_time (sequential processing)
- Single-cell unified mode: ~(N × single_file_time) / 2-3 (unified processing)

## Signal Types

### Stop Signals

Stop signals occur when reverse transcription is truncated at modification sites.

#### Characteristics
- **Detection**: Pipeline truncation events
- **Sensitivity**: High sensitivity for certain modifications
- **Specificity**: Good specificity with proper controls
- **Background**: Low background noise

#### Analysis
```bash
modtector count -b sample.bam -f reference.fa -o output.csv
```

### Mutation Signals

Mutation signals occur when reverse transcription introduces errors at modification sites.

#### Characteristics
- **Detection**: Base mutation events
- **Sensitivity**: Moderate sensitivity
- **Specificity**: Good specificity with proper controls
- **Background**: Moderate background noise

#### Analysis
```bash
modtector count -b sample.bam -f reference.fa -o output.csv
```

### Combined Analysis

modtector can analyze both signal types simultaneously:

```bash
modtector reactivity -M mod.csv -U unmod.csv -O reactivity.csv -t both
```

## Normalization Methods

### Purpose

Normalization removes systematic biases and noise from the data:
- **Technical noise**: Sequencing artifacts
- **Biological noise**: Background signal
- **Systematic bias**: Sample preparation effects

### Available Methods

#### 1. Percentile28 Normalization
- **Method**: 28th percentile scaling
- **Use case**: General purpose normalization
- **Advantages**: Robust to outliers
- **Disadvantages**: May be conservative

```bash
modtector norm -i input.csv -o output.csv -m percentile28
```

#### 2. Winsor90 Normalization
- **Method**: 90th percentile winsorization
- **Use case**: High-quality data
- **Advantages**: Preserves signal distribution
- **Disadvantages**: Sensitive to outliers

```bash
modtector norm -i input.csv -o output.csv -m winsor90
```

#### 3. Boxplot Normalization
- **Method**: Boxplot-based outlier removal
- **Use case**: Data with many outliers
- **Advantages**: Effective outlier removal
- **Disadvantages**: May remove valid signals

```bash
modtector norm -i input.csv -o output.csv -m boxplot
```

### Window-Based Normalization

#### Fixed Windows
```bash
modtector norm -i input.csv -o output.csv -m winsor90 --window 1000
```

#### Dynamic Windows
```bash
modtector norm -i input.csv -o output.csv -m winsor90 --dynamic
```

#### Sliding Windows
```bash
modtector norm -i input.csv -o output.csv -m winsor90 --window 500 --window-offset 100
```

## Reactivity Calculation

### Purpose

Reactivity calculation quantifies the difference between modified and unmodified samples to identify modification sites.

### Calculation Methods

#### Stop Signal Methods

1. **Current Method**:
   ```bash
   modtector reactivity -M mod.csv -U unmod.csv -O output.csv -s current
   ```

2. **Ding Method**:
   ```bash
   modtector reactivity -M mod.csv -U unmod.csv -O output.csv -s ding
   ```

3. **Rouskin Method**:
   ```bash
   modtector reactivity -M mod.csv -U unmod.csv -O output.csv -s rouskin
   ```

#### Mutation Signal Methods

1. **Current Method**:
   ```bash
   modtector reactivity -M mod.csv -U unmod.csv -O output.csv -m current
   ```

2. **Siegfried Method**:
   ```bash
   modtector reactivity -M mod.csv -U unmod.csv -O output.csv -m siegfried
   ```

3. **Zubradt Method**:
   ```bash
   modtector reactivity -M mod.csv -U unmod.csv -O output.csv -m zubradt
   ```

### Parameters

#### Pseudocount
- **Purpose**: Avoid zero values in logarithmic calculations
- **Default**: 1
- **Range**: 0.1 - 10
- **Recommendation**: Use 0.5 for sparse data

```bash
modtector reactivity -M mod.csv -U unmod.csv -O output.csv --pseudocount 0.5
```

#### Maximum Score
- **Purpose**: Limit upper bound of reactivity values
- **Default**: 10
- **Range**: 5 - 50
- **Recommendation**: Use 5 for conservative analysis

```bash
modtector reactivity -M mod.csv -U unmod.csv -O output.csv --maxscore 5
```

## Evaluation Metrics

### Purpose

Evaluation metrics assess the accuracy of modification detection using known secondary structures.

### Available Metrics

#### 1. Area Under the Curve (AUC)
- **Range**: 0 - 1
- **Interpretation**: 
  - > 0.8: Excellent
  - 0.7 - 0.8: Good
  - 0.6 - 0.7: Fair
  - < 0.6: Poor

#### 2. F1-Score
- **Range**: 0 - 1
- **Interpretation**:
  - > 0.8: Excellent
  - 0.6 - 0.8: Good
  - 0.4 - 0.6: Fair
  - < 0.4: Poor

#### 3. Accuracy
- **Range**: 0 - 1
- **Interpretation**:
  - > 0.9: Excellent
  - 0.8 - 0.9: Good
  - 0.7 - 0.8: Fair
  - < 0.7: Poor

#### 4. Sensitivity (Recall)
- **Range**: 0 - 1
- **Interpretation**: Proportion of true modifications detected

#### 5. Specificity
- **Range**: 0 - 1
- **Interpretation**: Proportion of true non-modifications correctly identified

#### 6. Positive Predictive Value (PPV)
- **Range**: 0 - 1
- **Interpretation**: Proportion of predicted modifications that are true

#### 7. Negative Predictive Value (NPV)
- **Range**: 0 - 1
- **Interpretation**: Proportion of predicted non-modifications that are true

### Evaluation Process

```bash
modtector evaluate -r reactivity.csv -s structure.dp -o results/ -g gene_id
```

### Output Files

1. **Comprehensive Results**: `*_comprehensive.txt`
2. **ROC Curves**: `*_roc.svg`
3. **PR Curves**: `*_pr.svg`
4. **Combined Plots**: `*_combined_roc.svg`

## Visualization

### Purpose

Visualization helps interpret results and identify patterns in the data.

### Available Plots

#### 1. Signal Distribution Plots
- **Purpose**: Show signal distribution across positions
- **Format**: SVG
- **Content**: Stop and mutation signals

#### 2. Reactivity Plots
- **Purpose**: Display reactivity scores
- **Format**: SVG
- **Content**: Modification sites and scores

#### 3. ROC Curves
- **Purpose**: Show classification performance
- **Format**: SVG
- **Content**: True positive rate vs false positive rate

#### 4. PR Curves
- **Purpose**: Show precision-recall performance
- **Format**: SVG
- **Content**: Precision vs recall

#### 5. Comparison Plots
- **Purpose**: Compare different samples or methods
- **Format**: SVG

#### 6. RNA Structure SVG Plots
- **Purpose**: Visualize reactivity data on RNA secondary structure
- **Format**: SVG
- **Content**: Colored circles mapped to structure positions
- **Features**: 
  - Multi-signal support (stop, mutation, etc.)
  - Strand selection (+, -, both)
  - Base filtering (A, T, C, G)
  - Alignment support
  - Color-coded reactivity scores
- **Content**: Side-by-side comparisons

#### 7. Interactive HTML Visualizations
- **Purpose**: Interactive web-based visualization of reactivity data on RNA structure
- **Format**: HTML (with embedded SVG and JavaScript)
- **Content**: Interactive RNA structure with reactivity data overlay
- **Features**:
  - **Hover Tooltips**: Display position, base type, and reactivity values on hover
  - **Zoom and Pan**: Mouse wheel zoom and click-drag panning
  - **Filtering Controls**:
    - Reactivity threshold filtering (min/max sliders)
    - Base type filtering (A, T, C, G, or all)
  - **Highlighting**: Click circles to highlight specific positions
  - **Export**: Export the current view as SVG
  - **Reset View**: Reset zoom and filters to default
- **Usage**: Add `--interactive` flag to the plot command to generate HTML instead of static SVG

### SVG Plotting

#### Overview
RNA structure SVG plots provide an intuitive way to visualize reactivity data directly on the RNA secondary structure. This helps researchers identify modification sites and understand the relationship between structure and reactivity.

#### SVG-Only Mode
When you only need SVG plots without regular distribution plots, use the simplified command:

```bash
modtector plot \
    -o svg_output/ \
    --svg-template rna_structure.svg \
    --reactivity reactivity_data.csv
```

#### Multi-Signal SVG Plotting
For data with multiple signal types (stop, mutation), plot all signals:

```bash
modtector plot \
    -o svg_output/ \
    --svg-template rna_structure.svg \
    --reactivity reactivity_data.csv \
    --svg-bases ATGC \
    --svg-signal all \
    --svg-strand +
```

#### Interactive HTML Visualization
Generate interactive HTML visualizations with zoom, pan, filtering, and tooltip features:

```bash
modtector plot \
    -o interactive_output/ \
    --svg-template rna_structure.svg \
    --reactivity reactivity_data.csv \
    --interactive \
    --svg-signal stop \
    --svg-bases ACGT
```

The interactive HTML file can be opened in any modern web browser and provides:
- **Mouse wheel zoom**: Scroll to zoom in/out
- **Click and drag**: Pan around the structure
- **Hover tooltips**: See position, base, and reactivity values
- **Filtering**: Adjust reactivity thresholds and filter by base type
- **Highlighting**: Click circles to highlight specific positions
- **Export**: Download the current view as SVG

#### Single Signal SVG Plotting
For single signal data, specify the signal type:

```bash
modtector plot \
    -o svg_output/ \
    --svg-template rna_structure.svg \
    --reactivity reactivity_data.csv \
    --svg-bases AC \
    --svg-signal stop \
    --svg-strand +
```

#### SVG Plotting with Alignment
When sequence alignment is needed:

```bash
modtector plot \
    -o svg_output/ \
    --svg-template rna_structure.svg \
    --reactivity reactivity_data.csv \
    --svg-bases ATGC \
    --svg-signal all \
    --svg-strand + \
    --svg-ref reference.fa \
    --svg-max-shift 10
```

#### SVG Output Files
- **Multi-signal**: `rna_structure_colored_[signal_type].svg`
- **Single-signal**: `rna_structure_colored_score.svg`

### Plotting Options

#### Basic Plotting
```bash
modtector plot -M mod.csv -U unmod.csv -o plots/
```

#### With Reactivity Data
```bash
modtector plot -M mod.csv -U unmod.csv -o plots/ -r reactivity.csv
```

#### Custom Thresholds
```bash
modtector plot -M mod.csv -U unmod.csv -o plots/ -c 0.3 -d 100
```

#### With Genome Annotation
```bash
modtector plot -M mod.csv -U unmod.csv -o plots/ -g annotation.gff
```

## Advanced Features

### Multi-threading

modtector supports parallel processing for improved performance:

```bash
# Use 8 threads
modtector count -b sample.bam -f reference.fa -o output.csv -t 8
modtector reactivity -M mod.csv -U unmod.csv -O output.csv -t 8
modtector plot -M mod.csv -U unmod.csv -o plots/ -t 8
```

### Window Analysis

#### Fixed Windows
```bash
modtector count -b sample.bam -f reference.fa -o output.csv -w 1000
```

#### Dynamic Windows
```bash
modtector norm -i input.csv -o output.csv -m winsor90 --dynamic
```

### Base-Specific Analysis

Target specific bases for analysis:

```bash
# Analyze only A and C bases
modtector norm -i input.csv -o output.csv -m winsor90 --bases AC

# Analyze only G and T bases
modtector norm -i input.csv -o output.csv -m winsor90 --bases GT
```

### Statistical Testing

Compare samples using different statistical tests:

```bash
# Student's t-test
modtector compare -M mod.csv -U unmod.csv -o comparison.csv -t t-test

# Mann-Whitney U test
modtector compare -M mod.csv -U unmod.csv -o comparison.csv -t mann-whitney

# Wilcoxon signed-rank test
modtector compare -M mod.csv -U unmod.csv -o comparison.csv -t wilcoxon
```

### Auto-shift Correction

Automatically correct for sequence length differences:

```bash
modtector evaluate -r reactivity.csv -s structure.dp -o results/ -g gene_id --auto-shift
```

### Base Matching

Use intelligent base matching for T/U equivalence:

```bash
modtector evaluate -r reactivity.csv -s structure.dp -o results/ -g gene_id --base-matching
```

## Best Practices

### Data Quality

1. **High-quality alignments**: Ensure BAM files are properly aligned
2. **Sufficient coverage**: Use at least 50x coverage
3. **Proper controls**: Include unmodified control samples
4. **Quality filtering**: Remove low-quality reads and positions

### Parameter Selection

1. **Normalization method**: Choose based on data characteristics
2. **Window size**: Balance between noise reduction and signal preservation
3. **Thresholds**: Adjust based on expected signal levels
4. **Statistical tests**: Select appropriate test for your data

### Performance Optimization

1. **Thread count**: Match to available CPU cores
2. **Memory usage**: Monitor RAM usage for large datasets
3. **Disk space**: Ensure sufficient storage for output files
4. **Batch processing**: Process multiple samples in parallel

### Result Interpretation

1. **Check evaluation metrics**: Ensure good performance scores
2. **Validate with known sites**: Compare with literature
3. **Consider biological context**: Interpret results in context
4. **Reproducibility**: Document parameters and methods

### Troubleshooting

1. **Low coverage**: Increase sequencing depth
2. **Poor normalization**: Try different methods
3. **Low evaluation scores**: Check data quality
4. **Memory issues**: Reduce thread count or dataset size

### Documentation

1. **Record parameters**: Keep track of all settings
2. **Version control**: Use version control for reproducibility
3. **Log files**: Save and review log files
4. **Backup results**: Keep copies of important results

## Common Workflows

### Basic Workflow
```bash
# 1. Generate pileup data
modtector count -b mod.bam -f reference.fa -o mod_count.csv -t 4
modtector count -b unmod.bam -f reference.fa -o unmod_count.csv -t 4

# 2. Normalize signals
modtector norm -i mod_count.csv -o mod_norm.csv -m winsor90 --bases AC
modtector norm -i unmod_count.csv -o unmod_norm.csv -m winsor90 --bases AC

# 3. Calculate reactivity
modtector reactivity -M mod_norm.csv -U unmod_norm.csv -O reactivity.csv -t 4

# 4. Generate plots
modtector plot -M mod_norm.csv -U unmod_norm.csv -o plots/ -r reactivity.csv -t 4

# 5. Evaluate accuracy
modtector evaluate -r reactivity.csv -s structure.dp -o results/ -g gene_id
```

### Advanced Workflow
```bash
# 1. Generate pileup data with windowing
modtector count -b mod.bam -f reference.fa -o mod_count.csv -w 1000 -t 8
modtector count -b unmod.bam -f reference.fa -o unmod_count.csv -w 1000 -t 8

# 2. Normalize with dynamic windows
modtector norm -i mod_count.csv -o mod_norm.csv -m winsor90 --dynamic --bases AC
modtector norm -i unmod_count.csv -o unmod_norm.csv -m winsor90 --dynamic --bases AC

# 3. Calculate reactivity with custom parameters
modtector reactivity -M mod_norm.csv -U unmod_norm.csv -O reactivity.csv \
    -s ding -m siegfried --pseudocount 0.5 --maxscore 5 -t 8

# 4. Compare samples
modtector compare -M mod_norm.csv -U unmod_norm.csv -o comparison.csv \
    -t mann-whitney -d 20 -f 1.5

# 5. Generate comprehensive plots
modtector plot -M mod_norm.csv -U unmod_norm.csv -o plots/ \
    -r reactivity.csv -c 0.3 -d 100 -t 8

# 6. Evaluate with auto-shift correction
modtector evaluate -r reactivity.csv -s structure.dp -o results/ \
    -g gene_id --auto-shift --base-matching
```

This user guide provides comprehensive information for using modtector effectively. For specific command details, refer to the [Command Reference](commands.md).
