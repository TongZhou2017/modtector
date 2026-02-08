# Command Reference

This page provides comprehensive documentation for all modtector commands and their options.

## Overview

modtector provides ten main commands:

- `count` - Generate pileup data from BAM files
- `norm` - Normalize and filter signals
- `reactivity` - Calculate reactivity scores
- `compare` - Compare samples and identify differences
- `evaluate` - Evaluate accuracy using secondary structure
- `plot` - Generate visualization plots
- `duet` - Decompose normalized reactivity into dynamic ensembles using read-level co-variation
- `extract` - Extract gene regions from count results using GTF annotation
- `convert` - Convert various input formats to modtector pileup CSV format
- `correct` - Apply PCR bias correction to pileup CSV file

## count - Data Processing

Generate pileup data from BAM files by counting stop and mutation signals.

### Usage
```bash
modtector count [OPTIONS] --bam <BAM> --fasta <FASTA> --output <OUTPUT>
```

### Options

| Option | Short | Type | Default | Description |
|--------|-------|------|---------|-------------|
| `--bam` | `-b` | String | Required | BAM file path (or glob pattern for batch/single-cell mode) |
| `--fasta` | `-f` | String | Required | Reference FASTA file path |
| `--output` | `-o` | String | Required | Output CSV path (or directory for batch/single-cell mode) |
| `--strand` | `-s` | String | "+/-" | Strand filter (`+`, `-`, or `+/-`) |
| `--threads` | `-t` | Integer | 1 | Number of parallel threads |
| `--log` | `-l` | String | None | Log file path (optional) |
| `--window` | `-w` | Integer | None | Window size for genome segmentation (bases) |
| `--batch` | | Flag | false | Enable batch mode: process multiple BAM files sequentially |
| `--single-cell` | | Flag | false | Enable single-cell unified mode: unified processing with cell labels |
| `--min-base-qual` | | Integer | 0 | Minimum base quality score (Phred score) to count a mutation. Only mutations with quality >= threshold are counted. Recommended: 20 (same as RNAFramework rf-count). Set to 0 to disable. |
| `--pcr-bias-correction` | | Flag | false | Enable PCR bias correction for depth calculation. When enabled, effective_depth will be calculated and used instead of raw depth. |
| `--weight-increase` | | Float | 1.0 | Weight for increasing depth correction (used when depth is too low) |
| `--weight-decrease` | | Float | 0.5 | Weight for decreasing depth correction (used when depth is too high) |

### Examples

```bash
# Basic usage (standard mode)
modtector count -b sample.bam -f reference.fa -o output.csv

# With multiple threads
modtector count -b sample.bam -f reference.fa -o output.csv -t 8

# With windowing
modtector count -b sample.bam -f reference.fa -o output.csv -w 1000

# With logging
modtector count -b sample.bam -f reference.fa -o output.csv -l count.log

# Batch mode: process multiple BAM files sequentially
modtector count --batch \
  -b "/path/to/bam/*sort.bam" \
  -f reference.fa \
  -o output_dir/ \
  -t 8 \
  -w 10000

# Single-cell unified mode: unified processing with cell labels
modtector count --single-cell \
  -b "/path/to/single_cell/*sort.bam" \
  -f reference.fa \
  -o output_dir/ \
  -t 8 \
  -w 10000 \
  -l batch.log
```

### Batch and Single-cell Modes

**Batch Mode (`--batch`)**:
- Processes multiple BAM files sequentially
- Each file is processed independently
- Suitable for scenarios requiring separate processing per file
- Output: Each file generates a separate CSV file in the output directory
- Example glob patterns: `*sort.bam`, `*_RHX*.bam`, `sample_*.bam`

**Single-cell Unified Mode (`--single-cell`)**:
- True unified processing strategy for single-cell data
- **Processing Flow**:
  1. For each window, collects all reads from all BAM files at once
  2. Performs unified pileup processing (all reads together)
  3. Tracks cell labels during processing, splits results by cell
  4. Outputs separate CSV file for each cell
- **Performance Benefits** (2-3x speedup):
  - Skips data distribution scanning (direct processing of all reference sequences)
  - Unified read collection and pileup (reduces I/O overhead significantly)
  - Better parallelization efficiency (cross-file parallel processing)
  - Window-based memory management (prevents memory accumulation)
- **Progress Reporting**:
  - Real-time display: chunks processed, percentage, speed (chunks/s), ETA
  - CPU usage recommendations based on BAM file count
- **Output**: Each cell generates a separate CSV file (cell label extracted from filename)
- **Cell Label Extraction**: Supports RHX pattern (e.g., `*RHX672.sort.bam` → `RHX672`) or last underscore-separated part

### Output Format

The output CSV contains the following columns:
- `ChrID`: Chromosome/contig identifier
- `pipe_truncation_Strand`: Strand (+ or -)
- `pipe_truncation_ChrPos`: Position (1-based)
- `rf_mutation_Base`: Reference base
- `rf_mutation_Count`: Mutation count
- `pipe_truncation_count`: Stop count
- `depth`: Total read depth
- `rf_mutation_ins`: Insertion count
- `rf_mutation_del`: Deletion count
- `base_A`, `base_C`, `base_G`, `base_T`: Base counts

## norm - Signal Normalization

Normalize and filter signals to remove noise and outliers.

### Usage
```bash
modtector norm [OPTIONS] --input <INPUT> --output <OUTPUT> --method <METHOD>
```

### Options

| Option | Short | Type | Default | Description |
|--------|-------|------|---------|-------------|
| `--input` | `-i` | String | Required | Input CSV file |
| `--output` | `-o` | String | Required | Output CSV file |
| `--method` | `-m` | String | Required | Normalization method |
| `--cutoff` | | Float | 0.25 | SNP threshold for filtering mutation signals |
| `--bases` | | String | "ACGT" | Reactive base types (e.g., "AC" for A and C) |
| `--window` | | Integer | 0 | Fixed window size (fixed sliding window mode) |
| `--window-offset` | | Integer | 0 | Window offset (fixed sliding window mode) |
| `--dynamic` | | Flag | False | Enable dynamic sliding window mode |
| `--window-size-dynamic` | | Integer | 50 | Number of reactive bases required per window |
| `--linear` | | Flag | False | Apply piecewise linear mapping |
| `--log` | `-l` | String | None | Log file path (optional) |

### Normalization Methods

1. **percentile28**: 28th percentile normalization
2. **winsor90**: 90th percentile winsorization
3. **boxplot**: Boxplot-based outlier removal

### Examples

```bash
# Basic normalization
modtector norm -i input.csv -o output.csv -m winsor90

# Target specific bases
modtector norm -i input.csv -o output.csv -m winsor90 --bases AC

# With dynamic windowing
modtector norm -i input.csv -o output.csv -m winsor90 --dynamic

# With linear mapping
modtector norm -i input.csv -o output.csv -m winsor90 --linear
```

## reactivity - Reactivity Calculation

Calculate reactivity scores by comparing modified and unmodified samples.

### Usage
```bash
modtector reactivity [OPTIONS] --mod <MOD_CSV> --unmod <UNMOD_CSV> --output <OUTPUT>
```

### Options

| Option | Short | Type | Default | Description |
|--------|-------|------|---------|-------------|
| `--mod` | `-M` | String | Required | Modified sample CSV |
| `--unmod` | `-U` | String | Optional | Unmodified sample CSV (optional for mod-only mode) |
| `--output` | `-O` | String | Required | Reactivity output file |
| `--stop-method` | `-s` | String | "kfactor" | Stop signal reactivity method: `kfactor`, `ding`, `rouskin` |
| `--mutation-method` | `-m` | String | "kfactor" | Mutation signal reactivity method: `kfactor`, `siegfried`, `zubradt` |
| `--pseudocount` | | Float | 1.0 | Pseudocount parameter for logarithmic calculations (Ding method) |
| `--maxscore` | | Float | 10.0 | Maximum score limit parameter (Ding method) |
| `--snp-cutoff` | | Float | 0.25 | SNP threshold for filtering positions (filter positions where unmod sample mutation rate >= cutoff) |
| `--k-prediction-method` | | String | "background" | K-factor prediction method: `background` (default), `distribution`, `recursive` |
| `--structure-file` | | String | None | Reference secondary structure file (required for recursive k-factor method) |
| `--k-background-gene-id` | | String | None | Gene ID to use for k-factor calculation (only positions with this gene_id will be used as background regions) |
| `--threads` | `-t` | Integer | 1 | Number of parallel threads |
| `--log` | `-l` | String | None | Log file path (optional) |

### Reactivity Methods

#### Stop Signal Methods
- **kfactor** (default): K-factor method using difference between modified and unmodified samples
- **ding**: Ding et al. method with pseudocount and maxscore parameters
- **rouskin**: Rouskin et al. method (uses only modified sample stop signal)

#### Mutation Signal Methods
- **kfactor** (default): K-factor method using difference between modified and unmodified samples
- **siegfried**: Siegfried et al. method (allows negative reactivity values)
- **zubradt**: Zubradt et al. method

### K-Factor Prediction Methods

- **background** (default): Predict k-factor from background regions (low mutation rate positions)
- **distribution**: Predict k-factor using statistical distribution analysis (improved accuracy for many sample types)
- **recursive**: Predict k-factor recursively using reference secondary structure (requires `--structure-file`)

### Examples

```bash
# Basic reactivity calculation (with unmodified sample)
modtector reactivity -M mod.csv -U unmod.csv -O reactivity.csv

# Mod-only mode (no unmodified sample, for smartSHAPE datasets)
modtector reactivity -M mod.csv -O reactivity.csv

# With specific methods
modtector reactivity -M mod.csv -U unmod.csv -O reactivity.csv -s ding -m siegfried

# With custom parameters
modtector reactivity -M mod.csv -U unmod.csv -O reactivity.csv --pseudocount 0.5 --maxscore 5

# Using distribution-based k-factor prediction
modtector reactivity -M mod.csv -U unmod.csv -O reactivity.csv --k-prediction-method distribution

# Using recursive k-factor prediction (requires structure file)
modtector reactivity -M mod.csv -U unmod.csv -O reactivity.csv \
    --k-prediction-method recursive \
    --structure-file structure.dp \
    --k-background-gene-id gene_id

# With SNP filtering
modtector reactivity -M mod.csv -U unmod.csv -O reactivity.csv --snp-cutoff 0.2
```

## compare - Sample Comparison

Compare modified and unmodified samples to identify significant differences.

### Usage
```bash
modtector compare [OPTIONS] --output <OUTPUT>
```

### Options

| Option | Short | Type | Default | Description |
|--------|-------|------|---------|-------------|
| `--mode` | `-m` | String | "mod-vs-unmod" | Comparison mode |
| `--mod` | `-M` | String | None | Modified sample CSV (mod-vs-unmod mode) |
| `--unmod` | `-U` | String | None | Unmodified sample CSV (mod-vs-unmod mode) |
| `--group1` | | String | None | Group 1 file list (biological-replicates mode) |
| `--group2` | | String | None | Group 2 file list (biological-replicates mode) |
| `--reactivity-group1` | | String | None | Group 1 reactivity files (reactivity-groups mode) |
| `--reactivity-group2` | | String | None | Group 2 reactivity files (reactivity-groups mode) |
| `--reactivity1` | | String | None | First reactivity result file (reactivity-results mode) |
| `--reactivity2` | | String | None | Second reactivity result file (reactivity-results mode) |
| `--output` | `-o` | String | Required | Output CSV path |
| `--min-depth` | `-d` | Integer | 10 | Minimum depth |
| `--min-fold` | `-f` | Float | 2 | Minimum fold change |
| `--test` | `-t` | String | "t-test" | Statistical test type |
| `--alpha` | `-a` | Float | 0.05 | Significance level |
| `--log` | `-l` | String | None | Log file path (optional) |

### Comparison Modes

1. **mod-vs-unmod**: Compare modified vs unmodified samples
2. **reactivity-groups**: Compare two groups of reactivity results
3. **biological-replicates**: Compare biological replicates
4. **reactivity-results**: Compare two reactivity result files

### Statistical Tests

- **t-test**: Student's t-test
- **mann-whitney**: Mann-Whitney U test
- **wilcoxon**: Wilcoxon signed-rank test
- **chi-square**: Chi-square test
- **continuity**: Continuity correction
- **diffscan**: Differential scanning
- **deltashape**: DeltaSHAPE method

### Examples

```bash
# Basic comparison
modtector compare -M mod.csv -U unmod.csv -o comparison.csv

# With custom thresholds
modtector compare -M mod.csv -U unmod.csv -o comparison.csv -d 20 -f 1.5

# Using different statistical test
modtector compare -M mod.csv -U unmod.csv -o comparison.csv -t mann-whitney

# Compare reactivity groups
modtector compare --mode reactivity-groups \
    --reactivity-group1 group1_reactivity.csv \
    --reactivity-group2 group2_reactivity.csv \
    -o comparison.csv
```

## evaluate - Accuracy Evaluation

Evaluate reactivity accuracy using known secondary structure.

### Usage
```bash
modtector evaluate [OPTIONS] --reactivity <REACTIVITY_FILE> --structure <STRUCTURE_FILE> --output <OUTPUT_DIR> --gene-id <GENE_ID>
```

### Options

| Option | Short | Type | Default | Description |
|--------|-------|------|---------|-------------|
| `--reactivity` | `-r` | String | Required | Reactivity signal file |
| `--structure` | `-s` | String | Required | Secondary structure file |
| `--output` | `-o` | String | Required | Output directory |
| `--log` | `-l` | String | None | Log file path (default to result directory) |
| `--signal-type` | `-t` | String | "stop" | Signal type (`stop`, `mutation`, or `both`) |
| `--gene-id` | `-g` | String | Required | Gene ID (for base matching) |
| `--strand` | `-S` | String | "+" | Strand information (`+` or `-`) |
| `--base-matching` | | Flag | True | Use base matching (default: true) |
| `--auto-shift` | | Flag | True | Use auto-shift correction (default: true) |
| `--reactive-bases` | | String | "ACGT" | Reactive bases for AUC/ROC calculation (e.g., `AC` for DMS; default `ACGT` for all bases) |
| `--optimized` | `-O` | Flag | False | Use optimized version |

### Examples

```bash
# Basic evaluation
modtector evaluate -r reactivity.csv -s structure.dp -o results/ -g gene_id

# Evaluate both signal types
modtector evaluate -r reactivity.csv -s structure.dp -o results/ -g gene_id -t both

# Without auto-shift correction
modtector evaluate -r reactivity.csv -s structure.dp -o results/ -g gene_id --no-auto-shift

# Using optimized version
modtector evaluate -r reactivity.csv -s structure.dp -o results/ -g gene_id -O
```

## plot - Visualization

Generate signal distribution plots and reactivity visualizations, including RNA structure SVG plots.

### Usage

#### Regular Mode
```bash
modtector plot [OPTIONS] --mod <MOD_CSV> --unmod <UNMOD_CSV> --output <OUTPUT>
```

#### SVG-Only Mode
```bash
modtector plot [OPTIONS] --output <OUTPUT> --svg-template <SVG_TEMPLATE> --reactivity <REACTIVITY_CSV>
```

### Options

| Option | Short | Type | Default | Description |
|--------|-------|------|---------|-------------|
| `--mod` | `-M` | String | Optional* | Modified sample CSV (required for regular mode) |
| `--unmod` | `-U` | String | Optional* | Unmodified sample CSV (required for regular mode) |
| `--output` | `-o` | String | Required | Output directory |
| `--coverage` | `-c` | Float | 0.2 | Coverage threshold |
| `--depth` | `-d` | Integer | 50 | Depth threshold (reads) |
| `--reactivity` | `-r` | String | None | Reactivity CSV file (optional) |
| `--gff` | `-g` | String | None | Genome annotation GFF/GTF file (optional) |
| `--threads` | `-t` | Integer | 1 | Number of parallel threads |
| `--log` | `-l` | String | None | Log file path (optional) |

### SVG Plotting Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--svg-template` | String | None | SVG template file for RNA structure visualization |
| `--svg-bases` | String | "ACGT" | Bases to include in SVG plot (e.g., "AC" for DMS-seq) |
| `--svg-signal` | String | "all" | Signal type to plot (`all`, `stop`, `mutation`, `score`) |
| `--svg-strand` | String | "+" | Strand to include (`+`, `-`, `both`) |
| `--svg-ref` | String | None | Reference sequence file for alignment (optional) |
| `--svg-max-shift` | Integer | 5 | Maximum shift for alignment search |
| `--svg-circle-filled` | Flag | False | Circle fill type: `true` for filled circles, `false` for hollow rings |
| `--svg-font-color` | String | None | Font color for text elements (e.g., "black", "red", "#000000") |
| `--svg-legend-width` | Float | 30.0 | Legend item width for horizontal layout |
| `--svg-legend-height` | Float | 15.0 | Legend item height for vertical layout |
| `--interactive` | Flag | False | Generate interactive HTML visualization instead of static SVG |

*Note: `-M` and `-U` parameters are optional when using SVG-only mode (when both `--svg-template` and `--reactivity` are provided).*

### Examples

#### Regular Plotting
```bash
# Basic plotting
modtector plot -M mod.csv -U unmod.csv -o plots/

# With reactivity data
modtector plot -M mod.csv -U unmod.csv -o plots/ -r reactivity.csv

# With custom thresholds
modtector plot -M mod.csv -U unmod.csv -o plots/ -c 0.3 -d 100

# With genome annotation
modtector plot -M mod.csv -U unmod.csv -o plots/ -g annotation.gff
```

#### SVG-Only Mode
```bash
# Basic SVG plotting (simplified command)
modtector plot -o output/ --svg-template template.svg --reactivity data.csv

# Multi-signal SVG plotting
modtector plot -o output/ \
    --svg-template template.svg \
    --reactivity data.csv \
    --svg-bases ATGC \
    --svg-signal all \
    --svg-strand +

# Single signal SVG plotting
modtector plot -o output/ \
    --svg-template template.svg \
    --reactivity data.csv \
    --svg-bases AC \
    --svg-signal stop \
    --svg-strand +

# SVG plotting with alignment
modtector plot -o output/ \
    --svg-template template.svg \
    --reactivity data.csv \
    --svg-bases ATGC \
    --svg-signal all \
    --svg-strand + \
    --svg-ref reference.fa \
    --svg-max-shift 10

# SVG plotting with custom styling
modtector plot -o output/ \
    --svg-template template.svg \
    --reactivity data.csv \
    --svg-circle-filled \
    --svg-font-color "black" \
    --svg-legend-width 40.0 \
    --svg-legend-height 20.0

# Interactive HTML visualization
modtector plot -o output/ \
    --svg-template template.svg \
    --reactivity data.csv \
    --interactive
```

#### Combined Mode (Regular + SVG)
```bash
# Generate both regular plots and SVG plots
modtector plot -M mod.csv -U unmod.csv -o output/ \
    --svg-template template.svg \
    --reactivity data.csv \
    --svg-bases ATGC \
    --svg-signal all
```

### SVG Plotting Features

#### Signal Type Support
- **Auto-detection**: Automatically detects signal types from CSV headers
- **Multiple signals**: Supports plotting multiple signal types (stop, mutation, etc.)
- **Single signal**: Works with single signal files
- **Custom naming**: Output files named by signal type (e.g., `rna_structure_colored_stop.svg`)

#### Strand Selection
- **Positive strand**: `--svg-strand +` (default)
- **Negative strand**: `--svg-strand -`
- **Both strands**: `--svg-strand both`

#### Base Filtering
- **All bases**: `--svg-bases ATGC` (default)
- **DMS-seq**: `--svg-bases AC` (A and C only)
- **Custom**: Any combination of A, T, C, G

#### Alignment Support
- **Reference sequence**: Use `--svg-ref` for sequence alignment
- **Shift calculation**: Automatic calculation of optimal shift values
- **Max shift**: Control maximum shift search range with `--svg-max-shift`

### Output Files

#### Regular Mode Output
- `overall_scatter.png`: Signal distribution plot
- `high/`: High coverage gene plots
- `low/`: Low coverage gene plots

#### SVG Mode Output
- `rna_structure_colored_[signal].svg`: SVG files for each signal type
- Multiple SVG files for multi-signal data
- Single SVG file for single-signal data

#### Combined Mode Output
- All regular mode files
- All SVG mode files

## duet - Dynamic Ensemble Decomposition

Duet performs a sliding-window analysis that merges normalized stop/mutation reactivity with read-level co-variation to detect alternative ensembles without predefining cluster counts. Within each window, read-level stop/mutation co-occurrence vectors are clustered with DBSCAN to highlight dynamic structural mixtures.

### Usage
```bash
modtector duet [OPTIONS] --input <INPUT> --bam <BAM> --fasta <FASTA> --output <OUTPUT>
```

### Options

| Option | Short | Type | Default | Description |
|--------|-------|------|---------|-------------|
| `--input` | `-i` | String | Required | Normalized reactivity CSV file |
| `--bam` | `-b` | String | Required | Sorted BAM file for co-variation analysis |
| `--fasta` | `-f` | String | Required | Reference FASTA used for alignment |
| `--output` | `-o` | String | Required | Window-level Duet results CSV |
| `--epsilon` |  | Float | 0.75 | DBSCAN ε (radius) in standardized feature space |
| `--min-samples` |  | Integer | 5 | Minimum neighbours required to form a DBSCAN core read |
| `--window-size` |  | Integer | 100 | Sliding-window size in nucleotides |
| `--window-step` |  | Integer | 50 | Sliding-window step in nucleotides |
| `--threads` | `-t` | Integer | System CPU count | Number of threads for parallel processing |
| `--summary-output` |  | String | `<output>_summary.csv` | Optional per-window/per-ensemble summary CSV |
| `--log` | `-l` | String | None | Log file path (optional) |

### Outputs
- **Window CSV**: Each window annotated with status (`OK`, `InsufficientReads`, or `LowDensity`), read counts, detected ensemble count, primary ensemble occupancy, noise fraction, dual-signal fraction, window-averaged reactivities, a concise cluster summary string (including per-cluster confidence), and a `GlobalEnsembles` mapping column.
- **Summary CSV**: Per window/per ensemble (including DBSCAN noise) details covering occupancy, cluster confidence, average stop/mutation counts per read, dual-signal fractions, and the associated global ensemble ID. Generated automatically unless `--summary-output` is provided.
- **Global Ensemble CSV**: (`<output>_global.csv`) summarises each global ensemble with aggregated read counts, stop/mutation totals, dual-signal read counts, contributing window/cluster counts, unique position counts, and a global confidence score.
- **Global Per-base CSV**: (`<output>_global_per_base.csv`) lists every position associated with each global ensemble alongside ensemble-specific read support, normalized reactivity values, and per-base confidence scores.

### Progress Reporting

During execution, duet provides detailed statistics:
- **Dual-signal reads**: Number and percentage of reads carrying both stop and mutation signals
- **High-confidence positions**: Number and percentage of positions assigned high confidence (reactivity >= 0.7)
- **Average confidence**: Mean position confidence score across all analyzed positions

These statistics help assess data quality and the reliability of ensemble detection.

### Examples

```bash
# Sliding-window Duet analysis (100-nt windows, 50-nt step, 16 threads)
modtector duet \
  -i signal/03_norm/iso-1_rep-1_norm.csv \
  -b signal/00_bam/iso-1_rep-1.sort.bam \
  -f reference/PDL1.fa \
  -o signal/04_duet/iso-1_rep-1_windows.csv \
  -t 16

# Custom density and window configuration with explicit summary output
modtector duet \
  -i signal/03_norm/sample_norm.csv \
  -b signal/00_bam/sample.sort.bam \
  -f reference/transcriptome.fa \
  -o signal/04_duet/sample_windows.csv \
  --epsilon 0.9 \
  --min-samples 8 \
  --window-size 120 \
  --window-step 30 \
  --summary-output signal/04_duet/sample_windows_summary.csv
```

## extract - Extract Gene Regions from Count Results

Extract count results for specific gene regions using GTF/GFF annotation files. This command allows you to filter count CSV files to include only positions within specified gene regions, such as 18S rRNA or other target genes.

### Usage
```bash
modtector extract [OPTIONS] --input <INPUT> --gtf <GTF> --target-gene <TARGET_GENE> --output <OUTPUT>
```

### Options

| Option | Short | Type | Default | Description |
|--------|-------|------|---------|-------------|
| `--input` | `-i` | String | Required | Count result CSV file from `count` command, or bam-readcount output file |
| `--gtf` | `-g` | String | Required | GTF/GFF annotation file |
| `--target-gene` | `-t` | String | Required | Target gene name or ID (e.g., "18S", "RN18S1", "rRNA") |
| `--output` | `-o` | String | Required | Output file prefix (each gene region saved as: `{prefix}_{gene_name}_depth{avg_depth}.csv`) |
| `--input-format` | | String | "csv" | Input format: "csv" (default) or "bam-readcount" |
| `--threads` | | Integer | 1 | Number of parallel threads for processing multiple gene regions |
| `--relative-position` | | Flag | false | Use relative position (1-based, relative to gene start) instead of absolute genomic position |
| `--log` | `-l` | String | None | Log file path (optional) |

### How It Works

1. **Parse GTF file**: Reads the GTF/GFF annotation file and identifies gene regions matching the target gene name or ID
2. **Match gene regions**: Searches for genes where `gene_id` or `gene_name` contains the target string
3. **Extract positions**: For each matched gene region, filters the count CSV to include only positions within that region
4. **Calculate statistics**: Computes average depth for each gene region
5. **Output separate files**: Creates a separate CSV file for each gene region with filename format: `{prefix}_{gene_name}_depth{avg_depth}.csv`

### Gene Matching

The command matches genes by:
- **gene_id**: Exact match or substring match
- **gene_name**: Exact match or substring match

For example, searching for "18S" will match:
- `RN18S1` (gene_id)
- `18S_rRNA` (gene_name)
- Any gene containing "18S" in its ID or name

### Examples

```bash
# Extract 18S rRNA gene regions from count results (CSV format)
# If multiple regions match, each will be saved as a separate file
modtector extract \
  -i benchmark/sc/count/single_cell_transcripts_in_hESC_NAIN3_RHE1581_b1.sort.csv \
  -g annotation.gtf \
  -t "18S" \
  -o 18S_rRNA \
  -l extract.log

# This will create files like:
# - 18S_rRNA_RN18S1_depth1234.5.csv
# - 18S_rRNA_RN18S2_depth987.2.csv
# (if multiple 18S genes are found)

# Extract from bam-readcount output format with multithreading
modtector extract \
  -i /data1/bioinfo/07_people/lxt/Modector/PRJNA946273_fixed_bamreadcount/in_vivo_single_cell_Mut_transfected_RNAs_in_HEK293T_DMSO_RHX672.sort.readcount.txt \
  -g annotation.gtf \
  -t "18S" \
  --input-format bam-readcount \
  --threads 8 \
  -o 18S_rRNA \
  -l extract.log

# Extract with multithreading (CSV format)
modtector extract \
  -i count_results.csv \
  -g annotation.gtf \
  -t "18S" \
  --threads 8 \
  -o 18S_output

# Extract with relative position (1-based, relative to gene start)
# Position will be converted from absolute genomic position to relative position
# Example: If gene starts at position 125931, absolute position 125931 becomes relative position 1
modtector extract \
  -i count_results.csv \
  -g annotation.gtf \
  -t "18S" \
  --relative-position \
  -o 18S_rRNA_relpos

# Extract specific gene by gene_id
modtector extract \
  -i count_results.csv \
  -g annotation.gtf \
  -t "RN18S1" \
  -o RN18S1 \
  -l extract.log

# Extract gene by partial name match (may match multiple genes)
modtector extract \
  -i count_results.csv \
  -g annotation.gtf \
  -t "rRNA" \
  -o rRNA_genes
```

### Output Format

Each output CSV file maintains the same format as the input count CSV:
```
ChrID,pipe_truncation_Strand,pipe_truncation_ChrPos,rf_mutation_Base,rf_mutation_Count,pipe_truncation_count,depth,rf_mutation_ins,rf_mutation_del,base_A,base_C,base_G,base_T
```

Only rows matching the specific gene region are included in each output file.

**Position Format**:
- **Default (absolute position)**: `pipe_truncation_ChrPos` contains the absolute genomic position (e.g., 125931)
- **With `--relative-position`**: `pipe_truncation_ChrPos` contains the relative position (1-based, relative to gene start). For example, if a gene starts at position 125931, the first position in the gene will be 1, the second will be 2, etc.

**Example**:
- Gene region: NT_187388.1:125931-127799
- Absolute position 125931 → Relative position 1 (with `--relative-position`)
- Absolute position 125932 → Relative position 2 (with `--relative-position`)

### File Naming

Output files are named using the pattern: `{prefix}_{gene_name}_depth{avg_depth}.csv`

- `{prefix}`: The output prefix specified with `-o`
- `{gene_name}`: The gene name from GTF (or gene_id if gene_name is not available)
- `{avg_depth}`: The average depth calculated for that gene region (rounded to 1 decimal place)

Example filenames:
- `18S_rRNA_RN18S1_depth1234.5.csv`
- `rRNA_genes_28S_rRNA_depth987.2.csv`

### Input Formats

#### CSV Format (default)
Standard CSV format from `modtector count` command with the following columns:
- `ChrID`, `pipe_truncation_Strand`, `pipe_truncation_ChrPos`, `rf_mutation_Base`, `rf_mutation_Count`, `pipe_truncation_count`, `depth`, `rf_mutation_ins`, `rf_mutation_del`, `base_A`, `base_C`, `base_G`, `base_T`

#### bam-readcount Format
Tab-separated format from bam-readcount software:
- Format: `chr position reference_base depth base:count:avg_mapping_quality:...`
- Each base (A, C, G, T, N, =, +, -) has colon-separated statistics
- **Note**: bam-readcount format does not include stop signals (`pipe_truncation_count`), so this value will be set to 0 in the output
- Strand information is taken from the GTF file (not available in bam-readcount output)

### Performance Optimization

- **Multithreading**: Use `--threads` to enable parallel processing of multiple gene regions
  - Recommended: 4-8 threads for optimal performance
  - Automatically falls back to sequential processing for single region or single thread
  - Provides 6-40× speedup depending on format and thread count
- **bam-readcount Format**: File is read once into memory, then processed in parallel
  - Significantly faster than sequential processing (40× speedup with 8 threads)
  - Memory usage: Entire file loaded into memory (consider file size)

### Notes

- The GTF file should contain `gene_id` and/or `gene_name` attributes
- When using bam-readcount format, stop signals are not available and will be set to 0
- Both GTF and GFF3 formats are supported
- Gene coordinates are 1-based (GTF standard)
- **Each matching gene region is saved to a separate file**
- If multiple genes match, each will have its own output file with its gene name and average depth in the filename
- Gene names in filenames are sanitized (special characters like `/`, `:`, spaces are replaced with `_`)
- Files with zero extracted rows are not created
- **Multithreading**: Parallel processing is most effective when extracting multiple gene regions

## convert - Format Conversion

Convert various input formats to modtector pileup CSV format. Supports multiple input formats including bamreadcount, rf-rctools, rf-norm, shapemapper2, samtools-mpileup, icSHAPE-rt, and bedgraph formats.

### Usage
```bash
modtector convert [OPTIONS] --input <INPUT> --output <OUTPUT>
```

### Options

| Option | Short | Type | Default | Description |
|--------|-------|------|---------|-------------|
| `--input` | `-i` | String | Required* | Input file (or use `--input-mutation`/`--input-stop` for dual input mode) |
| `--output` | `-o` | String | Required | Output pileup CSV file |
| `--format` | `-f` | String | Auto-detect | Input format: `bamreadcount`, `rf-rctools`, `rf-norm`, `rf-norm-xml`, `shapemapper2`, `shapemapper-profile`, `samtools-mpileup`, `icSHAPE-rt`, or `bedgraph` |
| `--strand` | `-s` | String | "+" | Strand orientation ('+' for forward, '-' for reverse, or 'None' for all strands) |
| `--rf-count-mutations` | | Flag | false | For rf-rctools format: if true, input was generated with `rf-count -m` (mutation mode) |
| `--ref-fasta` | | String | None | Reference FASTA file path (required for shapemapper2, icSHAPE-rt, and bedgraph formats) |
| `--chr-name` | | String | None | Chromosome/contig name (for formats with missing metadata) |
| `--input-mutation` | | String | None | Optional mutation input file (rf-rctools format with `-m` flag) |
| `--input-stop` | | String | None | Optional stop input file (rf-rctools format without `-m` flag) |
| `--filter-strand` | | String | "+" | Filter by strand orientation ('+' or '-' or 'None' to process all strands) |
| `--log` | `-l` | String | None | Log file path (optional) |

*Note: In dual input mode (`--input-mutation`/`--input-stop`), `--input` can be empty.

### Supported Input Formats

1. **bamreadcount**: Tab-separated format from bam-readcount software
2. **rf-rctools**: RNA Framework / rctools counting output (with optional `--rf-count-mutations` flag)
3. **rf-norm**: rf-norm CSV output (reactivity-level data)
4. **rf-norm-xml**: rf-norm XML output
5. **shapemapper2**: ShapeMapper2 profile files (requires `--ref-fasta`)
6. **shapemapper-profile**: ShapeMapper2 profile format with dynamic column detection
7. **samtools-mpileup**: Samtools mpileup output (10-column tab-separated format)
8. **icSHAPE-rt**: icSHAPE-pipe RT file format (requires `--ref-fasta`)
9. **bedgraph**: Bedgraph format (5-column: chr, start, end, coverage, strand; requires `--ref-fasta`)

### Examples

```bash
# Basic usage (bamreadcount format, auto-detected)
modtector convert -i bamreadcount.txt -o pileup.csv

# With strand specification
modtector convert -i bamreadcount.txt -o pileup.csv -s +

# Convert rf-rctools format (mutation mode)
modtector convert -i mutation.tsv -o pileup.csv -f rf-rctools --rf-count-mutations

# Convert rf-rctools format (stop mode)
modtector convert -i stop.tsv -o pileup.csv -f rf-rctools

# Dual input mode: merge mutation and stop files
modtector convert --input-mutation mutation.tsv --input-stop stop.tsv -o pileup.csv -f rf-rctools

# Convert shapemapper2 format (requires reference)
modtector convert -i shapemapper.profile -o pileup.csv -f shapemapper-profile --ref-fasta reference.fa

# Convert icSHAPE RT format
modtector convert -i input.rt -o pileup.csv -f icSHAPE-rt --ref-fasta reference.fa --filter-strand "+"

# Convert bedgraph format
modtector convert -i input.bedgraph -o pileup.csv -f bedgraph --ref-fasta reference.fa --filter-strand "+"

# Convert samtools mpileup format
modtector convert -i mpileup.txt -o pileup.csv -f samtools-mpileup -s +

# With logging
modtector convert -i input.txt -o pileup.csv -s + -l convert.log
```

### Input Format (bamreadcount)

Tab-separated format from bam-readcount software:
```
EC_16S	1	A	77	A:77:29.25:33.82:...	C:0:0.00:...	G:0:0.00:...	T:0:0.00:...
```

Fields:
- Column 1: Chromosome/contig name
- Column 2: Position (1-based)
- Column 3: Reference base
- Column 4: Total depth
- Column 5+: Base counts in format `base:count:quality:...`

### Output Format (modtector pileup)

CSV format compatible with modtector reactivity command:
```
ChrID,pipe_truncation_Strand,pipe_truncation_ChrPos,rf_mutation_Base,rf_mutation_Count,pipe_truncation_count,depth,rf_mutation_ins,rf_mutation_del,base_A,base_C,base_G,base_T
EC_16S,+,1,A,0,0,77,0,0,77,0,0,0
```

### Conversion Details

- **Mutation Count**: Calculated as sum of non-reference bases (A, C, G, T)
- **Stop Signals**: `pipe_truncation_count` is set to 0 (not available in bamreadcount format)
- **Strand Information**: Uses the specified strand parameter (default: '+')
- **Insertions/Deletions**: Extracted from bamreadcount format (`+` for insertions, `-` for deletions)

### Performance

- **Streaming Processing**: Processes files line-by-line without loading entire file into memory
- **Large File Support**: Tested with files up to 35GB
- **Processing Speed**: ~1-2 million lines/second (depends on I/O)
- **Memory Usage**: O(1) constant regardless of file size
- **NFS Optimization**: Uses 8MB buffer for better network filesystem performance

### Progress Reporting

The command reports progress every 1 million lines:
- Total lines processed
- Valid lines converted
- Error count (invalid lines skipped)

### Error Handling

- Invalid lines are skipped and error count is tracked
- First 10 errors are logged for debugging
- Processing continues despite errors
- Final statistics show total, valid, and error counts

### Use Cases

1. **Workflow Integration**: Convert bamreadcount output for use with modtector reactivity command
2. **Format Standardization**: Standardize bamreadcount data to modtector pileup format
3. **Pipeline Compatibility**: Enable bamreadcount-based workflows to use modtector reactivity features

### Notes

- Output file must have `.csv` extension
- Strand parameter must be either '+' or '-' (or 'None' for all strands)
- Input file must exist and be readable
- The conversion preserves all position and depth information
- Base counts (A, C, G, T) are extracted from input format
- Insertions and deletions are included if present in input format
- Format auto-detection works by checking file extension, headers, and data patterns
- NaN values are preserved in output for formats that support them (rf-norm)

## correct - PCR Bias Correction

Apply PCR bias correction to pileup CSV files using Chi-Square distribution fitting. This command corrects for PCR amplification bias where high depths cause mutation rate dilution, similar to shapemapper2's effective_depth calculation.

### Usage
```bash
modtector correct [OPTIONS] --input <INPUT> --output <OUTPUT>
```

### Options

| Option | Short | Type | Default | Description |
|--------|-------|------|---------|-------------|
| `--input` | `-i` | String | Required | Input pileup CSV file |
| `--output` | `-o` | String | Required | Output pileup CSV file (with corrected effective_depth) |
| `--weight-increase` | | Float | 1.0 | Weight for increasing depth correction (used when depth is too low) |
| `--weight-decrease` | | Float | 0.5 | Weight for decreasing depth correction (used when depth is too high) |
| `--species` | | String | None | Species name for filtering (optional) |

### How It Works

1. **Data Binning**: Groups depth and mutation rate data into bins
2. **Chi-Square Fitting**: Fits a Chi-Square distribution (df=2) to the depth-mutation rate relationship
3. **Correction Factor Calculation**: Computes correction factors based on the fitted distribution
4. **Effective Depth Calculation**: Applies correction factors to adjust `effective_depth` values
5. **Output**: Writes corrected pileup CSV with updated `effective_depth` column

### Method Details

The correction uses a Chi-Square distribution model:
- **Distribution**: Chi-Square with 2 degrees of freedom
- **Fitting**: Brent's method for 1D minimization to find optimal scale parameter
- **Correction**: Applies weights for increasing/decreasing depth based on fitted model

### Examples

```bash
# Basic usage
modtector correct -i pileup.csv -o pileup_corrected.csv

# With custom weights
modtector correct -i pileup.csv -o pileup_corrected.csv \
    --weight-increase 1.2 \
    --weight-decrease 0.4

# With species filtering
modtector correct -i pileup.csv -o pileup_corrected.csv \
    --species "human"
```

### Output Format

The output CSV maintains the same format as input, with the `effective_depth` column updated:
- Original `depth` column: Unchanged (total read depth)
- `effective_depth` column: Corrected depth based on PCR bias model

### Use Cases

1. **Quality Filtering Enhancement**: Combine with `--min-base-qual` in `count` command for comprehensive depth correction
2. **High-Depth Data**: Correct for mutation rate dilution in high-coverage regions
3. **Workflow Integration**: Optional step in analysis pipeline (can be disabled in default workflow)

### Notes

- This is an optional advanced feature
- Can be combined with quality filtering (`--min-base-qual`) for best results
- Correction factors are calculated per-species if `--species` is specified
- The method is based on Chi-Square distribution fitting, similar to shapemapper2's approach

## Global Options

All commands support these global options:

| Option | Short | Description |
|--------|-------|-------------|
| `--help` | `-h` | Print help information |
| `--version` | `-V` | Print version information |

## Output Files

### Common Output Formats

1. **CSV Files**: Comma-separated values with headers
2. **SVG Files**: Scalable vector graphics for plots
3. **Text Files**: Plain text for evaluation results
4. **Log Files**: Processing logs and error messages

### File Naming Conventions

- Input files: `*.bam`, `*.fa`, `*.csv`
- Output files: `*_count.csv`, `*_norm.csv`, `*_reactivity.csv`
- Plots: `*_distribution.svg`, `*_reactivity.svg`, `*_roc.svg`
- Evaluation: `*_evaluation.txt`, `*_comprehensive.txt`

## Best Practices

1. **Use appropriate thread counts**: Match `-t` to your CPU cores
2. **Monitor memory usage**: Large datasets may require significant RAM
3. **Check log files**: Always review logs for warnings and errors
4. **Validate inputs**: Ensure BAM files are properly aligned
5. **Use consistent parameters**: Keep parameters consistent across runs
6. **Backup results**: Save intermediate results for reproducibility

## Troubleshooting

### Common Issues

1. **Memory errors**: Reduce thread count or use smaller datasets
2. **File not found**: Check file paths and permissions
3. **Low coverage**: Increase sequencing depth or adjust thresholds
4. **Poor results**: Check data quality and parameter settings

### Getting Help

- Use `--help` for command-specific help
- Check log files for detailed error messages
- Refer to the [Troubleshooting Guide](troubleshooting.md)
- Open an issue on the [GitHub repository](https://github.com/TongZhou2017/modtector/issues)
