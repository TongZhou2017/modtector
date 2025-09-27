# Command Reference

This page provides comprehensive documentation for all modtector commands and their options.

## Overview

modtector provides six main commands:

- `count` - Generate pileup data from BAM files
- `norm` - Normalize and filter signals
- `reactivity` - Calculate reactivity scores
- `compare` - Compare samples and identify differences
- `evaluate` - Evaluate accuracy using secondary structure
- `plot` - Generate visualization plots

## count - Data Processing

Generate pileup data from BAM files by counting stop and mutation signals.

### Usage
```bash
modtector count [OPTIONS] --bam <BAM> --fasta <FASTA> --output <OUTPUT>
```

### Options

| Option | Short | Type | Default | Description |
|--------|-------|------|---------|-------------|
| `--bam` | `-b` | String | Required | BAM file path |
| `--fasta` | `-f` | String | Required | Reference FASTA file path |
| `--output` | `-o` | String | Required | Output CSV path |
| `--threads` | `-t` | Integer | 1 | Number of parallel threads |
| `--log` | `-l` | String | None | Log file path (optional) |
| `--window` | `-w` | Integer | None | Window size for genome segmentation (bases) |

### Examples

```bash
# Basic usage
modtector count -b sample.bam -f reference.fa -o output.csv

# With multiple threads
modtector count -b sample.bam -f reference.fa -o output.csv -t 8

# With windowing
modtector count -b sample.bam -f reference.fa -o output.csv -w 1000

# With logging
modtector count -b sample.bam -f reference.fa -o output.csv -l count.log
```

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
| `--unmod` | `-U` | String | Required | Unmodified sample CSV |
| `--output` | `-O` | String | Required | Reactivity output file |
| `--stop-method` | `-s` | String | "current" | Stop signal reactivity method |
| `--mutation-method` | `-m` | String | "current" | Mutation signal reactivity method |
| `--pseudocount` | | Float | 1 | Pseudocount parameter for logarithmic calculations |
| `--maxscore` | | Float | 10 | Maximum score limit parameter |
| `--threads` | `-t` | Integer | 1 | Number of parallel threads |
| `--log` | `-l` | String | None | Log file path (optional) |

### Reactivity Methods

#### Stop Signal Methods
- **current**: Current implementation
- **ding**: Ding et al. method
- **rouskin**: Rouskin et al. method

#### Mutation Signal Methods
- **current**: Current implementation
- **siegfried**: Siegfried et al. method
- **zubradt**: Zubradt et al. method

### Examples

```bash
# Basic reactivity calculation
modtector reactivity -M mod.csv -U unmod.csv -O reactivity.csv

# With specific methods
modtector reactivity -M mod.csv -U unmod.csv -O reactivity.csv -s ding -m siegfried

# With custom parameters
modtector reactivity -M mod.csv -U unmod.csv -O reactivity.csv --pseudocount 0.5 --maxscore 5
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
| `--signal-type` | `-t` | String | "stop" | Signal type (stop, mutation, or both) |
| `--gene-id` | `-g` | String | Required | Gene ID (for base matching) |
| `--strand` | `-S` | String | "+" | Strand information (+ or -) |
| `--base-matching` | | Flag | True | Use base matching |
| `--auto-shift` | | Flag | True | Use auto-shift correction |
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
| `--svg-signal` | String | "all" | Signal type to plot (all, stop, mutation, score) |
| `--svg-strand` | String | "+" | Strand to include (+, -, both) |
| `--svg-ref` | String | None | Reference sequence file for alignment (optional) |
| `--svg-max-shift` | Integer | 5 | Maximum shift for alignment search |

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
- Open an issue on the project repository
