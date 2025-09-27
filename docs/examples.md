# Examples and Tutorials

This page provides practical examples and tutorials for using modtector in various scenarios.

## Table of Contents

1. [Basic Example](#basic-example)
2. [Advanced Example](#advanced-example)
3. [SVG Plotting Examples](#svg-plotting-examples)
4. [Batch Processing](#batch-processing)
5. [Custom Analysis](#custom-analysis)
6. [Troubleshooting Examples](#troubleshooting-examples)
7. [Performance Optimization](#performance-optimization)

## Basic Example

### Scenario
You have BAM files from modified and unmodified samples and want to detect RNA modifications.

### Data
- `mod_sample.bam`: Modified sample BAM file
- `unmod_sample.bam`: Unmodified sample BAM file
- `reference.fa`: Reference FASTA file
- `structure.dp`: Secondary structure file

### Step-by-Step Process

#### 1. Generate Pileup Data
```bash
# Process modified sample
modtector count \
    -b mod_sample.bam \
    -f reference.fa \
    -o mod_count.csv \
    -t 4

# Process unmodified sample
modtector count \
    -b unmod_sample.bam \
    -f reference.fa \
    -o unmod_count.csv \
    -t 4
```

#### 2. Calculate Reactivity
```bash
modtector reactivity \
    -M mod_count.csv \
    -U unmod_count.csv \
    -O reactivity.csv \
    -t 4
```

#### 3. Normalize Reactivity Signals
```bash
modtector norm \
    -i reactivity.csv \
    -o normalized_reactivity.csv \
    -m winsor90
```

#### 4. Generate Plots
```bash
modtector plot \
    -M mod_count.csv \
    -U unmod_count.csv \
    -o plots/ \
    -r normalized_reactivity.csv \
    -t 4
```

#### 5. Evaluate Accuracy
```bash
modtector evaluate \
    -r normalized_reactivity.csv \
    -s structure.dp \
    -o results/ \
    -g gene_id
```

### Expected Outputs
```
output/
├── mod_count.csv          # Raw pileup data (modified)
├── unmod_count.csv        # Raw pileup data (unmodified)
├── reactivity.csv         # Reactivity scores
├── normalized_reactivity.csv  # Normalized reactivity data
├── plots/                 # Visualization plots
│   ├── signal_distribution.svg
│   ├── reactivity_plot.svg
│   └── comparison_plot.svg
└── results/               # Evaluation results
    ├── evaluation_comprehensive.txt
    ├── roc_curve.svg
    └── pr_curve.svg
```

## Advanced Example

### Scenario
You have multiple samples and want to perform comprehensive analysis with different methods and parameters.

### Data
- Multiple BAM files for different conditions
- High-quality reference sequences
- Experimentally determined secondary structures

### Advanced Workflow

#### 1. Batch Processing with Windowing
```bash
# Process multiple samples with windowing
for sample in mod1 mod2 mod3 unmod1 unmod2 unmod3; do
    modtector count \
        -b ${sample}.bam \
        -f reference.fa \
        -o ${sample}_count.csv \
        -w 1000 \
        -t 8
done
```

#### 2. Advanced Normalization
```bash
# Use dynamic windowing for normalization
for sample in mod1 mod2 mod3 unmod1 unmod2 unmod3; do
    modtector norm \
        -i ${sample}_count.csv \
        -o ${sample}_norm.csv \
        -m winsor90 \
        --dynamic \
        --bases AC \
        --linear
done
```

#### 3. Multiple Reactivity Methods
```bash
# Compare different reactivity calculation methods
modtector reactivity \
    -M mod1_norm.csv \
    -U unmod1_norm.csv \
    -O reactivity_current.csv \
    -s current \
    -m current

modtector reactivity \
    -M mod1_norm.csv \
    -U unmod1_norm.csv \
    -O reactivity_ding.csv \
    -s ding \
    -m siegfried \
    --pseudocount 0.5 \
    --maxscore 5
```

#### 4. Statistical Comparison
```bash
# Compare different samples
modtector compare \
    -M mod1_norm.csv \
    -U unmod1_norm.csv \
    -o comparison_t_test.csv \
    -t t-test \
    -d 20 \
    -f 1.5

modtector compare \
    -M mod1_norm.csv \
    -U unmod1_norm.csv \
    -o comparison_mann_whitney.csv \
    -t mann-whitney \
    -d 20 \
    -f 1.5
```

#### 5. Comprehensive Evaluation
```bash
# Evaluate with different parameters
modtector evaluate \
    -r reactivity_current.csv \
    -s structure.dp \
    -o results_current/ \
    -g gene_id \
    --auto-shift \
    --base-matching

modtector evaluate \
    -r reactivity_ding.csv \
    -s structure.dp \
    -o results_ding/ \
    -g gene_id \
    --auto-shift \
    --base-matching
```

## SVG Plotting Examples

### Scenario
You want to visualize reactivity data on RNA secondary structure using SVG plots.

### Data Requirements
- Reactivity CSV file with position, base, and score data
- SVG template file with RNA structure
- Optional: Reference sequence for alignment

### Basic SVG Plotting

#### 1. Simple SVG Plot (SVG-Only Mode)
```bash
# Basic SVG plotting without regular plots
modtector plot \
    -o svg_output/ \
    --svg-template rna_structure.svg \
    --reactivity reactivity_data.csv
```

#### 2. Multi-Signal SVG Plotting
```bash
# Plot all signal types (stop, mutation, etc.)
modtector plot \
    -o svg_output/ \
    --svg-template rna_structure.svg \
    --reactivity reactivity_data.csv \
    --svg-bases ATGC \
    --svg-signal all \
    --svg-strand +
```

#### 3. Single Signal SVG Plotting
```bash
# Plot only stop signals
modtector plot \
    -o svg_output/ \
    --svg-template rna_structure.svg \
    --reactivity reactivity_data.csv \
    --svg-bases AC \
    --svg-signal stop \
    --svg-strand +
```

### Advanced SVG Plotting

#### 4. SVG Plotting with Alignment
```bash
# Use reference sequence for alignment
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

#### 5. Negative Strand SVG Plotting
```bash
# Plot negative strand data
modtector plot \
    -o svg_output/ \
    --svg-template rna_structure.svg \
    --reactivity reactivity_data.csv \
    --svg-bases ATGC \
    --svg-signal stop \
    --svg-strand -
```

#### 6. Both Strands SVG Plotting
```bash
# Plot both positive and negative strands
modtector plot \
    -o svg_output/ \
    --svg-template rna_structure.svg \
    --reactivity reactivity_data.csv \
    --svg-bases ATGC \
    --svg-signal all \
    --svg-strand both
```

### Combined Mode (Regular + SVG)

#### 7. Generate Both Regular and SVG Plots
```bash
# Generate both regular plots and SVG plots
modtector plot \
    -M mod_sample.csv \
    -U unmod_sample.csv \
    -o combined_output/ \
    --svg-template rna_structure.svg \
    --reactivity reactivity_data.csv \
    --svg-bases ATGC \
    --svg-signal all
```

### Real-World Examples

#### Example 1: Human 18S rRNA Structure
```bash
# Plot Human 18S rRNA with multiple signals
modtector plot \
    -o human_18s_output/ \
    --svg-template Human_18S.svg \
    --reactivity md_HEK293_Human_norm.csv \
    --svg-bases ATGC \
    --svg-signal all \
    --svg-strand +
```

**Output files:**
- `rna_structure_colored_stop.svg` - Stop signal visualization
- `rna_structure_colored_mutation.svg` - Mutation signal visualization

#### Example 2: E.coli 16S rRNA Structure
```bash
# Plot E.coli 16S rRNA with alignment
modtector plot \
    -o ecoli_16s_output/ \
    --svg-template Ecoli_RNA_structure.svg \
    --reactivity EC_DMS-seq.csv \
    --svg-bases AC \
    --svg-signal all \
    --svg-strand + \
    --svg-ref Ecoli_RNA_structure.svg \
    --svg-max-shift 5
```

**Output files:**
- `rna_structure_colored_score.svg` - Single signal visualization

#### Example 3: DMS-seq Specific Analysis
```bash
# DMS-seq analysis (A and C bases only)
modtector plot \
    -o dms_seq_output/ \
    --svg-template rna_structure.svg \
    --reactivity dms_data.csv \
    --svg-bases AC \
    --svg-signal stop \
    --svg-strand +
```

### SVG Output Files

#### File Naming Convention
- **Multi-signal**: `rna_structure_colored_[signal_type].svg`
  - `rna_structure_colored_stop.svg`
  - `rna_structure_colored_mutation.svg`
- **Single-signal**: `rna_structure_colored_score.svg`

#### SVG File Structure
- **Original template**: Complete RNA structure template
- **Colored circles**: Reactivity scores mapped to positions
- **Color bars**: Legend showing score ranges
- **Multiple SVG elements**: One for each signal type

### Troubleshooting SVG Plotting

#### Common Issues

1. **No circles in SVG output**
   ```bash
   # Check CSV format and position mapping
   head -5 reactivity_data.csv
   
   # Verify SVG template has position labels
   grep -c "title>" rna_structure.svg
   ```

2. **Wrong signal type plotted**
   ```bash
   # Check CSV headers for signal columns
   head -1 reactivity_data.csv
   
   # Use specific signal type
   --svg-signal stop  # or mutation, score
   ```

3. **Missing bases in output**
   ```bash
   # Check base filtering
   --svg-bases ATGC  # All bases
   --svg-bases AC    # A and C only
   ```

#### Quality Control
```bash
# Check SVG file size and content
ls -la svg_output/
wc -l svg_output/*.svg
grep -c "<circle" svg_output/*.svg
```

### Best Practices

1. **Use appropriate base filtering**: Match your experimental method
   - DMS-seq: `--svg-bases AC`
   - SHAPE-seq: `--svg-bases ATGC`

2. **Choose correct signal type**: Based on your data
   - Single signal: `--svg-signal score`
   - Multiple signals: `--svg-signal all`

3. **Use alignment when needed**: For sequence mismatches
   - Reference sequence: `--svg-ref reference.fa`
   - Max shift: `--svg-max-shift 10`

4. **Optimize for your data**: Adjust parameters
   - Strand selection: `--svg-strand +` (default)
   - Base filtering: Match experimental method

## Batch Processing

### Scenario
You have many samples and want to process them efficiently.

### Batch Script Example
```bash
#!/bin/bash

# Configuration
BAM_DIR="/path/to/bam/files"
REFERENCE="/path/to/reference.fa"
OUTPUT_DIR="/path/to/output"
THREADS=8

# Create output directories
mkdir -p ${OUTPUT_DIR}/{count,norm,reactivity,plots,results}

# Process all BAM files
for bam_file in ${BAM_DIR}/*.bam; do
    sample_name=$(basename ${bam_file} .bam)
    
    echo "Processing ${sample_name}..."
    
    # Generate pileup data
    modtector count \
        -b ${bam_file} \
        -f ${REFERENCE} \
        -o ${OUTPUT_DIR}/count/${sample_name}_count.csv \
        -t ${THREADS}
    
    # Normalize signals
    modtector norm \
        -i ${OUTPUT_DIR}/count/${sample_name}_count.csv \
        -o ${OUTPUT_DIR}/norm/${sample_name}_norm.csv \
        -m winsor90 \
        --bases AC
done

# Calculate reactivity for all pairs
for mod_sample in ${OUTPUT_DIR}/norm/*_mod_*.csv; do
    unmod_sample=${mod_sample/_mod_/_unmod_}
    
    if [ -f "${unmod_sample}" ]; then
        sample_name=$(basename ${mod_sample} _norm.csv)
        
        echo "Calculating reactivity for ${sample_name}..."
        
        modtector reactivity \
            -M ${mod_sample} \
            -U ${unmod_sample} \
            -O ${OUTPUT_DIR}/reactivity/${sample_name}_reactivity.csv \
            -t ${THREADS}
    fi
done
```

### Parallel Processing
```bash
#!/bin/bash

# Use GNU parallel for maximum efficiency
parallel -j 4 modtector count \
    -b {} \
    -f reference.fa \
    -o {.}_count.csv \
    -t 2 ::: *.bam
```

## Custom Analysis

### Scenario
You want to analyze specific regions or use custom parameters.

### Region-Specific Analysis
```bash
# Extract specific regions first
samtools view -b mod_sample.bam chr1:1000-2000 > mod_region.bam
samtools view -b unmod_sample.bam chr1:1000-2000 > unmod_region.bam

# Process specific regions
modtector count \
    -b mod_region.bam \
    -f reference.fa \
    -o mod_region_count.csv \
    -t 4

modtector count \
    -b unmod_region.bam \
    -f reference.fa \
    -o unmod_region_count.csv \
    -t 4
```

### Custom Normalization Parameters
```bash
# Use custom window size and offset
modtector norm \
    -i input.csv \
    -o output.csv \
    -m winsor90 \
    --window 500 \
    --window-offset 100 \
    --bases AC
```

### Custom Reactivity Parameters
```bash
# Use custom pseudocount and maxscore
modtector reactivity \
    -M mod.csv \
    -U unmod.csv \
    -O reactivity.csv \
    --pseudocount 0.1 \
    --maxscore 20
```

## Troubleshooting Examples

### Low Coverage Issue
```bash
# Check coverage
samtools depth mod_sample.bam | awk '{sum+=$3} END {print sum/NR}'

# If coverage is low, increase depth threshold
modtector plot \
    -M mod_norm.csv \
    -U unmod_norm.csv \
    -o plots/ \
    -d 20  # Lower depth threshold
```

### Memory Issues
```bash
# Reduce thread count
modtector count \
    -b large_sample.bam \
    -f reference.fa \
    -o output.csv \
    -t 2  # Use fewer threads

# Use windowing to reduce memory usage
modtector count \
    -b large_sample.bam \
    -f reference.fa \
    -o output.csv \
    -w 1000  # Use smaller windows
```

### Poor Normalization
```bash
# Try different normalization methods
modtector norm \
    -i input.csv \
    -o output_percentile.csv \
    -m percentile28

modtector norm \
    -i input.csv \
    -o output_boxplot.csv \
    -m boxplot

# Compare results
modtector plot \
    -M output_percentile.csv \
    -U output_boxplot.csv \
    -o comparison_plots/
```

### Low Evaluation Scores
```bash
# Check data quality
modtector evaluate \
    -r reactivity.csv \
    -s structure.dp \
    -o results/ \
    -g gene_id \
    --no-auto-shift  # Try without auto-shift

# Use different signal types
modtector evaluate \
    -r reactivity.csv \
    -s structure.dp \
    -o results/ \
    -g gene_id \
    -t mutation  # Use only mutation signals
```

## Performance Optimization

### Large Dataset Processing
```bash
# Use maximum threads
modtector count \
    -b large_sample.bam \
    -f reference.fa \
    -o output.csv \
    -t 16  # Use all available cores

# Use windowing for large genomes
modtector count \
    -b large_sample.bam \
    -f reference.fa \
    -o output.csv \
    -w 5000  # Use larger windows
```

### Memory Optimization
```bash
# Process in chunks
split -l 10000 large_input.csv chunk_
for chunk in chunk_*; do
    modtector norm \
        -i ${chunk} \
        -o ${chunk}_norm.csv \
        -m winsor90
done
```

### Disk Space Management
```bash
# Compress intermediate files
gzip *.csv

# Remove intermediate files after processing
rm *_count.csv *_norm.csv

# Keep only final results
mkdir final_results
mv reactivity.csv plots/ results/ final_results/
```

## Real-World Examples

### Example 1: m6A Detection
```bash
# m6A-specific analysis
modtector norm \
    -i input.csv \
    -o output.csv \
    -m winsor90 \
    --bases A  # Focus on A bases

modtector reactivity \
    -M mod.csv \
    -U unmod.csv \
    -O m6A_reactivity.csv \
    -s ding \
    --pseudocount 0.5
```

### Example 2: m1A Detection
```bash
# m1A-specific analysis
modtector norm \
    -i input.csv \
    -o output.csv \
    -m winsor90 \
    --bases A  # Focus on A bases

modtector reactivity \
    -M mod.csv \
    -U unmod.csv \
    -O m1A_reactivity.csv \
    -s current \
    --maxscore 5
```

### Example 3: Pseudouridine Detection
```bash
# Pseudouridine-specific analysis
modtector norm \
    -i input.csv \
    -o output.csv \
    -m winsor90 \
    --bases U  # Focus on U bases

modtector reactivity \
    -M mod.csv \
    -U unmod.csv \
    -O psi_reactivity.csv \
    -m siegfried
```

## Quality Control Examples

### Data Quality Assessment
```bash
# Check BAM file quality
samtools flagstat mod_sample.bam
samtools view -c mod_sample.bam

# Check pileup data quality
head -100 mod_count.csv
wc -l mod_count.csv

# Check normalization results
modtector plot \
    -M mod_norm.csv \
    -U unmod_norm.csv \
    -o qc_plots/ \
    -c 0.1 \
    -d 10
```

### Result Validation
```bash
# Compare different methods
modtector reactivity \
    -M mod.csv \
    -U unmod.csv \
    -O reactivity_method1.csv \
    -s current

modtector reactivity \
    -M mod.csv \
    -U unmod.csv \
    -O reactivity_method2.csv \
    -s ding

# Compare results
modtector compare \
    --mode reactivity-results \
    --reactivity1 reactivity_method1.csv \
    --reactivity2 reactivity_method2.csv \
    -o method_comparison.csv
```

## Integration Examples

### Snakemake Workflow
```python
rule count:
    input:
        bam = "data/{sample}.bam",
        fasta = "ref/reference.fa"
    output:
        csv = "results/count/{sample}_count.csv"
    shell:
        "modtector count -b {input.bam} -f {input.fasta} -o {output.csv} -t {threads}"

rule norm:
    input:
        csv = "results/count/{sample}_count.csv"
    output:
        csv = "results/norm/{sample}_norm.csv"
    shell:
        "modtector norm -i {input.csv} -o {output.csv} -m winsor90 --bases AC"

rule reactivity:
    input:
        mod = "results/norm/{sample}_mod_norm.csv",
        unmod = "results/norm/{sample}_unmod_norm.csv"
    output:
        csv = "results/reactivity/{sample}_reactivity.csv"
    shell:
        "modtector reactivity -M {input.mod} -U {input.unmod} -O {output.csv}"
```

### Nextflow Workflow
```groovy
process count {
    input:
    file bam
    file fasta
    
    output:
    file "*.csv"
    
    script:
    """
    modtector count -b ${bam} -f ${fasta} -o output.csv -t ${task.cpus}
    """
}

process norm {
    input:
    file csv
    
    output:
    file "*.csv"
    
    script:
    """
    modtector norm -i ${csv} -o output.csv -m winsor90 --bases AC
    """
}
```

These examples demonstrate various ways to use modtector for different scenarios and requirements. Choose the approach that best fits your specific needs and data characteristics.
