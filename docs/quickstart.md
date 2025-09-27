# Quick Start Guide

This guide will help you get started with modtector quickly. We'll walk through a complete example from raw data to final results.

## Prerequisites

- modtector installed (see [Installation Guide](installation.md))
- Sample BAM files (modified and unmodified samples)
- Reference FASTA file
- Secondary structure file (for evaluation)

## Example Dataset

For this quick start, we'll use the example data provided in the repository:

```bash
# Navigate to the example directory
cd example/

# Check available data
ls data/
ls ref/
```

## Step 1: Generate Pileup Data

First, we'll process the BAM files to generate pileup data:

```bash
# Process modified sample
modtector count \
    -b data/HEK293_mod_Human.bam \
    -f ref/Human_18S.fa \
    -o signal/01_count/mod_sample.csv \
    -t 4

# Process unmodified sample
modtector count \
    -b data/HEK293_unmod_Human.bam \
    -f ref/Human_18S.fa \
    -o signal/01_count/unmod_sample.csv \
    -t 4
```

**Parameters explained:**
- `-b`: Input BAM file
- `-f`: Reference FASTA file
- `-o`: Output CSV file
- `-t`: Number of threads for parallel processing

## Step 2: Calculate Reactivity Scores

Calculate reactivity scores by comparing modified and unmodified samples:

```bash
modtector reactivity \
    -M signal/01_count/mod_sample.csv \
    -U signal/01_count/unmod_sample.csv \
    -O signal/02_reactivity/reactivity.csv \
    -s current \
    -m current \
    -t 4
```

**Parameters explained:**
- `-M`: Modified sample CSV
- `-U`: Unmodified sample CSV
- `-O`: Output reactivity file
- `-s`: Stop signal method
- `-m`: Mutation signal method
- `-t`: Number of threads

## Step 3: Normalize Reactivity Signals

Normalize the reactivity data to remove noise and outliers:

```bash
modtector norm \
    -i signal/02_reactivity/reactivity.csv \
    -o signal/03_norm/normalized_reactivity.csv \
    -m winsor90
```

**Parameters explained:**
- `-M`: Modified sample CSV
- `-U`: Unmodified sample CSV
- `-O`: Output reactivity file
- `-s`: Stop signal method
- `-m`: Mutation signal method
- `-t`: Number of threads

## Step 4: Generate Visualizations

Create plots to visualize the results:

```bash
modtector plot \
    -M signal/01_count/mod_sample.csv \
    -U signal/01_count/unmod_sample.csv \
    -o signal/04_plot/ \
    -r signal/03_norm/normalized_reactivity.csv \
    -t 4
```

**Parameters explained:**
- `-M`: Modified sample CSV
- `-U`: Unmodified sample CSV
- `-o`: Output directory for plots
- `-r`: Reactivity file (optional)
- `-t`: Number of threads

## Step 5: Evaluate Accuracy

Evaluate the accuracy of your results using known secondary structure:

```bash
modtector evaluate \
    -r signal/03_norm/normalized_reactivity.csv \
    -s ref/Human_18S.dp \
    -o signal/05_evaluate/ \
    -g Human_18S \
    -S +
```

**Parameters explained:**
- `-r`: Reactivity file
- `-s`: Secondary structure file (.dp format)
- `-o`: Output directory
- `-g`: Gene ID
- `-S`: Strand information (+ or -)

## Expected Outputs

After running all steps, you should have:

```
signal/
├── 01_count/          # Raw pileup data
├── 02_reactivity/     # Reactivity scores
├── 03_norm/           # Normalized reactivity data
├── 04_plot/           # Visualization plots
└── 05_evaluate/       # Accuracy evaluation
```

### Key Output Files

1. **Pileup Data** (`01_count/*.csv`):
   - Raw signal counts for each position
   - Stop and mutation signals
   - Coverage information

2. **Reactivity Scores** (`02_reactivity/reactivity.csv`):
   - Calculated reactivity values
   - Signal differences between samples

3. **Normalized Reactivity** (`03_norm/normalized_reactivity.csv`):
   - Filtered and normalized reactivity signals
   - Outlier-corrected reactivity values

4. **Plots** (`04_plot/*.svg`):
   - Signal distribution plots
   - Reactivity visualization
   - Comparison charts

5. **Evaluation Results** (`05_evaluate/*.txt`):
   - AUC scores
   - F1-scores
   - Accuracy metrics
   - ROC/PR curves

## Understanding the Results

### Reactivity Scores
- **High positive values**: Likely modification sites
- **Near zero**: No significant modification
- **Negative values**: Possible artifacts or noise

### Evaluation Metrics
- **AUC > 0.7**: Good performance
- **F1-score > 0.6**: Reasonable accuracy
- **Accuracy > 0.8**: High confidence

## Common Issues and Solutions

### Issue 1: Low Coverage
**Problem**: Insufficient read depth
**Solution**: 
- Increase sequencing depth
- Adjust coverage thresholds in plotting

### Issue 2: Poor Normalization
**Problem**: High background noise
**Solution**:
- Try different normalization methods
- Adjust window sizes
- Filter low-quality positions

### Issue 3: Low Evaluation Scores
**Problem**: Poor accuracy metrics
**Solution**:
- Check data quality
- Verify reference sequence
- Ensure proper sample preparation

## Next Steps

Now that you've completed the quick start:

1. **Explore Advanced Features**:
   - Try different normalization methods
   - Experiment with various reactivity calculation methods
   - Use different statistical tests in comparison

2. **Analyze Your Own Data**:
   - Prepare your BAM files
   - Obtain reference sequences
   - Get secondary structure information

3. **Read the Full Documentation**:
   - [User Guide](user-guide.md) for detailed explanations
   - [Command Reference](commands.md) for all options
   - [Examples](examples.md) for more use cases

## Getting Help

If you encounter issues:

1. Check the [Troubleshooting Guide](troubleshooting.md)
2. Review the [Command Reference](commands.md)
3. Look at the [Examples](examples.md)
4. Open an issue on the project repository

## Tips for Success

1. **Start Small**: Begin with a small dataset to test your pipeline
2. **Check Data Quality**: Ensure your BAM files are properly aligned
3. **Use Appropriate References**: Match your reference sequence to your data
4. **Monitor Resources**: Large datasets may require significant memory
5. **Validate Results**: Always check your results against known modifications
