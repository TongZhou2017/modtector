# Overview

modtector is a comprehensive tool for detecting RNA modifications from sequencing data. It provides a complete workflow from raw BAM files to final evaluation results, supporting multiple signal types, normalization methods, and accuracy assessment.

## Key Features

- **Multi-signal Analysis**: Simultaneous analysis of stop signals (pipeline truncation) and mutation signals (base mutations)
- **High-performance Pileup**: Robust pileup traversal based on htslib, supporting large file processing
- **Batch Processing**: Process multiple BAM files sequentially with glob pattern matching
- **Single-cell Unified Processing**: Unified processing strategy for single-cell data with automatic cell label extraction (v0.11.0)
- **Data Distribution Optimization**: Smart data distribution scanning for efficient batch processing (v0.11.0)
- **Data Normalization**: Signal validity filtering and outlier handling
- **Reactivity Calculation**: Calculate signal differences between modified and unmodified samples with multi-threading support
- **Duet Ensemble Decomposition**: Sliding-window inference of dynamic RNA ensembles by combining normalized reactivity with read-level stop/mutation co-variation
- **Accuracy Assessment**: AUC, F1-score and other metrics evaluation based on secondary structure with auto-shift correction
- **Rich Visualization**: Signal distribution plots, reactivity plots, ROC/PR curves, RNA structure SVG plots, and interactive HTML visualizations
- **Multi-threading**: Parallel processing and plotting for improved efficiency
- **Complete Workflow**: One-stop solution from raw data to final evaluation
- **Base Matching**: Intelligent base matching algorithm handling T/U equivalence
- **Auto-alignment**: Auto-shift correction for sequence length differences

## Quick Start

### Installation

1. Ensure Rust and Cargo are installed on your system
2. Clone the repository and build:

```bash
git clone https://github.com/TongZhou2017/modtector.git
cd modtector
cargo build --release
```

Or install directly from [crates.io](https://crates.io/crates/modtector):

```bash
cargo install modtector
```

### Example Dataset

To help you quickly get started with modtector, we provide a minimal example dataset available at [Zenodo (10.5281/zenodo.17316476)](https://doi.org/10.5281/zenodo.17316476). This dataset contains a small subset of data that can be processed quickly to demonstrate modtector's functionality.

1. Download and extract the example dataset:
   ```bash
   wget https://zenodo.org/record/17316476/files/modtector_example_dataset.zip
   unzip modtector_example_dataset.zip
   cd modtector_example_dataset
   ```

2. Run the example script:
   ```bash
   bash test_modtector_v0.10.0.sh
   ```

This will execute a complete workflow on the example data, including:
- Pileup processing of BAM files
- Reactivity calculation with multiple methods
- Data normalization
- Duet ensemble analysis of normalized reactivity
- Visualization generation
- Accuracy evaluation

### Basic Usage

```bash
# Generate pileup data from BAM files
modtector count -b sample.bam -f reference.fa -o output.csv

# Calculate reactivity scores
modtector reactivity -M mod_sample.csv -U unmod_sample.csv -O reactivity.csv
# For experiments without an unmodified control (e.g., smartSHAPE), omit -U to run in mod-only mode
# modtector reactivity -M smartshape_mod.csv -O reactivity.csv

# Normalize reactivity signals
modtector norm -i reactivity.csv -o normalized_reactivity.csv -m winsor90

# Compare samples and identify differences
modtector compare -m mod_sample.csv -u unmod_sample.csv -o comparison.csv

# Generate visualizations
modtector plot -M mod_sample.csv -U unmod_sample.csv -o plots/ -r normalized_reactivity.csv

# Generate RNA structure SVG plots
modtector plot -o svg_output/ --svg-template rna_structure.svg --reactivity normalized_reactivity.csv

# Evaluate accuracy
modtector evaluate -r normalized_reactivity.csv -s structure.dp -o results/ -g gene_id

# Duet ensemble decomposition (Stop + Mutation co-variation, 100-nt windows)
modtector duet \
  -i normalized_reactivity.csv \
  -b sample.sort.bam \
  -f reference.fa \
  -o duet_windows.csv \
  --epsilon 0.8 \
  --min-samples 8 \
  --window-size 100 \
  --window-step 50
# Produces `<output>.csv`, `<output>_summary.csv`, `<output>_global.csv`, and `<output>_global_per_base.csv`
```

## Workflow Overview

modtector provides a complete workflow from raw data to evaluation results:

![Workflow Diagram](images/workflow.png)

1. **Data Input**: BAM files (modified/unmodified samples), FASTA reference sequences, secondary structure files
2. **Statistical Analysis**: Pileup traversal, counting stop and mutation signals
3. **Reactivity Calculation**: Calculate signal differences between modified and unmodified samples, generate reactivity data
4. **Data Normalization**: Normalize reactivity signals, signal filtering, outlier handling, background correction
5. **Duet Analysis**: Sliding-window decomposition of normalized reactivity into dynamic ensembles using read-level co-variation
6. **Comparative Analysis**: Compare modified vs unmodified samples, identify differential modification sites
7. **Visualization**: Generate signal distribution plots, reactivity plots, and RNA structure SVG plots
8. **Accuracy Assessment**: Performance evaluation based on secondary structure

## Version Information

- **Current Version**: v0.15.4
- **Release Date**: 2026-02-08
- **License**: MIT License

## Support

For questions, issues, or contributions, please refer to the [Contributing Guide](contributing.md) or open an issue on the [GitHub repository](https://github.com/TongZhou2017/modtector).

## Project Links

- **GitHub Repository**: [https://github.com/TongZhou2017/modtector](https://github.com/TongZhou2017/modtector)
- **Crates.io**: [https://crates.io/crates/modtector](https://crates.io/crates/modtector)
- **Documentation**: [https://modtector.readthedocs.io/](https://modtector.readthedocs.io/)

## Citation

If you use modtector in your research, please cite:

```bibtex
[Add citation information when available]
```

## Documentation Contents

This section provides a complete overview of the modtector documentation structure.

### Documentation Sections

#### Getting Started
- **[Installation](installation.md)** - Detailed installation instructions and system requirements
- **[Quick Start Guide](quickstart.md)** - Get started with modtector in minutes

#### User Documentation
- **[User Guide](user-guide.md)** - Comprehensive usage guide and workflow explanation
- **[Command Reference](commands.md)** - Complete command documentation and options
- **[Examples](examples.md)** - Practical examples and tutorials

#### Support & Development
- **[Troubleshooting](troubleshooting.md)** - Common issues and solutions
- **[Contributing](contributing.md)** - How to contribute to modtector

### Quick Navigation

- **New to modtector?** Start with [Installation](installation.md) and [Quick Start](quickstart.md)
- **Need help with commands?** Check the [Command Reference](commands.md)
- **Looking for examples?** Browse the [Examples](examples.md) section
- **Having issues?** See [Troubleshooting](troubleshooting.md)
- **Want to contribute?** Read the [Contributing Guide](contributing.md)