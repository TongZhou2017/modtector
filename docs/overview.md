# Overview

modtector is a comprehensive tool for detecting RNA modifications from sequencing data. It provides a complete workflow from raw BAM files to final evaluation results, supporting multiple signal types, normalization methods, and accuracy assessment.

## Key Features

- **Multi-signal Analysis**: Simultaneous analysis of stop signals (pipeline truncation) and mutation signals (base mutations)
- **High-performance Pileup**: Robust pileup traversal based on htslib, supporting large file processing
- **Data Normalization**: Signal validity filtering and outlier handling
- **Reactivity Calculation**: Calculate signal differences between modified and unmodified samples with multi-threading support
- **Accuracy Assessment**: AUC, F1-score and other metrics evaluation based on secondary structure with auto-shift correction
- **Rich Visualization**: Signal distribution plots, reactivity plots, ROC/PR curves, and RNA structure SVG plots
- **Multi-threading**: Parallel processing and plotting for improved efficiency
- **Complete Workflow**: One-stop solution from raw data to final evaluation
- **Base Matching**: Intelligent base matching algorithm handling T/U equivalence
- **Auto-alignment**: Auto-shift correction for sequence length differences

## Quick Start

### Installation

1. Ensure Rust and Cargo are installed on your system
2. Clone the repository and build:

```bash
git clone <repository-url>
cd modtector
cargo build --release
```

### Basic Usage

```bash
# Generate pileup data from BAM files
modtector count -b sample.bam -f reference.fa -o output.csv

# Calculate reactivity scores
modtector reactivity -M mod_sample.csv -U unmod_sample.csv -O reactivity.csv

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
```

## Workflow Overview

modtector provides a complete workflow from raw data to evaluation results:

1. **Data Input**: BAM files (modified/unmodified samples), FASTA reference sequences, secondary structure files
2. **Statistical Analysis**: Pileup traversal, counting stop and mutation signals
3. **Reactivity Calculation**: Calculate signal differences between modified and unmodified samples, generate reactivity data
4. **Data Normalization**: Normalize reactivity signals, signal filtering, outlier handling, background correction
5. **Comparative Analysis**: Compare modified vs unmodified samples, identify differential modification sites
6. **Visualization**: Generate signal distribution plots, reactivity plots, and RNA structure SVG plots
7. **Accuracy Assessment**: Performance evaluation based on secondary structure

## Version Information

- **Current Version**: v0.9.5
- **Release Date**: 2025-09-27
- **License**: [Add license information]

## Support

For questions, issues, or contributions, please refer to the [Contributing Guide](contributing.md) or open an issue on the project repository.

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
