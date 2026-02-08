# modtector Documentation

Welcome to modtector, a high-performance RNA modification detection tool built with Rust.

```{toctree}
:maxdepth: 2
:hidden:

overview
installation
quickstart
user-guide
commands
examples
troubleshooting
contributing
```

## Welcome

Welcome to the modtector documentation! This is your starting point for learning about our high-performance RNA modification detection tool.

## Getting Started

- **New to modtector?** Start with the [Overview](overview.md) to understand the tool's capabilities and browse the documentation structure
- **Need installation help?** Check the [Installation Guide](installation.md)
- **Ready to get started?** Follow the [Quick Start Guide](quickstart.md)

### Quick Links

- [Installation](installation.md) - Detailed installation instructions
- [Quick Start Guide](quickstart.md) - Get started with modtector
- [User Guide](user-guide.md) - Comprehensive usage guide
- [Command Reference](commands.md) - Complete command documentation
- [Examples](examples.md) - Practical examples and tutorials
- [Troubleshooting](troubleshooting.md) - Common issues and solutions
- [Contributing](contributing.md) - How to contribute to modtector

## Workflow Overview

modtector provides a complete workflow from raw data to evaluation results:

![Workflow Diagram](images/workflow.png)

1. **Data Input**: BAM files (modified/unmodified samples), FASTA reference sequences, secondary structure files
2. **Statistical Analysis**: Pileup traversal, counting stop and mutation signals
3. **Data Normalization**: Signal filtering, outlier handling, background correction
4. **Comparative Analysis**: Modified vs unmodified sample comparison, identifying modification sites
5. **Reactivity Calculation**: Calculate signal differences, generate reactivity data
6. **Duet Analysis**: Sliding-window inference of dynamic ensembles from normalized reactivity and read-level co-variation
7. **Visualization**: Generate signal distribution plots and reactivity plots
8. **Accuracy Assessment**: Performance evaluation based on secondary structure

## Version Information

- **Current Version**: v0.15.4
- **Release Date**: 2026-02-08
- **License**: MIT License

## Project Links

- **GitHub Repository**: [https://github.com/TongZhou2017/modtector](https://github.com/TongZhou2017/modtector)
- **Crates.io**: [https://crates.io/crates/modtector](https://crates.io/crates/modtector)
- **Documentation**: [https://modtector.readthedocs.io/](https://modtector.readthedocs.io/)

### Recent Highlights

- **PCR Bias Correction**: Chi-Square distribution-based depth correction (`correct` command)
- **Base Quality Filtering**: Per-base quality filtering for mutation detection (`--min-base-qual`)
- **Extended Format Support**: Support for multiple input formats in `convert` command (rf-rctools, shapemapper2, icSHAPE-rt, bedgraph, etc.)
- **Zarringhalam Remap Improvements**: Uses actual maximum values instead of fixed 1.0 for better high-value region mapping
- **Distribution-Based K-Factor Prediction**: Advanced k-factor prediction method using statistical distribution analysis
- **Batch Processing**: Process multiple BAM files sequentially with glob patterns
- **Single-cell Unified Processing**: Unified processing with cell label extraction for 2-3x performance improvement

## Example Dataset

To help you quickly evaluate modtector, we provide a minimal example dataset available at [Zenodo (10.5281/zenodo.17316476)](https://doi.org/10.5281/zenodo.17316476). This dataset contains a small subset of data that can be processed quickly to demonstrate modtector's functionality.

### Quick Start with Example Data

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
- Duet sliding-window ensemble decomposition using normalized reactivity and BAM co-variation
- Visualization generation
- Accuracy evaluation

The results will be organized in the `signal_v0.5.6/` directory with subdirectories for each processing step.

## Support

For questions, issues, or contributions, please refer to the [Contributing Guide](contributing.md) or open an issue on the [GitHub repository](https://github.com/TongZhou2017/modtector).

## Citation

If you use modtector in your research, please cite:

```bibtex
[Add citation information when available]
```