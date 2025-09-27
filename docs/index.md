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

1. **Data Input**: BAM files (modified/unmodified samples), FASTA reference sequences, secondary structure files
2. **Statistical Analysis**: Pileup traversal, counting stop and mutation signals
3. **Data Normalization**: Signal filtering, outlier handling, background correction
4. **Comparative Analysis**: Modified vs unmodified sample comparison, identifying modification sites
5. **Reactivity Calculation**: Calculate signal differences, generate reactivity data
6. **Visualization**: Generate signal distribution plots and reactivity plots
7. **Accuracy Assessment**: Performance evaluation based on secondary structure

## Version Information

- **Current Version**: v0.9.6
- **Release Date**: 2025-01-27
- **License**: MIT License

## Support

For questions, issues, or contributions, please refer to the [Contributing Guide](contributing.md) or open an issue on the project repository.

## Citation

If you use modtector in your research, please cite:

```bibtex
[Add citation information when available]
```
