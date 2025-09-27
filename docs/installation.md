# Installation Guide

This guide provides detailed instructions for installing modtector on different operating systems.

## System Requirements

### Minimum Requirements
- **Operating System**: Linux, macOS, or Windows
- **Memory**: 4 GB RAM (8 GB recommended for large datasets)
- **Storage**: 2 GB free space
- **CPU**: Multi-core processor recommended for parallel processing

### Required Dependencies
- **Rust**: Version 1.70 or higher
- **Cargo**: Included with Rust installation

## Installing Rust

### Linux/macOS

1. **Install Rust using rustup (recommended)**:
```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
source ~/.cargo/env
```

2. **Verify installation**:
```bash
rustc --version
cargo --version
```

### Windows

1. **Download and run rustup installer**:
   - Visit https://rustup.rs/
   - Download and run `rustup-init.exe`
   - Follow the installation prompts

2. **Verify installation**:
```cmd
rustc --version
cargo --version
```

## Installing modtector

### Method 1: Build from Source (Recommended)

1. **Clone the repository**:
```bash
git clone https://github.com/your-org/modtector.git
cd modtector
```

2. **Build the project**:
```bash
# Debug build (faster compilation)
cargo build

# Release build (optimized performance)
cargo build --release
```

3. **Verify installation**:
```bash
./target/release/modtector --version
./target/release/modtector --help
```

### Method 2: Install via Cargo (if published)

```bash
cargo install modtector
```

## Dependencies

modtector uses the following Rust crates:

```toml
[dependencies]
rust-htslib = "0.47.0"    # BAM file processing
bio = "2.0.1"             # Bioinformatics utilities
plotters = "0.3"          # Plotting and visualization
csv = "1.3"               # CSV file handling
rayon = "1.8"             # Parallel processing
clap = "4.5"              # Command-line argument parsing
chrono = "0.4"            # Date and time handling
rand = "0.8"              # Random number generation
```

These dependencies will be automatically downloaded and compiled during the build process.

## Platform-Specific Notes

### Linux

- **Ubuntu/Debian**: May require additional packages:
```bash
sudo apt-get update
sudo apt-get install build-essential pkg-config libssl-dev
```

- **CentOS/RHEL**: May require:
```bash
sudo yum groupinstall "Development Tools"
sudo yum install pkgconfig openssl-devel
```

### macOS

- **Xcode Command Line Tools**: Required for compilation:
```bash
xcode-select --install
```

### Windows

- **Visual Studio Build Tools**: Required for compilation
- **Git**: Required for cloning the repository

## Verification

After installation, verify that modtector is working correctly:

```bash
# Check version
modtector --version

# Check available commands
modtector --help

# Test with a simple command
modtector count --help
```

## Troubleshooting

### Common Issues

1. **Rust not found**:
   - Ensure Rust is properly installed and in your PATH
   - Run `source ~/.cargo/env` (Linux/macOS) or restart terminal

2. **Build failures**:
   - Update Rust: `rustup update`
   - Clean build: `cargo clean && cargo build --release`

3. **Permission errors**:
   - Ensure you have write permissions in the installation directory
   - Use `sudo` if necessary (not recommended for development)

4. **Memory issues during compilation**:
   - Increase swap space or available RAM
   - Use `cargo build` instead of `cargo build --release` for initial testing

### Getting Help

If you encounter issues during installation:

1. Check the [Troubleshooting Guide](troubleshooting.md)
2. Search existing issues on the project repository
3. Create a new issue with detailed error information

## Next Steps

After successful installation:

1. Read the [Quick Start Guide](quickstart.md)
2. Explore the [User Guide](user-guide.md)
3. Try the [Examples](examples.md)

## Uninstallation

To remove modtector:

1. **If installed via cargo**:
```bash
cargo uninstall modtector
```

2. **If built from source**:
```bash
# Simply delete the cloned directory
rm -rf modtector
```

3. **Remove Rust (optional)**:
```bash
rustup self uninstall
```
