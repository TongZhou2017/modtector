# Troubleshooting Guide

This guide helps you diagnose and resolve common issues when using modtector.

## Table of Contents

1. [Installation Issues](#installation-issues)
2. [Data Processing Issues](#data-processing-issues)
3. [Performance Issues](#performance-issues)
4. [Output Quality Issues](#output-quality-issues)
5. [Memory and Resource Issues](#memory-and-resource-issues)
6. [File Format Issues](#file-format-issues)
7. [Getting Help](#getting-help)

## Installation Issues

### Rust Not Found

**Problem**: `rustc: command not found` or `cargo: command not found`

**Solutions**:
1. **Install Rust**:
   ```bash
   curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
   source ~/.cargo/env
   ```

2. **Check PATH**:
   ```bash
   echo $PATH
   export PATH="$HOME/.cargo/bin:$PATH"
   ```

3. **Verify installation**:
   ```bash
   rustc --version
   cargo --version
   ```

### Build Failures

**Problem**: `cargo build` fails with compilation errors

**Solutions**:
1. **Update Rust**:
   ```bash
   rustup update
   ```

2. **Clean and rebuild**:
   ```bash
   cargo clean
   cargo build --release
   ```

3. **Check dependencies**:
   ```bash
   cargo check
   ```

4. **Install build tools** (Linux):
   ```bash
   sudo apt-get install build-essential pkg-config libssl-dev
   ```

### Permission Errors

**Problem**: Permission denied during installation

**Solutions**:
1. **Check permissions**:
   ```bash
   ls -la target/
   ```

2. **Fix permissions**:
   ```bash
   chmod -R 755 target/
   ```

3. **Use user directory**:
   ```bash
   cargo install --path . --force
   ```

## Data Processing Issues

### BAM File Issues

**Problem**: BAM file cannot be read or processed

**Solutions**:
1. **Check BAM file**:
   ```bash
   samtools view -H sample.bam
   samtools flagstat sample.bam
   ```

2. **Validate BAM file**:
   ```bash
   samtools quickcheck sample.bam
   ```

3. **Check file format**:
   ```bash
   file sample.bam
   ```

4. **Re-index BAM file**:
   ```bash
   samtools index sample.bam
   ```

### Reference Sequence Issues

**Problem**: Reference FASTA file cannot be read

**Solutions**:
1. **Check FASTA file**:
   ```bash
   head -20 reference.fa
   samtools faidx reference.fa
   ```

2. **Validate format**:
   ```bash
   file reference.fa
   ```

3. **Check sequence names**:
   ```bash
   grep ">" reference.fa
   ```

### Low Coverage

**Problem**: Insufficient read depth for analysis

**Solutions**:
1. **Check coverage**:
   ```bash
   samtools depth sample.bam | awk '{sum+=$3} END {print sum/NR}'
   ```

2. **Filter low coverage**:
   ```bash
   modtector plot -M mod.csv -U unmod.csv -o plots/ -d 10
   ```

3. **Increase sequencing depth**:
   - Re-sequence with higher depth
   - Use more sensitive parameters

### Alignment Quality Issues

**Problem**: Poor alignment quality affecting results

**Solutions**:
1. **Check alignment quality**:
   ```bash
   samtools view sample.bam | head -100
   ```

2. **Filter alignments**:
   ```bash
   samtools view -q 20 -b sample.bam > sample_filtered.bam
   ```

3. **Check mapping rate**:
   ```bash
   samtools flagstat sample.bam
   ```

## Performance Issues

### Slow Processing

**Problem**: modtector runs very slowly

**Solutions**:
1. **Increase threads**:
   ```bash
   modtector count -b sample.bam -f reference.fa -o output.csv -t 8
   ```

2. **Use windowing**:
   ```bash
   modtector count -b sample.bam -f reference.fa -o output.csv -w 1000
   ```

3. **Optimize system**:
   - Use SSD storage
   - Increase RAM
   - Close other applications

### High CPU Usage

**Problem**: modtector uses too much CPU

**Solutions**:
1. **Reduce threads**:
   ```bash
   modtector count -b sample.bam -f reference.fa -o output.csv -t 2
   ```

2. **Monitor usage**:
   ```bash
   top -p $(pgrep modtector)
   ```

3. **Use nice command**:
   ```bash
   nice -n 10 modtector count -b sample.bam -f reference.fa -o output.csv
   ```

### Disk I/O Issues

**Problem**: High disk usage or slow I/O

**Solutions**:
1. **Use faster storage**:
   - Move data to SSD
   - Use local storage instead of network

2. **Optimize file operations**:
   ```bash
   modtector count -b sample.bam -f reference.fa -o /tmp/output.csv
   ```

3. **Monitor disk usage**:
   ```bash
   iostat -x 1
   ```

## Output Quality Issues

### Poor Normalization

**Problem**: Normalization results look incorrect

**Solutions**:
1. **Try different methods**:
   ```bash
   modtector norm -i input.csv -o output_percentile.csv -m percentile28
   modtector norm -i input.csv -o output_winsor.csv -m winsor90
   modtector norm -i input.csv -o output_boxplot.csv -m boxplot
   ```

2. **Adjust parameters**:
   ```bash
   modtector norm -i input.csv -o output.csv -m winsor90 --bases AC
   ```

3. **Check input data**:
   ```bash
   head -100 input.csv
   ```

### Low Reactivity Scores

**Problem**: Reactivity scores are too low or inconsistent

**Solutions**:
1. **Check data quality**:
   ```bash
   modtector plot -M mod.csv -U unmod.csv -o plots/ -c 0.1
   ```

2. **Try different methods**:
   ```bash
   modtector reactivity -M mod.csv -U unmod.csv -O output.csv -s ding
   modtector reactivity -M mod.csv -U unmod.csv -O output.csv -s rouskin
   ```

3. **Adjust parameters**:
   ```bash
   modtector reactivity -M mod.csv -U unmod.csv -O output.csv --pseudocount 0.5
   ```

### Poor Evaluation Scores

**Problem**: AUC, F1-score, or other metrics are low

**Solutions**:
1. **Check structure file**:
   ```bash
   head -20 structure.dp
   ```

2. **Try different signal types**:
   ```bash
   modtector evaluate -r reactivity.csv -s structure.dp -o results/ -g gene_id -t stop
   modtector evaluate -r reactivity.csv -s structure.dp -o results/ -g gene_id -t mutation
   modtector evaluate -r reactivity.csv -s structure.dp -o results/ -g gene_id -t both
   ```

3. **Disable auto-shift**:
   ```bash
   modtector evaluate -r reactivity.csv -s structure.dp -o results/ -g gene_id --no-auto-shift
   ```

4. **Check base matching**:
   ```bash
   modtector evaluate -r reactivity.csv -s structure.dp -o results/ -g gene_id --no-base-matching
   ```

## Memory and Resource Issues

### Out of Memory

**Problem**: modtector runs out of memory

**Solutions**:
1. **Reduce threads**:
   ```bash
   modtector count -b sample.bam -f reference.fa -o output.csv -t 1
   ```

2. **Use windowing**:
   ```bash
   modtector count -b sample.bam -f reference.fa -o output.csv -w 500
   ```

3. **Increase swap space**:
   ```bash
   sudo swapon -s
   sudo fallocate -l 8G /swapfile
   sudo chmod 600 /swapfile
   sudo mkswap /swapfile
   sudo swapon /swapfile
   ```

4. **Monitor memory usage**:
   ```bash
   free -h
   top -p $(pgrep modtector)
   ```

### Disk Space Issues

**Problem**: Not enough disk space for output files

**Solutions**:
1. **Check disk space**:
   ```bash
   df -h
   du -sh *
   ```

2. **Clean up old files**:
   ```bash
   rm -f *.tmp
   rm -f old_*.csv
   ```

3. **Use compression**:
   ```bash
   gzip *.csv
   ```

4. **Move to different location**:
   ```bash
   modtector count -b sample.bam -f reference.fa -o /path/to/large/disk/output.csv
   ```

## File Format Issues

### CSV Format Problems

**Problem**: CSV files cannot be read or have incorrect format

**Solutions**:
1. **Check CSV format**:
   ```bash
   head -5 file.csv
   file file.csv
   ```

2. **Validate CSV**:
   ```bash
   python -c "import pandas as pd; pd.read_csv('file.csv')"
   ```

3. **Fix encoding**:
   ```bash
   iconv -f utf-8 -t utf-8 file.csv > file_fixed.csv
   ```

### Structure File Issues

**Problem**: Secondary structure file cannot be read

**Solutions**:
1. **Check structure format**:
   ```bash
   head -20 structure.dp
   ```

2. **Validate dot-bracket notation**:
   - Ensure proper parentheses matching
   - Check for invalid characters
   - Verify sequence length matches

3. **Convert format if needed**:
   ```bash
   # Convert from other formats to dot-bracket
   # (Use appropriate conversion tool)
   ```

### BAM Index Issues

**Problem**: BAM file index is missing or corrupted

**Solutions**:
1. **Check for index**:
   ```bash
   ls -la sample.bam*
   ```

2. **Create index**:
   ```bash
   samtools index sample.bam
   ```

3. **Recreate index**:
   ```bash
   rm sample.bam.bai
   samtools index sample.bam
   ```

## Getting Help

### Log Files

**Check log files for detailed error information**:
```bash
# Check if log file was created
ls -la *.log

# View log content
cat modtector.log

# Search for errors
grep -i error modtector.log
grep -i warning modtector.log
```

### Debug Mode

**Enable verbose output for debugging**:
```bash
# Use log file
modtector count -b sample.bam -f reference.fa -o output.csv -l debug.log

# Check system resources
htop
iostat -x 1
```

### Common Error Messages

#### "BAM file not found"
- Check file path and permissions
- Verify file exists: `ls -la sample.bam`

#### "Reference sequence not found"
- Check FASTA file path
- Verify file format: `head -5 reference.fa`

#### "Insufficient memory"
- Reduce thread count
- Use windowing
- Increase system memory

#### "Low coverage"
- Check sequencing depth
- Adjust coverage thresholds
- Filter low-quality positions

#### "Poor alignment quality"
- Check BAM file quality
- Filter alignments by quality
- Re-align if necessary

### Reporting Issues

When reporting issues, include:

1. **System information**:
   ```bash
   uname -a
   rustc --version
   cargo --version
   ```

2. **Command used**:
   ```bash
   modtector count -b sample.bam -f reference.fa -o output.csv -t 4
   ```

3. **Error messages**:
   ```bash
   cat modtector.log
   ```

4. **Data information**:
   - BAM file size and type
   - Reference sequence information
   - System resources (RAM, CPU, disk)

5. **Expected vs actual behavior**:
   - What you expected to happen
   - What actually happened
   - Any error messages or warnings

### Community Support

- **GitHub Issues**: Report bugs and request features
- **Documentation**: Check this guide and other docs
- **Examples**: Look at example workflows
- **Forums**: Join community discussions

### Professional Support

For critical issues or commercial support:
- Contact the development team
- Consider professional bioinformatics services
- Check for commercial support options

## Prevention Tips

### Best Practices

1. **Validate inputs**:
   - Check BAM file quality
   - Verify reference sequences
   - Ensure proper file formats

2. **Monitor resources**:
   - Check available memory
   - Monitor disk space
   - Watch CPU usage

3. **Use appropriate parameters**:
   - Match thread count to system
   - Choose suitable normalization methods
   - Set appropriate thresholds

4. **Keep backups**:
   - Save intermediate results
   - Document parameters used
   - Version control your workflows

5. **Test with small datasets**:
   - Start with subset of data
   - Verify pipeline works
   - Scale up gradually

### Regular Maintenance

1. **Update software**:
   ```bash
   rustup update
   cargo update
   ```

2. **Clean up files**:
   ```bash
   cargo clean
   rm -f *.tmp
   ```

3. **Check system health**:
   ```bash
   df -h
   free -h
   ```

This troubleshooting guide should help you resolve most common issues. If you encounter problems not covered here, please refer to the [Getting Help](#getting-help) section for additional support options.
