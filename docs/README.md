# modtector Documentation

This directory contains the documentation for modtector, a high-performance RNA modification detection tool.

## Building the Documentation

### Prerequisites

- Python 3.7 or higher
- pip

### Installation

1. **Install dependencies**:
   ```bash
   pip install -r requirements.txt
   ```

2. **Build the documentation**:
   ```bash
   make html
   ```

3. **View the documentation**:
   ```bash
   make serve
   ```
   Then open http://localhost:8000 in your browser.

### Development

For development, you can use:

```bash
# Build with warnings
make dev

# Watch for changes and rebuild
make watch

# Check for errors
make check

# Clean build directory
make clean
```

## Documentation Structure

```
docs/
├── index.md              # Main documentation page
├── installation.md       # Installation guide
├── quickstart.md         # Quick start guide
├── user-guide.md         # Comprehensive user guide
├── commands.md           # Command reference
├── examples.md           # Examples and tutorials
├── troubleshooting.md    # Troubleshooting guide
├── contributing.md       # Contributing guidelines
├── conf.py              # Sphinx configuration
├── requirements.txt     # Python dependencies
├── Makefile            # Build automation
├── _static/            # Static assets
│   └── custom.css      # Custom styling
└── README.md           # This file
```

## Writing Documentation

### Markdown Format

The documentation is written in Markdown with MyST extensions:

- Use `#` for headers
- Use ````bash` for code blocks
- Use `**bold**` and `*italic*` for emphasis
- Use `[link text](url)` for links
- Use `![alt text](image.png)` for images

### Code Examples

Use fenced code blocks with language specification:

````markdown
```bash
modtector count -b sample.bam -f reference.fa -o output.csv
```
````

### Tables

Use Markdown tables for structured data:

```markdown
| Option | Description | Default |
|--------|-------------|---------|
| `-b` | BAM file path | Required |
| `-f` | FASTA file path | Required |
| `-o` | Output file path | Required |
```

### Admonitions

Use MyST admonitions for important information:

````markdown
```{note}
This is a note with important information.
```

```{warning}
This is a warning about potential issues.
```

```{tip}
This is a tip for better usage.
```
````

### Cross-References

Link to other sections:

```markdown
See the [Installation Guide](installation.md) for setup instructions.
```

## Customization

### Styling

Custom CSS is in `_static/custom.css`. Key features:

- Custom color scheme
- Improved code block styling
- Better table formatting
- Responsive design
- Print styles

### Configuration

Sphinx configuration is in `conf.py`. Key settings:

- Project information
- Extensions
- Theme settings
- Custom options

## Deployment

### ReadTheDocs

For ReadTheDocs deployment:

1. Connect your repository to ReadTheDocs
2. Set the documentation type to "Sphinx"
3. Set the configuration file to `docs/conf.py`
4. Set the requirements file to `docs/requirements.txt`

### Local Deployment

For local deployment:

```bash
# Build for production
make deploy

# The built documentation will be in _build/html/
```

## Contributing

When contributing to the documentation:

1. Follow the existing style and structure
2. Test your changes locally
3. Ensure all links work
4. Check for typos and grammar
5. Update the table of contents if needed

### Guidelines

- Use clear, concise language
- Provide examples where helpful
- Include error conditions and solutions
- Keep formatting consistent
- Test all code examples

## Tools and Extensions

### Sphinx Extensions

- **MyST Parser**: Markdown support
- **sphinx-copybutton**: Copy code blocks
- **sphinx-tabs**: Tabbed content
- **sphinx-rtd-theme**: ReadTheDocs theme

### Development Tools

- **pre-commit**: Code quality hooks
- **black**: Code formatting
- **flake8**: Linting

## Troubleshooting

### Common Issues

1. **Build errors**: Check Python version and dependencies
2. **Missing extensions**: Install required packages
3. **Theme issues**: Check CSS and configuration
4. **Link errors**: Verify file paths and references

### Getting Help

- Check the Sphinx documentation
- Review the MyST parser documentation
- Look at the ReadTheDocs theme documentation
- Ask for help in project discussions

## License

The documentation is licensed under the same license as the modtector project.
