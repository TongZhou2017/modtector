# Contributing to modtector

Thank you for your interest in contributing to modtector! This guide will help you get started with contributing to the project.

## Table of Contents

1. [Getting Started](#getting-started)
2. [Development Setup](#development-setup)
3. [Contributing Guidelines](#contributing-guidelines)
4. [Code Style](#code-style)
5. [Testing](#testing)
6. [Documentation](#documentation)
7. [Submitting Changes](#submitting-changes)
8. [Issue Reporting](#issue-reporting)

## Getting Started

### Prerequisites

- Rust 1.70 or higher
- Cargo (included with Rust)
- Git
- Basic knowledge of bioinformatics and RNA modifications

### Fork and Clone

1. **Fork the repository** on [GitHub](https://github.com/TongZhou2017/modtector)
2. **Clone your fork**:
   ```bash
   git clone https://github.com/your-username/modtector.git
   cd modtector
   ```

3. **Add upstream remote**:
   ```bash
   git remote add upstream https://github.com/TongZhou2017/modtector.git
   ```

### Development Setup

1. **Install dependencies**:
   ```bash
   cargo build
   ```

2. **Run tests**:
   ```bash
   cargo test
   ```

3. **Check code style**:
   ```bash
   cargo fmt
   cargo clippy
   ```

## Contributing Guidelines

### Types of Contributions

We welcome various types of contributions:

1. **Bug fixes**: Fix existing issues
2. **New features**: Add new functionality
3. **Documentation**: Improve documentation
4. **Tests**: Add or improve tests
5. **Performance**: Optimize existing code
6. **Examples**: Add usage examples

### Contribution Process

1. **Check existing issues**: Look for existing issues or discussions
2. **Create an issue**: If no existing issue, create one to discuss
3. **Fork and branch**: Create a feature branch
4. **Make changes**: Implement your changes
5. **Test**: Ensure all tests pass
6. **Document**: Update documentation if needed
7. **Submit PR**: Create a pull request

### Branch Naming

Use descriptive branch names:
- `fix/issue-123-description`
- `feature/new-normalization-method`
- `docs/update-user-guide`
- `test/add-integration-tests`

## Code Style

### Rust Style Guidelines

Follow Rust community standards:

1. **Formatting**: Use `cargo fmt`
2. **Linting**: Use `cargo clippy`
3. **Naming**: Use snake_case for variables and functions
4. **Documentation**: Document public functions and structs

### Code Formatting

```bash
# Format code
cargo fmt

# Check formatting
cargo fmt -- --check
```

### Linting

```bash
# Run clippy
cargo clippy

# Run clippy with all warnings
cargo clippy -- -W clippy::all
```

### Documentation

Document all public functions and structs:

```rust
/// Calculate reactivity scores between modified and unmodified samples.
///
/// # Arguments
///
/// * `mod_data` - Modified sample data
/// * `unmod_data` - Unmodified sample data
/// * `method` - Reactivity calculation method
///
/// # Returns
///
/// * `Result<Vec<f64>, Error>` - Calculated reactivity scores
///
/// # Examples
///
/// ```
/// let reactivity = calculate_reactivity(mod_data, unmod_data, "current")?;
/// ```
pub fn calculate_reactivity(
    mod_data: &[f64],
    unmod_data: &[f64],
    method: &str,
) -> Result<Vec<f64>, Box<dyn Error>> {
    // Implementation
}
```

## Testing

### Test Structure

Organize tests by module:

```rust
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_calculate_reactivity() {
        // Test implementation
    }

    #[test]
    fn test_normalize_signals() {
        // Test implementation
    }
}
```

### Running Tests

```bash
# Run all tests
cargo test

# Run specific test
cargo test test_calculate_reactivity

# Run tests with output
cargo test -- --nocapture

# Run integration tests
cargo test --test integration_tests
```

### Test Data

Use small, representative test data:

```rust
#[test]
fn test_pileup_analysis() {
    let test_bam = "tests/data/sample.bam";
    let test_fasta = "tests/data/reference.fa";
    
    // Test implementation
}
```

### Integration Tests

Create integration tests in `tests/` directory:

```rust
// tests/integration_tests.rs
use modtector;

#[test]
fn test_full_workflow() {
    // Test complete workflow
}
```

## Documentation

### Code Documentation

Document all public APIs:

```rust
/// Normalize signals using specified method.
///
/// # Arguments
///
/// * `signals` - Input signal data
/// * `method` - Normalization method
///
/// # Returns
///
/// * `Result<Vec<f64>, Error>` - Normalized signals
pub fn normalize_signals(signals: &[f64], method: &str) -> Result<Vec<f64>, Box<dyn Error>> {
    // Implementation
}
```

### User Documentation

Update user-facing documentation:

1. **README.md**: Update main documentation
2. **docs/**: Update user guides and examples
3. **CHANGELOG.md**: Document changes
4. **API documentation**: Update API docs

### Documentation Standards

- Use clear, concise language
- Provide examples where helpful
- Include error conditions
- Document parameter ranges and defaults

## Submitting Changes

### Pull Request Process

1. **Create feature branch**:
   ```bash
   git checkout -b feature/new-feature
   ```

2. **Make changes**:
   - Implement your changes
   - Add tests
   - Update documentation

3. **Commit changes**:
   ```bash
   git add .
   git commit -m "Add new normalization method"
   ```

4. **Push to fork**:
   ```bash
   git push origin feature/new-feature
   ```

5. **Create pull request**:
   - Go to GitHub
   - Click "New Pull Request"
   - Fill out the template

### Pull Request Template

```markdown
## Description
Brief description of changes

## Type of Change
- [ ] Bug fix
- [ ] New feature
- [ ] Documentation update
- [ ] Performance improvement
- [ ] Test addition

## Testing
- [ ] All tests pass
- [ ] New tests added
- [ ] Manual testing completed

## Documentation
- [ ] Code documented
- [ ] User documentation updated
- [ ] API documentation updated

## Checklist
- [ ] Code follows style guidelines
- [ ] Self-review completed
- [ ] No breaking changes
- [ ] Backward compatibility maintained
```

### Review Process

1. **Automated checks**: CI/CD pipeline runs
2. **Code review**: Maintainers review code
3. **Testing**: Verify functionality
4. **Documentation**: Check documentation updates
5. **Approval**: Maintainer approves
6. **Merge**: Changes merged to main branch

## Issue Reporting

### Bug Reports

When reporting bugs, include:

1. **Environment**:
   - Operating system
   - Rust version
   - modtector version

2. **Reproduction steps**:
   - Exact command used
   - Input data description
   - Expected vs actual behavior

3. **Error messages**:
   - Full error output
   - Log files
   - Stack traces

4. **Additional context**:
   - Related issues
   - Workarounds tried
   - Impact assessment

### Feature Requests

When requesting features, include:

1. **Use case**: Why is this feature needed?
2. **Proposed solution**: How should it work?
3. **Alternatives**: What alternatives were considered?
4. **Implementation**: Any implementation ideas?

### Issue Template

```markdown
## Bug Report / Feature Request

### Description
Brief description of the issue or feature request

### Environment
- OS: [e.g., Ubuntu 20.04]
- Rust version: [e.g., 1.70.0]
- modtector version: [e.g., v0.9.2]

### Steps to Reproduce (for bugs)
1. Step 1
2. Step 2
3. Step 3

### Expected Behavior
What you expected to happen

### Actual Behavior
What actually happened

### Additional Context
Any other relevant information
```

## Development Workflow

### Daily Workflow

1. **Sync with upstream**:
   ```bash
   git fetch upstream
   git checkout main
   git merge upstream/main
   ```

2. **Create feature branch**:
   ```bash
   git checkout -b feature/your-feature
   ```

3. **Make changes**:
   - Write code
   - Add tests
   - Update documentation

4. **Test changes**:
   ```bash
   cargo test
   cargo clippy
   cargo fmt
   ```

5. **Commit and push**:
   ```bash
   git add .
   git commit -m "Descriptive commit message"
   git push origin feature/your-feature
   ```

### Release Process

1. **Version bump**: Update version in Cargo.toml
2. **Changelog**: Update CHANGELOG.md
3. **Tag release**: Create git tag
4. **Build release**: Create release build
5. **Publish**: Publish to crates.io (if applicable)

## Code of Conduct

### Our Pledge

We are committed to providing a welcoming and inclusive environment for all contributors.

### Expected Behavior

- Be respectful and inclusive
- Accept constructive criticism
- Focus on what's best for the community
- Show empathy towards others

### Unacceptable Behavior

- Harassment or discrimination
- Trolling or inflammatory comments
- Personal attacks
- Inappropriate language or imagery

### Enforcement

Violations will be addressed by project maintainers and may result in:
- Warning
- Temporary ban
- Permanent ban

## Getting Help

### Communication Channels

1. **GitHub Issues**: For bugs and feature requests
2. **GitHub Discussions**: For questions and discussions
3. **Email**: For sensitive issues
4. **Documentation**: Check existing docs first

### Mentorship

New contributors can:
- Ask for help in GitHub Discussions
- Request code review
- Ask for guidance on implementation
- Seek clarification on requirements

### Resources

- [Rust Book](https://doc.rust-lang.org/book/)
- [Rust API Guidelines](https://rust-lang.github.io/api-guidelines/)
- [Bioinformatics Best Practices](https://www.bioconductor.org/developers/how-to/best-practices-for-package-development/)

## Recognition

Contributors will be recognized in:
- CONTRIBUTORS.md file
- Release notes
- Project documentation
- Community acknowledgments

Thank you for contributing to modtector! Your contributions help make RNA modification detection more accessible and reliable for the scientific community.
