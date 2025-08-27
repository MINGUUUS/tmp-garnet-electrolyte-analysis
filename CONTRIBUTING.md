# Contributing to Garnet Electrolyte Analysis

We welcome contributions to the Garnet Electrolyte Analysis package! This document provides guidelines for contributing.

## Getting Started

1. Fork the repository on GitHub
2. Clone your fork locally
3. Create a development environment:
   ```bash
   conda env create -f environment.yml
   conda activate garnet-analysis
   pip install -e ".[dev]"
   ```

## Development Workflow

1. Create a new branch for your feature or bug fix:
   ```bash
   git checkout -b feature/your-feature-name
   ```

2. Make your changes and add tests if applicable

3. Run the test suite:
   ```bash
   pytest tests/
   ```

4. Check code formatting:
   ```bash
   black garnet_analysis/
   flake8 garnet_analysis/
   ```

5. Commit your changes with a descriptive commit message

6. Push to your fork and submit a pull request

## Code Style Guidelines

- Follow PEP 8 style guide
- Use type hints where possible
- Include docstrings for all public functions and classes
- Use meaningful variable and function names
- Keep functions focused and modular

## Documentation

- Update docstrings for any new or modified functions
- Add example usage in docstrings where helpful
- Update README.md if adding new features
- Add or update tutorial notebooks for major features

## Testing

- Add unit tests for new functionality
- Ensure existing tests continue to pass
- Aim for good test coverage of critical functions
- Use pytest fixtures for common test setup

## Submitting Changes

### Pull Request Guidelines

- Provide a clear description of the changes
- Reference any related issues
- Include tests for new functionality
- Update documentation as needed
- Ensure CI tests pass

### Commit Message Format

Use descriptive commit messages:
```
Add temperature-dependent conductivity calculation

- Implement Nernst-Einstein equation for conductivity
- Add unit tests for conversion functions
- Update example notebook with conductivity analysis
```

## Types of Contributions

### Bug Reports

When filing bug reports, please include:
- Python version and environment details
- Minimal code example that reproduces the issue
- Expected vs actual behavior
- Error messages and stack traces

### Feature Requests

For feature requests, please:
- Describe the use case and motivation
- Provide examples of how the feature would be used
- Consider backwards compatibility

### Code Contributions

Areas where contributions are especially welcome:
- Additional analysis methods
- Performance optimizations
- Better error handling and validation
- Documentation improvements
- Example notebooks and tutorials
- Support for additional file formats

## Code Organization

```
garnet_analysis/
├── core/           # Core analysis algorithms
├── io/            # File I/O functionality
├── plotting/      # Visualization functions
├── utils/         # Utilities and constants
└── config.py      # Configuration management
```

## Dependencies

- Keep dependencies minimal
- Use established scientific Python packages (numpy, scipy, matplotlib)
- Make heavy dependencies optional when possible
- Update requirements.txt and environment.yml as needed

## Release Process

Releases are managed by project maintainers:

1. Update version in `__init__.py`
2. Update CHANGELOG.md
3. Create release tag
4. Build and upload to PyPI

## Getting Help

- Open an issue for questions about development
- Reach out to maintainers for guidance on major changes
- Check existing issues and PRs for similar work

## Recognition

All contributors will be acknowledged in the project. Significant contributions may be recognized with co-authorship on related publications.

Thank you for contributing to advancing computational materials science!