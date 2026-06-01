# Contributing to Apple MAS Toolkit

Thank you for your interest in contributing to the Apple MAS Toolkit! This guide will help you get started.

## How to Contribute

### Reporting Issues

- Use the GitHub Issues tab to report bugs
- Include a clear description of the issue
- Provide steps to reproduce the problem
- Include your Python version and operating system

### Suggesting Features

- Open a GitHub Issue with the "enhancement" label
- Describe the feature and its use case
- Explain how it would benefit apple breeding research

### Submitting Changes

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/my-feature`)
3. Make your changes
4. Add tests for new functionality
5. Ensure all tests pass (`pytest tests/ -v`)
6. Commit your changes (`git commit -m "Add my feature"`)
7. Push to your branch (`git push origin feature/my-feature`)
8. Open a Pull Request

## Development Setup

```bash
# Clone the repository
git clone https://github.com/Mosiuropu/apple-mas-toolkit.git
cd apple-mas-toolkit

# Create a virtual environment
python -m venv venv
source venv/bin/activate  # Linux/Mac
# or
venv\Scripts\activate  # Windows

# Install development dependencies
pip install -e ".[dev,notebooks]"
```

## Code Style

- Follow PEP 8 guidelines
- Use type hints where appropriate
- Write docstrings for all public functions and classes
- Keep functions focused and under 50 lines when possible

## Testing

```bash
# Run all tests
pytest tests/ -v

# Run with coverage
pytest tests/ --cov=apple_mas --cov-report=html

# Run specific test file
pytest tests/test_all.py -v
```

## Adding New Markers

To add a new marker to the database:

1. Edit `data/marker_database/apple_markers.json`
2. Add the marker entry with all required fields
3. Update `trait_groups` and `chromosomal_locations`
4. Add tests for the new marker
5. Update documentation in `docs/marker_reference.md`

## Documentation

- Update README.md if adding new features
- Add docstrings to all new functions
- Update `docs/marker_reference.md` for new markers
- Include examples in docstrings

## Pull Request Guidelines

- Keep PRs focused on a single change
- Include tests for new functionality
- Update documentation as needed
- Follow the existing code style
- Write a clear PR description

## Questions?

If you have questions about contributing, feel free to open an Issue with the "question" label.
