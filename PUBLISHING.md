# Publishing to PyPI

This guide explains how to publish `synth-pdb` to PyPI for easier installation and better discoverability.

## Prerequisites

1. Create accounts:
   - PyPI: https://pypi.org/account/register/
   - TestPyPI (for testing): https://test.pypi.org/account/register/

2. Install build tools:
```bash
pip install build twine
```

## Step 1: Build the Package

```bash
# Clean previous builds
rm -rf dist/ build/ *.egg-info

# Build the package
python -m build
```

This creates:
- `dist/synth_pdb-1.0.0-py3-none-any.whl` (wheel)
- `dist/synth-pdb-1.0.0.tar.gz` (source distribution)

## Step 2: Test on TestPyPI (Recommended)

```bash
# Upload to TestPyPI
python -m twine upload --repository testpypi dist/*

# Test installation
pip install --index-url https://test.pypi.org/simple/ synth-pdb
```

## Step 3: Upload to PyPI

```bash
# Upload to PyPI
python -m twine upload dist/*
```

You'll be prompted for your PyPI username and password.

## Step 4: Verify Installation

```bash
# Install from PyPI
pip install synth-pdb

# Test it works
synth-pdb --help
```

## Using API Tokens (Recommended)

Instead of username/password, use API tokens:

1. Go to https://pypi.org/manage/account/token/
2. Create a new API token
3. Create `~/.pypirc`:

```ini
[pypi]
username = __token__
password = pypi-YOUR-TOKEN-HERE

[testpypi]
username = __token__
password = pypi-YOUR-TESTPYPI-TOKEN-HERE
```

## Automation with GitHub Actions

Create `.github/workflows/publish.yml`:

```yaml
name: Publish to PyPI

on:
  release:
    types: [published]

jobs:
  deploy:
    runs-on: ubuntu-latest
    steps:
    - uses: actions/checkout@v3
    - name: Set up Python
      uses: actions/setup-python@v4
      with:
        python-version: '3.x'
    - name: Install dependencies
      run: |
        python -m pip install --upgrade pip
        pip install build twine
    - name: Build package
      run: python -m build
    - name: Publish to PyPI
      env:
        TWINE_USERNAME: __token__
        TWINE_PASSWORD: ${{ secrets.PYPI_API_TOKEN }}
      run: twine upload dist/*
```

Add your PyPI API token as a GitHub secret named `PYPI_API_TOKEN`.

## Version Management

Update version in `pyproject.toml` before each release:

```toml
version = "1.0.1"  # Increment for each release
```

Follow semantic versioning:
- `1.0.0` → `1.0.1` for bug fixes
- `1.0.0` → `1.1.0` for new features
- `1.0.0` → `2.0.0` for breaking changes

## After Publishing

1. Verify on PyPI: https://pypi.org/project/synth-pdb/
2. Test installation: `pip install synth-pdb`
3. Update README with PyPI badge:
   ```markdown
   [![PyPI version](https://badge.fury.io/py/synth-pdb.svg)](https://badge.fury.io/py/synth-pdb)
   ```

## Troubleshooting

- **"File already exists"**: You can't re-upload the same version. Increment version number.
- **"Invalid distribution"**: Check `pyproject.toml` syntax
- **"README rendering failed"**: Validate your README.md markdown
