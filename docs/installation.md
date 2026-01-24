# Installation

## From PyPI (Recommended)

Install the latest stable release from PyPI:

```bash
pip install synth-pdb
```

This installs the `synth-pdb` package and makes the `synth-pdb` command available system-wide.

## From Bioconda

Bioinformatics users can install via Conda:

```bash
conda install -c bioconda synth-pdb
```

## From Source

Install directly from the project directory for development:

```bash
git clone https://github.com/elkins/synth-pdb.git
cd synth-pdb
pip install .
```

## Requirements

- Python 3.8+
- NumPy
- Biotite
- OpenMM (optional, for energy minimization)
- 3Dmol.js (optional, for visualization)
