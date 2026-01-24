# synth-pdb

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/elkins/synth-pdb/blob/main/demo.ipynb)
[![PyPI version](https://badge.fury.io/py/synth-pdb.svg)](https://badge.fury.io/py/synth-pdb)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/latestdoi/123456789)](https://zenodo.org/badge/latestdoi/123456789)

**Generate realistic PDB files with mixed secondary structures for bioinformatics testing, education, and tool development.**

> ⚠️ **Important**: The generated structures use idealized geometries and may contain violations of standard structural constraints. These files are intended for **testing computational tools** and **educational demonstrations**, not for simulation or experimental validation.

---

## Why synth-pdb?

In the fields of structural biology and bioinformatics, researchers frequently require datasets of protein structures to test algorithms, train machine learning models, or validatate analytical pipelines. While the Protein Data Bank (PDB) contains over 200,000 experimental structures, relying solely on experimental data has limitations:

1.  **Bias**: PDB data is biased toward crystallizable or stable proteins.
2.  **Complexity**: Experimental files often contain artifacts, missing atoms, or non-standard residues that complicate initial testing.
3.  **Lack of Ground Truth**: When developing algorithms for NMR assignment or structure calculation, "perfect" synthetic data is essential for unit testing.

``synth-pdb`` fills this gap by providing a lightweight, deterministic generator that produces chemically valid, full-atom PDB files with user-defined secondary structures (helices, sheets) in seconds.

## Key Features

✨ **Structure Generation**
- Full atomic representation with backbone and side-chain heavy atoms + hydrogens
- Customizable sequence (1-letter or 3-letter amino acid codes)
- **Conformational diversity**: Generate alpha helices, beta sheets, extended chains, or random conformations
- **Rotamer-based side-chain placement** for all 20 standard amino acids (Dunbrack library)

🔬 **Validation Suite**
- Bond length & angle validation
- Ramachandran angle checking
- Steric clash detection & refinement
- Sequence improbability detection

⚙️ **Quality Control**
- `--best-of-N`: Generate multiple structures and select the one with fewest violations
- `--minimize`: Relax structures using OpenMM (Implicit Solvent / AMBER forcefield)

## Quick Visual Demo

Want to see the **Physics + Visualization** capabilities in action?

Run this command to generate a **Leucine Zipper** (classic alpha helix), **minimize** its energy using OpenMM, and immediately **visualize** it in your browser:

```bash
synth-pdb --sequence "LKELEKELEKELEKELEKELEKEL" --conformation alpha --minimize --visualize
```

## Citation

If you use `synth-pdb` in your research, please cite it using the metadata in the `CITATION.cff` file in the repository.
