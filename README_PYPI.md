# synth-pdb

**Generate realistic PDB files with mixed secondary structures for bioinformatics testing, education, and tool development.**

## Key Features

- 🧬 **Mixed Secondary Structures**: Create helix-turn-sheet motifs, zinc fingers, immunoglobulins
- 🔬 **Biologically Realistic**: Ramachandran-based sampling, rotamer libraries, chirality validation
- ✅ **Comprehensive Validation**: Bond geometry, steric clashes, peptide planarity
- 🎯 **Quality Control**: Best-of-N selection, iterative clash refinement
- 📚 **Educational**: Extensive code comments explaining biophysical principles

## Quick Start

```bash
# Install
pip install synth-pdb

# Generate a simple peptide
synth-pdb --length 20

# Create a helix-turn-helix motif
synth-pdb --length 25 --structure "1-10:alpha,11-15:random,16-25:alpha"

# Generate and validate
synth-pdb --sequence "ACDEFGHIKLMNPQRSTVWY" --validate
```

## Use Cases

- 🧪 Test bioinformatics tools and structure analysis software
- 📚 Teach protein structure concepts and secondary structure
- 🔬 Generate training data for machine learning models
- 🧬 Create test datasets for structure prediction algorithms
- 🎯 Benchmark validation and refinement protocols

## Documentation

See [README.md](README.md) for full documentation.

## Citation

If you use this software in your research, please cite:

```bibtex
@software{synth_pdb,
  author = {Elkins, George},
  title = {synth-pdb: Realistic Protein Structure Generator},
  year = {2026},
  url = {https://github.com/elkins/synth-pdb}
}
```

## License

MIT License - see [LICENSE](LICENSE) for details.
