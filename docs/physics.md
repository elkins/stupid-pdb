# Biophysics 101

## Understanding Energy Minimization

**Energy Minimization** is the process of moving atoms "downhill" to find the nearest stable shape.

```text
      High Energy
      (Unstable)
          |
         / \       Forces push atoms "downhill"
        /   \     (Gradient Descent)
       /     \
      /       \___
     /            \
    /              \__ Low Energy
   /                  (Stable / Minimized)
```

`synth-pdb` uses **Implicit Solvent (OBC2)** to simulate the effect of water without the computational cost of thousands of explicit water molecules.

## The Generation Pipeline

```text
[User] -> [Generator] -> [Geometry Builder] -> [Sidechain Packer] -> [Energy Minimizer] -> [PDB File]
             ^                  |                    |                      |
             |              (N-CA-C-O)           (Rotamers)             (OpenMM)
             |                                       |                      |
             +---------------------------------------+----------------------+
```
