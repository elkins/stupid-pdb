import logging
import openmm.app as app
import openmm as mm
from openmm import unit
import sys
import os

logger = logging.getLogger(__name__)

class DockingPrep:
    """
    Utilities for preparing PDB structures for molecular docking.
    """
    
    def __init__(self, forcefield_name='amber14-all.xml'):
        self.forcefield_name = forcefield_name
        try:
            self.forcefield = app.ForceField(self.forcefield_name, 'amber14/tip3pfb.xml')
        except Exception as e:
            logger.error(f"Failed to load forcefield '{forcefield_name}': {e}")
            raise
            
    def write_pqr(self, input_pdb: str, output_pqr: str) -> bool:
        """
        Converts a PDB file to PQR format (adding partial charges and radii).
        
        Uses OpenMM to assign charges based on the selected forcefield.
        Radii are derived from Lennard-Jones sigma (sigma / 2).
        
        Args:
            input_pdb: Path to input PDB.
            output_pqr: Path to output PQR.
            
        Returns:
            bool: True if successful.
        """
        try:
            # 1. Load PDB
            pdb = app.PDBFile(input_pdb)
            
            # 2. Add Hydrogens (Crucial for correct charge assignment)
            modeller = app.Modeller(pdb.topology, pdb.positions)
            modeller.addHydrogens(self.forcefield, pH=7.4) # Physiological pH
            
            # 3. Create System to get forces (charges/sigmas)
            system = self.forcefield.createSystem(
                modeller.topology,
                nonbondedMethod=app.NoCutoff,
                constraints=None,
                rigidWater=False
            )
            
            # 4. Extract NonbondedForce
            nonbonded = None
            for force in system.getForces():
                if isinstance(force, mm.NonbondedForce):
                    nonbonded = force
                    break
            
            if not nonbonded:
                raise ValueError("Forcefield does not contain NonbondedForce (no charges found).")
            
            # 5. Write PQR
            # Standard PDB columns, but Occupancy -> Charge, B-factor -> Radius
            # We can use Biotite or manual writing. Manual gives full control over PQR format quirks.
            # But converting OpenMM Topology -> Biotite AtomArray is verbose.
            # Let's write manually by iterating OpenMM topology.
            
            with open(output_pqr, 'w') as f:
                atom_idx = 0
                for chain in modeller.topology.chains():
                    for residue in chain.residues():
                        for atom in residue.atoms():
                            # Get parameters (Charge, Sigma, Epsilon)
                            charge, sigma, epsilon = nonbonded.getParticleParameters(atom.index)
                            
                            # Convert units
                            # charge is in elementary charge (float)
                            # sigma is in nanometers (Quantity)
                            # radius = sigma / 2
                            
                            q = charge.value_in_unit(unit.elementary_charge)
                            r_nm = sigma.value_in_unit(unit.nanometer) * 0.5
                            r_angstrom = r_nm * 10.0
                            
                            pos = modeller.positions[atom.index]
                            x = pos[0].value_in_unit(unit.angstrom)
                            y = pos[1].value_in_unit(unit.angstrom)
                            z = pos[2].value_in_unit(unit.angstrom)
                            
                            atom_idx += 1
                            
                            # PQR Format (Standard)
                            # Field_name Atom_number Atom_name Residue_name Chain_ID Residue_number X Y Z Charge Radius
                            # Note: Chain ID is often optional or skipped in some PQR variants (like PDB2PQR whitespace separated)
                            # We'll stick to fixed width PDB-like if possible, or whitespace
                            # ATOM      1  N   ALA A   1      27.525  26.046  14.628  0.1340 1.8240
                            
                            line = f"ATOM  {atom_idx:>5} {atom.name:<4} {residue.name:<3} {chain.id:<1} {int(residue.id):>4}    {x:>8.3f}{y:>8.3f}{z:>8.3f} {q:>6.4f} {r_angstrom:>6.4f}\n"
                            f.write(line)
                            
            logger.info(f"Successfully wrote PQR to {output_pqr}")
            return True
            
        except Exception as e:
            logger.error(f"PQR generation failed: {e}")
            return False
