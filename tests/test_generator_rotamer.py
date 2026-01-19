import unittest
import logging
import numpy as np
from typing import List, Dict, Any

from stupid_pdb.generator import generate_pdb_content

# Suppress logging during tests to keep output clean
logging.getLogger().setLevel(logging.CRITICAL)

def _parse_pdb_atoms(pdb_content: str) -> List[Dict[str, Any]]:
    """
    Parses the PDB content and extracts atom information, specifically coordinates.
    Returns a list of dictionaries, each representing an atom with residue and chain info.
    """
    parsed_atoms = []
    for line in pdb_content.splitlines():
        stripped_line = line.strip()
        if stripped_line.startswith("ATOM") or stripped_line.startswith("HETATM"):
            try:
                atom_number = int(stripped_line[6:11].strip())
                atom_name = stripped_line[12:16].strip()
                alt_loc = stripped_line[16].strip()
                residue_name = stripped_line[17:20].strip()
                chain_id = stripped_line[21].strip()
                residue_number = int(stripped_line[22:26].strip())
                insertion_code = stripped_line[26].strip()
                x = float(stripped_line[30:38])
                y = float(stripped_line[38:46])
                z = float(stripped_line[46:54])
                occupancy = float(stripped_line[54:60].strip())
                temp_factor = float(stripped_line[60:66].strip())
                element = stripped_line[76:78].strip()
                charge = stripped_line[78:80].strip()

                parsed_atoms.append(
                    {
                        "atom_number": atom_number,
                        "atom_name": atom_name,
                        "alt_loc": alt_loc,
                        "residue_name": residue_name,
                        "chain_id": chain_id,
                        "residue_number": residue_number,
                        "insertion_code": insertion_code,
                        "coords": np.array([x, y, z]),
                        "occupancy": occupancy,
                        "temp_factor": temp_factor,
                        "element": element,
                        "charge": charge,
                    }
                )
            except (ValueError, IndexError) as e:
                logging.warning(
                    f"Could not parse PDB ATOM/HETATM line: {line.strip()} - {e}"
                )
    return parsed_atoms

def _calculate_dihedral_angle(
    p1: np.ndarray, p2: np.ndarray, p3: np.ndarray, p4: np.ndarray
) -> float:
    """
    Calculates the dihedral angle (in degrees) defined by four points (p1, p2, p3, p4).
    """
    b1 = p2 - p1
    b2 = p3 - p2
    b3 = p4 - p3

    # Normalize b2 for better numerical stability
    b2_norm = np.linalg.norm(b2)
    if b2_norm < 1e-6:
        return 0.0
    b2 = b2 / b2_norm

    v = b1 - np.dot(b1, b2) * b2
    w = b3 - np.dot(b3, b2) * b2

    x = np.dot(v, w)
    y = np.dot(np.cross(b2, v), w)

    return np.degrees(np.arctan2(y, x))


class TestGeneratorWithRotamers(unittest.TestCase):

    def test_generate_pdb_with_rotamers_runs_without_error(self):
        """
        Test that PDB generation with side-chain rotamers completes without error.
        This is the initial failing test for the TDD process.
        """
        try:
            generate_pdb_content(sequence_str="ALA-LEU-ARG", full_atom=True, use_rotamers=True)
        except Exception as e:
            self.fail(f"generate_pdb_content with rotamers raised an exception: {e}")

    def test_rotamer_usage_changes_chi_angle(self):
        """
        Test that using rotamers results in a non-zero chi-1 angle for LEU,
        different from the hardcoded geometry.
        """
        # Generate a peptide with a LEU residue
        pdb_content = generate_pdb_content(sequence_str="LEU", full_atom=True, use_rotamers=True)
        atoms = _parse_pdb_atoms(pdb_content)

        # Get the atoms for the chi-1 angle of LEU
        leu_atoms = {atom['atom_name']: atom for atom in atoms if atom['residue_name'] == 'LEU'}
        
        self.assertIn('N', leu_atoms, "N atom not found for LEU")
        self.assertIn('CA', leu_atoms, "CA atom not found for LEU")
        self.assertIn('CB', leu_atoms, "CB atom not found for LEU")
        self.assertIn('CG', leu_atoms, "CG atom not found for LEU")

        # Get coordinates
        n_coords = leu_atoms['N']['coords']
        ca_coords = leu_atoms['CA']['coords']
        cb_coords = leu_atoms['CB']['coords']
        cg_coords = leu_atoms['CG']['coords']

        # Calculate the chi-1 angle
        chi1_angle = _calculate_dihedral_angle(n_coords, ca_coords, cb_coords, cg_coords)

        # The hardcoded geometry places atoms linearly, resulting in a dihedral angle
        # of ~0 or ~180 degrees. This test will fail if the angle is close to 0.
        self.assertNotAlmostEqual(chi1_angle, 0.0, places=1, msg="Chi-1 angle for LEU is close to 0.0, indicating hardcoded geometry is likely still in use.")
        self.assertNotAlmostEqual(abs(chi1_angle), 180.0, places=1, msg="Chi-1 angle for LEU is close to 180.0, indicating hardcoded geometry is likely still in use.")


if __name__ == '__main__':
    unittest.main()