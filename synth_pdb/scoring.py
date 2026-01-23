import numpy as np
import biotite.structure as struc
import logging
from .data import VAN_DER_WAALS_RADII

logger = logging.getLogger(__name__)

def calculate_clash_score(atom_array: struc.AtomArray) -> float:
    """
    Calculate a simple clash score based on Van der Waals overlaps.
    Lower is better.
    
    Args:
        atom_array: Biotite AtomArray containing the structure
        
    Returns:
        float: Systematic score representing the severity of steric clashes
    
    # EDUCATIONAL NOTE - Steric Repulsion and Forces
    # Atoms are not hard billiard balls, but they do have a "Van der Waals radius".
    # When two atoms get too close, their electron clouds repel each other.
    #
    # The standard "Lennard-Jones Potential" models this energy as:
    # E = 4ε [ (σ/r)^12 - (σ/r)^6 ]
    #
    # - (σ/r)^12 term: Repulsion. Rises steeply as atoms overlap (Pauli exclusion).
    # - (σ/r)^6 term: Attraction. Weak "dispersion" forces that hold non-bonded atoms together.
    #
    # In this function, we use a simplified "Soft Sphere" or "Cubic Penalty" approach.
    # Instead of shooting to infinity (which breaks optimization math), we assume
    # the energy rises polynomially once distance < sum_of_radii.
    """
    if atom_array.array_length() < 2:
        return 0.0

    # Get coordinates and elements
    coords = atom_array.coord
    elements = atom_array.element
    
    # Get Cell List for efficient neighbor search
    # Cell size should be max VdW radius * 2 roughly
    cell_list = struc.CellList(atom_array, cell_size=5.0)
    
    clash_score = 0.0
    
    # Iterate through all atoms using the cell list to find neighbors
    # This is much faster than N^2 for large structures
    for i in range(len(atom_array)):
        # Get potential neighbors within 5A
        # Note: get_atoms returns indices
        indices = cell_list.get_atoms(coords[i], radius=5.0)
        
        # Filter indices to only look at j > i to avoid double counting and self-interaction
        indices = indices[indices > i]
        
        if len(indices) == 0:
            continue
            
        atom1 = atom_array[i]
        
        # Get VdW radius for atom 1
        r1 = VAN_DER_WAALS_RADII.get(atom1.element, 1.5)
        
        for j in indices:
            atom2 = atom_array[j]
            
            # Skip atoms in same residue (simplified exclusion)
            # A full forcefield excludes 1-2, 1-3, and scaled 1-4 interactions
            # Here we just blindly skip intra-residue to avoid self-clashes from bond geometry
            if atom1.res_id == atom2.res_id:
                continue
                
            # Skip peptide bond connections (adjacent residues)
            # This is a heuristic: adjacent residues have bonded atoms that are close
            if abs(atom1.res_id - atom2.res_id) == 1:
                # Still check for severe clashes, but ignore backbone-backbone closeness
                if atom1.atom_name in ['C', 'O', 'N', 'CA'] and atom2.atom_name in ['C', 'O', 'N', 'CA']:
                    continue
            
            r2 = VAN_DER_WAALS_RADII.get(atom2.element, 1.5)
            
            # Simple Lennard-Jones-like repulsion term
            # Energy ~ (Rmin / r)^12 
            # We want a soft-ish repulsion to guide optimization
            
            dist_sq = np.sum((coords[i] - coords[j])**2)
            dist = np.sqrt(dist_sq)
            
            optimal_dist = r1 + r2
            
            if dist < optimal_dist * 0.8: # Overlap threshold
                # Severe clash
                overlap = (optimal_dist * 0.8) - dist
                # Cubic penalty for smoothness
                clash_score += (overlap * 10) ** 2
                
    return clash_score

def calculate_energy_score(atom_array: struc.AtomArray) -> float:
    """
    Calculate a total pseudo-energy score.
    Currently just wraps clash_score, but can be extended with electrostatics/solvation.
    """
    return calculate_clash_score(atom_array)
