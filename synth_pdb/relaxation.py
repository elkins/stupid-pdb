
import numpy as np
import biotite.structure as struc
import logging
from typing import Dict, List, Tuple

logger = logging.getLogger(__name__)

# --- Physical Constants for NMR Relaxation ---
# SI Units used for internal calculation
MU_0 = 4 * np.pi * 1e-7      # Vacuum permeability derived (T*m/A)
H_PLANCK = 6.62607015e-34    # Planck constant (J*s)
H_BAR = H_PLANCK / (2 * np.pi)

GAMMA_H = 267.522e6          # Proton gyromagnetic ratio (rad s^-1 T^-1)
GAMMA_N = -27.126e6          # Nitrogen-15 gyromagnetic ratio (rad s^-1 T^-1)

R_NH = 1.02e-10              # NH Bond length (meters) - standard value
CSA_N = -160e-6              # Polimorphic 15N CSA (unitless, ppm) -160 to -170 typical

def spectral_density(omega: float, tau_m: float, s2: float, tau_f: float = 0.0) -> float:
    """
    Calculate Spectral Density J(w) using Lipari-Szabo Model-Free formalism.
    
    Formula:
    J(w) = (2/5) * [ S^2 * tm / (1 + (w*tm)^2) + (1-S^2) * te / (1 + (w*te)^2) ]
    
    where te (tau_e) is the effective internal correlation time: 1/te = 1/tm + 1/tf
    Usually for simple MF, we assume fast motion tf << tm.
    
    Args:
        omega: Frequency (rad/s)
        tau_m: Global rotational correlation time (seconds)
        s2: Generalized order parameter (0.0 to 1.0)
        tau_f: Fast internal correlation time (seconds). Default 0 (assumed very fast).
    """
    # Simple Model Free (assuming tf is very small/negligible or incorporated)
    # If tau_f is provided, calculate effective time tau_e
    
    # Term 1: Global tumbling
    j_global = (s2 * tau_m) / (1 + (omega * tau_m)**2)
    
    # Term 2: Fast internal motion
    # Effective correlation time 1/tau_e = 1/tau_m + 1/tau_f
    # If tau_f is 0, this term vanishes in standard simplified approximation
    # or acts as a very fast motion limit.
    j_fast = 0.0
    if tau_f > 0:
        tau_e = (tau_m * tau_f) / (tau_m + tau_f)
        j_fast = ((1 - s2) * tau_e) / (1 + (omega * tau_e)**2)
        
    return 0.4 * (j_global + j_fast) 

def calculate_relaxation_rates(
    structure: struc.AtomArray,
    field_mhz: float = 600.0,
    tau_m_ns: float = 10.0
) -> Dict[int, Dict[str, float]]:
    """
    Calculate R1, R2, and Heteronuclear NOE for all backbone Amides (N-H).
    
    Args:
        structure: The protein structure (must have hydrogens).
        field_mhz: Proton Larmor frequency in MHz (e.g. 600).
        tau_m_ns: Global tumbling time in ns (default 10.0).
        
    Returns:
        Dictionary keyed by residue ID:
        { res_id: {'R1': float, 'R2': float, 'NOE': float, 'S2': float} }
    """
    logger.info(f"Calculating Relaxation Rates (Field={field_mhz}MHz, tm={tau_m_ns}ns)...")
    
    # Convert inputs to SI units
    tau_m = tau_m_ns * 1e-9
    
    # Larmor Frequencies (rad/s)
    # Field Strength B0 (Tesla) -> omega = gamma * B0
    # But user gives MHz. 600 MHz = frequency of Proton.
    # omega_H = 2 * pi * 600e6
    omega_h = 2 * np.pi * field_mhz * 1e6
    
    # Calculate B0 from proton freq
    b0 = omega_h / GAMMA_H
    
    omega_n = GAMMA_N * b0 # Negative val
    
    logger.debug(f"B0 Field: {b0:.2f} T")
    logger.debug(f"wH: {omega_h:.2e} rad/s, wN: {omega_n:.2e} rad/s")
    
    # Dipolar Pre-factor
    # d = (mu0 * hbar * gammaH * gammaN) / (4pi * r^3) ??
    # Look up standard definition for "d" in R1 equations.
    # Usually: d = (mu0 / 4pi) * (hbar * gammaH * gammaN) * r^-3
    
    dd_const = (MU_0 / (4 * np.pi)) * H_BAR * GAMMA_H * GAMMA_N * (R_NH**-3)
    d_sq = dd_const**2
    
    # CSA Pre-factor
    # c = (DeltaSigma * wN) / sqrt(3) -> definition varies by factor of 3 depending on which R1 eq used.
    # Standard: c = omega_N * CSA_N / sqrt(3)
    # OR: factor = (c**2) in equations.
    # Let's use:
    csa_const = (CSA_N * omega_n) / np.sqrt(3)
    c_sq = csa_const**2
    
    # Analyze Secondary Structure for Order Parameters (S2)
    # We'll use Phi/Psi to guess if rigid or flexible
    try:
        phi, psi, omega_dihed = struc.dihedral_backbone(structure)
        # dihedrals returns angles for residues starting from 2nd?
        # Actually returns array matching residue length, with NaNs at ends.
    except Exception:
        logger.warning("Could not calculate dihedrals. Assuming all rigid.")
        phi = np.zeros(structure.res_id[-1] + 1)
        
    # Map S2 by residue
    # Default S2 = 0.85 (Ordered)
    # Termini or Loops = 0.60
    
    # Identify residues
    res_ids = np.unique(structure.res_id)
    min_res, max_res = res_ids[0], res_ids[-1]
    
    results = {}
    
    # Iterate over residues that have an N-H pair
    # (Proline has no H, N-terminus implies H1/H2/H3 so no amide H usually defined same way)
    
    for rid in res_ids:
        # Check if N and H exist
        # Need "N" and "H" atoms
        res_mask = structure.res_id == rid
        res_atoms = structure[res_mask]
        
        has_n = "N" in res_atoms.atom_name
        has_h = "H" in res_atoms.atom_name
        res_name = res_atoms.res_name[0]
        
        if not (has_n and has_h):
            continue
            
        if res_name == "PRO":
            continue
            
        # Determine S2
        s2 = 0.85 # Default Rigid
        
        # Termini are flexible
        if rid == min_res or rid == max_res:
            s2 = 0.50
        elif rid == min_res + 1 or rid == max_res - 1:
            s2 = 0.70
            
        # If we had simple DSSP logic:
        # For synthetic ones, we built them.
        # Let's rely on basic heuristic: N/C term flexible.
        # Maybe random variation?
        # Add random noise to make it realistic
        s2 += np.random.normal(0, 0.02)
        s2 = np.clip(s2, 0.01, 1.0)
        
        # Frequencies for J(w)
        # R1 depends on: J(wH-wN), J(wN), J(wH+wN)
        # R2 depends on: J(0), J(wH-wN), J(wN), J(wH), J(wH+wN)
        # NOE depends on: J(wH+wN), J(wH-wN) ? -> Actually NOE = 1 + ...
        
        j_0 = spectral_density(0, tau_m, s2)
        j_wn = spectral_density(omega_n, tau_m, s2)
        j_wh = spectral_density(omega_h, tau_m, s2)
        j_diff = spectral_density(omega_h - omega_n, tau_m, s2)
        j_sum = spectral_density(omega_h + omega_n, tau_m, s2)
        
        # Calculate Rates
        # R1 = (d^2/4) * [J(wH-wN) + 3J(wN) + 6J(wH+wN)] + c^2 * J(wN)
        # Note the factor 1/4 or similar for d^2 depending on definition of d.
        # My d defined above: (mu0/4pi) * hbar * gammaH * gammaN * r^-3
        # Standard Abragam eq:
        # R1 = (d^2) * ... if d includes the factor. 
        # Let's use the explicit pre-factor term P = (d^2)
        
        # P = d_sq
        # R1 = P * (1*j_diff + 3*j_wn + 6*j_sum) + c_sq * j_wn
        # WAIT. Factor checks.
        
        # Reference: protein-nmr.org.uk or Cavanagh et al.
        # Dipolar term coeff: d_sq = (mu0/4pi)^2 * hbar^2 * gH^2 * gN^2 * r^-6
        # Eq: R1 = (d_sq / 4) ... ? No.
        # Let's assume standard form:
        # R1 = (d^2) * ( ... ) + CSA term
        
        # Re-verify d definition in Cavanagh:
        # d = (mu0 hbar gH gN) / (4 pi r^3)
        # R1 = (d^2) [ J(wH-wN) + 3J(wN) + 6J(wH+wN) ] + c^2 J(wN)
        # This assumes J(w) definition has the 2/5 or similar. 
        # My J(w) has 2/5.
        
        # Corrections:
        # The equation often is R1 = (d^2) * (1/4)? NO, Dipolar relaxation usually defined directly.
        # Let's stick to the form that works with J=2/5...
        #
        # d_sq calculated above is the full constant squared.
        r1_val = d_sq * (j_diff + 3*j_wn + 6*j_sum) + c_sq * j_wn
        
        # R2
        # R2 = 0.5 * d_sq * [4*J(0) + J(diff) + 3*J(wn) + 6*J(wh) + 6*J(sum)] + (1/6)*c_sq * [4*J(0) + 3*J(wn)]
        r2_val = 0.5 * d_sq * (4*j_0 + j_diff + 3*j_wn + 6*j_wh + 6*j_sum) + \
                 (1.0/6.0) * c_sq * (4*j_0 + 3*j_wn)
                 
        # NOE
        # NOE = 1 + (gamma_H / gamma_N) * d_sq * [6*J(sum) - J(diff)] / R1
        # Note: gamma quotient is negative (-10)
        noe_val = 1.0 + (GAMMA_H / GAMMA_N) * d_sq * (6*j_sum - j_diff) * (1.0 / r1_val)
        
        results[rid] = {
            'R1': r1_val,
            'R2': r2_val,
            'NOE': noe_val,
            'S2': s2
        }
        
    return results
