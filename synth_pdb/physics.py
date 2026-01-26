import logging
try:
    import openmm.app as app
    import openmm as mm
    from openmm import unit
    HAS_OPENMM = True
except ImportError:
    HAS_OPENMM = False
    app = None
    mm = None
    unit = None
import sys
import os

logger = logging.getLogger(__name__)

class EnergyMinimizer:
    """
    Performs energy minimization on molecular structures using OpenMM.
    
    ### Educational Note: What is Energy Minimization?
    --------------------------------------------------
    Proteins fold into specific 3D shapes to minimize their "Gibbs Free Energy".
    A generated structure (like one built from simple geometry) often has "clashes"
    where atoms are too close (high Van der Waals repulsion) or bond angles are strained.
    
    Energy Minimization is like rolling a ball down a hill. The "Energy Landscape"
    represents the potential energy of the protein as a function of all its atom coordinates.
    The algorithm moves atoms slightly to reduce this energy, finding a local minimum
    where the structure is physically relaxed.

    ### NMR Perspective:
    In NMR structure calculation (e.g., CYANA, XPLOR-NIH), minimization is often part of
    "Simulated Annealing". Structures are calculated to satisfy experimental restraints
    (NOEs, J-couplings) and then energy-minimized to ensure good geometry.
    This module performs that final "geometry regularization" step.
    """
    
    def __init__(self, forcefield_name='amber14-all.xml', solvent_model=None):
        """
        Initialize the Minimizer with a Forcefield and Solvent Model.
        
        Args:
            forcefield_name: The "rulebook" for how atoms interact.
                             'amber14-all.xml' describes protein atoms (parameters for bond lengths,
                             angles, charges, and VdW radii).
            solvent_model:   How water is simulated. 
                             'app.OBC2' is an "Implicit Solvent" model. Instead of simulating
                             thousands of individual water molecules (Explicit Solvent),
                             it uses a mathematical continuum to approximate water's dielectric 
                             shielding and hydrophobic effects. This is much faster.
                             
                             ### NMR Note:
                             Since NMR is performed in solution (not crystals), implicit solvent 
                             aims to approximate that solution environment, distinct from the
                             vacuum or crystal lattice assumptions of other methods.
        """
        if not HAS_OPENMM:
             # We allow initialization to pass so the module can be imported
             # but methods will check HAS_OPENMM
             pass

        # Set default solvent model
        if solvent_model is None and HAS_OPENMM:
            solvent_model = app.OBC2

        self.forcefield_name = forcefield_name
        self.water_model = 'amber14/tip3pfb.xml' 
        self.solvent_model = solvent_model
        
        if not HAS_OPENMM:
            return

        # Map solvent models to their parameter files in OpenMM
        # These need to be loaded alongside the main forcefield
        # Standard OpenMM paths often have 'implicit/' at the root
        solvent_xml_map = {
            app.OBC2: 'implicit/obc2.xml',
            app.OBC1: 'implicit/obc1.xml',
            app.GBn:  'implicit/gbn.xml',
            app.GBn2: 'implicit/gbn2.xml',
            app.HCT:  'implicit/hct.xml',
        }

        # Build list of XML files to load
        ff_files = [self.forcefield_name, self.water_model]
        
        if self.solvent_model in solvent_xml_map:
            ff_files.append(solvent_xml_map[self.solvent_model])
        
        try:
            # The ForceField object loads the definitions of atom types and parameters.
            self.forcefield = app.ForceField(*ff_files)
        except Exception as e:
            logger.error(f"Failed to load forcefield {ff_files}: {e}")
            raise

    def minimize(self, pdb_file_path: str, output_path: str, max_iterations: int = 0, tolerance: float = 10.0) -> bool:
        """
        Minimizes the energy of a structure already containing correct atoms (including Hydrogens).
        
        Args:
            pdb_file_path: Input PDB path.
            output_path: Output PDB path.
            max_iterations: Limit steps (0 = until convergence).
            tolerance: Target energy convergence threshold (kJ/mol).
        """
        if not HAS_OPENMM:
             logger.error("Cannot minimize: OpenMM not found.")
             return False

        # This method assumes the input PDB is perfect (has Hydrogens, correct names).
        # See 'add_hydrogens_and_minimize' for the robust version used by synth-pdb.
        return self._run_simulation(pdb_file_path, output_path, max_iterations, tolerance, add_hydrogens=False)

    def add_hydrogens_and_minimize(self, pdb_file_path: str, output_path: str) -> bool:
        """
        Robust minimization pipeline: Adds Hydrogens -> Creates/Minimizes System -> Saves Result.
        
        ### Why Add Hydrogens?
        X-ray crystallography often doesn't resolve hydrogen atoms because they have very few electrons.
        However, Molecular Dynamics forcefields (like Amber) are explicitly "All-Atom". They REQUIRE
        hydrogens to calculate bond angles and electrostatics (h-bonds) correctly.
        
        ### NMR Perspective:
        Unlike X-ray, NMR relies entirely on the magnetic spin of protons (H1). Hydrogens are
        the "eyes" of NMR. Correctly placing them is critical not just for physics but for
        predicting NOEs (Nuclear Overhauser Effects) which depend on H-H distances.
        We use `app.Modeller` to "guess" the standard positions of hydrogens at specific pH (7.0).
        """
        if not HAS_OPENMM:
             logger.error("Cannot add hydrogens: OpenMM not found.")
             return False
             
        return self._run_simulation(pdb_file_path, output_path, add_hydrogens=True)

    def equilibrate(self, pdb_file_path: str, output_path: str, steps: int = 1000) -> bool:
        """
        Run Thermal Equilibration (MD) at 300K.
        
        Args:
            pdb_file_path: Input PDB/File path.
            output_path: Output PDB path.
            steps: Number of MD steps (2 fs per step). 1000 steps = 2 ps.
            
        Returns:
            True if successful.
        """
        if not HAS_OPENMM:
             logger.error("Cannot equilibrate: OpenMM not found.")
             return False
        
        # We always add hydrogens for MD
        return self._run_simulation(pdb_file_path, output_path, add_hydrogens=True, equilibration_steps=steps)

    def _run_simulation(self, input_path, output_path, max_iterations=0, tolerance=10.0, add_hydrogens=True, equilibration_steps=0):
        """Internal worker for setting up and running the OpenMM Simulation."""
        logger.info(f"Processing physics for {input_path}...")
        
        try:
            # 1. Load the PDB file structure (topology and positions)
            pdb = app.PDBFile(input_path)
            
            # 2. Prepare the Topology (Add Hydrogens if requested)
            if add_hydrogens:
                logger.info("Adding missing hydrogens (protonation state @ pH 7.0)...")
                modeller = app.Modeller(pdb.topology, pdb.positions)
                modeller.addHydrogens(self.forcefield, pH=7.0)
                topology = modeller.topology
                positions = modeller.positions
            else:
                topology = pdb.topology
                positions = pdb.positions

            # 3. Create the 'System'
            # The System object connects the Topology (atoms/bonds) to the Forcefield (physics rules).
            # It calculates all forces acting on every atom.
            logger.debug("Creating OpenMM System (applying forcefield parameters)...")
            try:
                system = self.forcefield.createSystem(
                    topology,
                    nonbondedMethod=app.NoCutoff, # No cutoff for vacuum/implicit (calculates ALL pairwise interactions)
                    constraints=app.HBonds,       # Constrain Hydrogen bond lengths (allows larger time steps in MD)
                    implicitSolvent=self.solvent_model # Continuum water model
                )
            except ValueError as ve:
                # Fallback logic for when implicit solvent fails (common with some forcefield combos)
                if "implicitSolvent" in str(ve):
                    logger.warning(f"Implicit Solvent parameters not found for this forcefield configuration. Using Vacuum electrostatics instead (standard fallback).")
                    system = self.forcefield.createSystem(
                        topology,
                        nonbondedMethod=app.NoCutoff,
                        constraints=app.HBonds
                        # implicitSolvent defaults to None (Vacuum)
                    )
                else:
                    raise ve
            
            # 4. Create the Integrator
            # An Integrator is the math engine that moves atoms based on forces (F=ma).
            # 'LangevinIntegrator' simulates a heat bath (friction + random collisions) to maintain temperature.
            #
            # Educational Note:
            # For pure energy minimization (finding the nearest valley), we technically don't need a
            # thermostat like Langevin because we aren't simulating time-resolved motion yet.
            # However, OpenMM requires an Integrator to define the Simulation context.
            # If we were to run "simulation.step(1000)", this integrator would generate
            # realistic Brownian motion, simulating thermal fluctuations.
            integrator = mm.LangevinIntegrator(
                300 * unit.kelvin,       # Temperature (Room temp)
                1.0 / unit.picosecond,   # Friction coefficient
                2.0 * unit.femtoseconds  # Time step
            )
            
            # 5. Create the Simulation context
            simulation = app.Simulation(topology, system, integrator)
            simulation.context.setPositions(positions)
            
            # Report Energy BEFORE Minimization
            state_initial = simulation.context.getState(getEnergy=True)
            e_init = state_initial.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
            try:
                logger.info(f"Initial Potential Energy: {e_init:.2f} kJ/mol")
            except Exception:
                # Value might be non-numeric (e.g. during testing/mocking) or NaN
                logger.info(f"Initial Potential Energy: {e_init} kJ/mol")
            
            try:
                if e_init > 1e6:
                    logger.info("  -> High initial energy detected due to steric clashes. Minimization will resolve this.")
            except Exception:
                pass
            
            # 6. Run Energy Minimization (Gradient Descent)
            logger.info("Minimizing energy...")
            # Tolerance is Force threshold (kJ/mol/nm)
            simulation.minimizeEnergy(maxIterations=max_iterations, tolerance=tolerance*unit.kilojoule/(unit.mole*unit.nanometer))
            
            # 7. Run MD Equilibration (Optional)
            if equilibration_steps > 0:
                logger.info(f"Running Thermal Equilibration ({equilibration_steps} steps at 300K)...")
                # EDUCATIONAL NOTE - Molecular Dynamics (MD):
                # We now "unfreeze" the protein. The LangevinIntegrator simulates 
                # collisions with a solvent bath at 300 Kelvin.
                # This shakes the atoms out of shallow local minima and enables
                # realistic breathing motions.
                simulation.step(equilibration_steps)
            
            # 8. Report Final Energy
            state_final = simulation.context.getState(getEnergy=True, getPositions=True)
            e_final = state_final.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
            try:
                logger.info(f"Final Potential Energy:   {e_final:.2f} kJ/mol")
            except Exception:
                logger.info(f"Final Potential Energy:   {e_final} kJ/mol")
            
            # 7. Save Result
            with open(output_path, 'w') as f:
                app.PDBFile.writeFile(simulation.topology, state_final.getPositions(), f)
                
            return True

        except Exception as e:
            logger.error(f"Physics simulation failed: {e}")
            if "template" in str(e).lower():
                logger.error("Error Hint: OpenMM couldn't match residues to the forcefield. This usually means atoms are missing or named incorrectly.")
            return False
