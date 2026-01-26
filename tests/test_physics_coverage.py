
import pytest
from unittest.mock import patch, MagicMock
import sys
import synth_pdb.physics as physics

class TestPhysicsCoverage:

    def test_missing_openmm_init(self):
        """Test graceful initialization when OpenMM is not available."""
        with patch.object(physics, 'HAS_OPENMM', False):
            minimizer = physics.EnergyMinimizer()
            # Should initialize but likely set attributes to None or generic
            assert not hasattr(minimizer, 'forcefield') or minimizer.forcefield is None

    def test_missing_openmm_methods(self):
        """Test methods return False/fail gracefully without OpenMM."""
        with patch.object(physics, 'HAS_OPENMM', False):
            minimizer = physics.EnergyMinimizer()
            
            assert minimizer.minimize("dummy", "dummy") is False
            assert minimizer.add_hydrogens_and_minimize("dummy", "dummy") is False
            assert minimizer.equilibrate("dummy", "dummy") is False

    def test_forcefield_loading_error(self):
        """Test exception when forcefield fails to load."""
        with patch.object(physics, 'HAS_OPENMM', True):
            # We mock app to raise exception
            mock_app = MagicMock()
            mock_app.ForceField.side_effect = Exception("Forcefield missing")
            mock_app.OBC2 = "OBC2" 
            
            with patch.object(physics, 'app', mock_app):
                with pytest.raises(Exception, match="Forcefield missing"):
                    physics.EnergyMinimizer(forcefield_name='bad_ff')

    def test_raw_minimize_execution(self):
        """Test the _run_simulation internal worker direct calls via minimize()."""
        # This mocks the internal run to verify minimize() calls it with add_hydrogens=False
        with patch.object(physics, 'HAS_OPENMM', True):
            minimizer = physics.EnergyMinimizer()
            with patch.object(minimizer, '_run_simulation', return_value=True) as mock_run:
                minimizer.minimize("utils.pdb", "out.pdb")
                mock_run.assert_called_with("utils.pdb", "out.pdb", 0, 10.0, add_hydrogens=False)

    def test_simulation_failure_handling(self):
        """Test that _run_simulation handles exceptions gracefully."""
        with patch.object(physics, 'HAS_OPENMM', True):
            minimizer = physics.EnergyMinimizer()
            
            # Mock app.PDBFile to raise exception
            mock_app = MagicMock()
            mock_app.PDBFile.side_effect = Exception("Corrupt PDB")
            
            with patch.object(physics, 'app', mock_app):
                result = minimizer._run_simulation("bad.pdb", "out.pdb")
                assert result is False
