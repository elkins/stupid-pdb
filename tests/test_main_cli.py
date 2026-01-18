import pytest
import sys
import io
import logging # Import logging to set level for caplog
from unittest import mock
import numpy as np # Import numpy for array creation

from stupid_pdb import main
from stupid_pdb.validator import PDBValidator
from stupid_pdb.generator import create_atom_line # Import the helper function

# --- Tests for --guarantee-valid ---
def test_guarantee_valid_success(mocker, caplog):
    # Set caplog level to DEBUG to capture all relevant messages, especially the one about current_violations
    caplog.set_level(logging.DEBUG)

    # PDB content that causes steric clashes (2 violations: min_distance and VdW overlap)
    # Manually crafted to ensure it's a raw string
    clashing_pdb_content = (
        "HEADER    clashing_peptide\n" +
        create_atom_line(1, "CA", "ALA", "A", 1, 0.0, 0.0, 0.0, "C", alt_loc="", insertion_code="") + "\n" +
        create_atom_line(2, "CA", "ALA", "A", 2, 0.500, 0.0, 0.0, "C", alt_loc="", insertion_code="")
    )
        
        # Valid PDB content
        # Manually crafted to ensure it's a raw string
    valid_pdb_content = (
        "HEADER    valid_peptide\n" +
        create_atom_line(1, "CA", "GLY", "A", 1, 0.0, 0.0, 0.0, "C", alt_loc="", insertion_code="") + "\n" +
        create_atom_line(2, "CA", "GLY", "A", 2, 3.0, 0.0, 0.0, "C", alt_loc="", insertion_code="") + "\n" +
        create_atom_line(3, "CA", "GLY", "A", 3, 6.0, 0.0, 0.0, "C", alt_loc="", insertion_code="")
    )    
    mocker.patch("stupid_pdb.main.generate_pdb_content", side_effect=[
        clashing_pdb_content,
        clashing_pdb_content,
        valid_pdb_content
    ])
    
    # No need to mock PDBValidator or its methods, let the real one run.

    # Mock sys.argv to simulate CLI arguments, including --log-level DEBUG
    test_args = ["stupid_pdb", "--length", "1", "--guarantee-valid", "--max-attempts", "3", "--output", "test_gv_success.pdb", "--log-level", "DEBUG"]
    mocker.patch("sys.argv", test_args)

    # Mock sys.exit to prevent actual exit
    mocker.patch("sys.exit")

    main.main()

    # Assert that expected log messages are present
    assert "PDB generated in attempt 1 has 2 violations. Retrying..." in caplog.text
    assert "PDB generated in attempt 2 has 2 violations. Retrying..." in caplog.text
    assert "Successfully generated a valid PDB file after 3 attempts." in caplog.text
    assert "test_gv_success.pdb" in caplog.text
    sys.exit.assert_not_called() # Should not exit with error

def test_guarantee_valid_failure(mocker, caplog):
    caplog.set_level(logging.INFO) # Set to INFO to capture relevant messages
    
    # PDB content that causes steric clashes (2 violations: min_distance and VdW overlap)
    clashing_pdb_content = (
        "HEADER    clashing_peptide\n" +
        create_atom_line(1, "CA", "ALA", "A", 1, 0.0, 0.0, 0.0, "C", alt_loc="", insertion_code="") + "\n" +
        create_atom_line(2, "CA", "ALA", "A", 2, 0.500, 0.0, 0.0, "C", alt_loc="", insertion_code="")
    )

    mocker.patch("stupid_pdb.main.generate_pdb_content", return_value=clashing_pdb_content) # Always return clashing PDB

    # No need to mock PDBValidator or its methods, let the real one run.

    # Mock sys.argv
    test_args = ["stupid_pdb", "--length", "1", "--guarantee-valid", "--max-attempts", "2"]
    mocker.patch("sys.argv", test_args)

    # Mock sys.exit to check for error exit
    mock_sys_exit = mocker.patch("sys.exit")

    main.main()

    assert "PDB generated in attempt 1 has 2 violations. Retrying..." in caplog.text
    assert "PDB generated in attempt 2 has 2 violations. Retrying..." in caplog.text
    assert "Failed to generate a suitable PDB file after 2 attempts." in caplog.text
    mock_sys_exit.assert_called_once_with(1)

# --- Tests for --best-of-N ---
def test_best_of_N_selection(mocker, caplog):
    caplog.set_level(logging.INFO) # Set to INFO to capture relevant messages

    # PDB content with 2 violations (steric clash, VdW overlap)
    pdb_content_2_violations = (
        "HEADER    two_violations\n" +
        create_atom_line(1, "CA", "ALA", "A", 1, 0.0, 0.0, 0.0, "C", alt_loc="", insertion_code="") + "\n" +
        create_atom_line(2, "CA", "ALA", "A", 2, 0.500, 0.0, 0.0, "C", alt_loc="", insertion_code="")
    )
    # PDB content with 1 violation (e.g., a less severe steric clash)
    pdb_content_1_violation = (
        "HEADER    one_violation\n" +
        create_atom_line(1, "CA", "ALA", "A", 1, 0.0, 0.0, 0.0, "C", alt_loc="", insertion_code="") + "\n" +
        create_atom_line(2, "CA", "ALA", "A", 2, 1.0, 0.0, 0.0, "C", alt_loc="", insertion_code="")
    )
    # PDB content with 0 violations
    pdb_content_0_violations = (
        "HEADER    no_violations\n" +
        create_atom_line(1, "CA", "GLY", "A", 1, 0.0, 0.0, 0.0, "C", alt_loc="", insertion_code="") + "\n" +
        create_atom_line(2, "CA", "GLY", "A", 2, 3.0, 0.0, 0.0, "C", alt_loc="", insertion_code="")
    )

    mocker.patch("stupid_pdb.main.generate_pdb_content", side_effect=[
        pdb_content_2_violations, # First generated PDB will have 2 violations
        pdb_content_1_violation,  # Second generated PDB will have 1 violation
        pdb_content_0_violations, # Third generated PDB will have 0 violations
    ])

    # No need to mock PDBValidator or its methods, let the real one run.

    # Mock sys.argv
    test_args = ["stupid_pdb", "--length", "1", "--best-of-N", "3", "--output", "test_best_of_N.pdb"]
    mocker.patch("sys.argv", test_args)
    mocker.patch("sys.exit") # Should not exit with error

    main.main()

    # The actual violations found by the real PDBValidator (based on content above)
    assert "Attempt 1 yielded 2 violations" in caplog.text
    assert "Attempt 2 yielded 2 violations. Current minimum is 2." in caplog.text
    assert "Attempt 3 yielded 0 violations (new minimum)." in caplog.text
    assert "No violations found in the final PDB for" in caplog.text # Because the 0-violation PDB was chosen
    sys.exit.assert_not_called()
# --- Tests for --refine-clashes ---
def test_refine_clashes_reduces_violations(mocker, caplog):
    caplog.set_level(logging.INFO) # Set to INFO to capture relevant messages

    # PDB content that causes steric clashes (2 violations: min_distance and VdW overlap)
    initial_clashing_pdb_content = (
        "HEADER    clashing_peptide\n" +
        create_atom_line(1, "CA", "ALA", "A", 1, 0.0, 0.0, 0.0, "C", alt_loc="", insertion_code="") + "\n" +
        create_atom_line(2, "CA", "ALA", "A", 2, 0.500, 0.0, 0.0, "C", alt_loc="", insertion_code="")
    )

    # PDB content that has no steric clashes (parsed atoms for mocking tweak result)
    # Ensure coords are numpy arrays as the validator expects them for calculation
    non_clashing_parsed_atoms = [
        {"atom_number": 1, "atom_name": "CA", "alt_loc": "", "residue_name": "ALA", "chain_id": "A", "residue_number": 1, "insertion_code": "", "coords": np.array([0.0, 0.0, 0.0]), "occupancy": 1.0, "temp_factor": 0.0, "element": "C", "charge": ""},
        {"atom_number": 2, "atom_name": "CA", "alt_loc": "", "residue_name": "ALA", "chain_id": "A", "residue_number": 2, "insertion_code": "", "coords": np.array([3.0, 0.0, 0.0]), "occupancy": 1.0, "temp_factor": 0.0, "element": "C", "charge": ""},
    ]

    mocker.patch("stupid_pdb.main.generate_pdb_content", return_value=initial_clashing_pdb_content)

    # Mock _apply_steric_clash_tweak to simulate it working
    mocker.patch.object(PDBValidator, "_apply_steric_clash_tweak", return_value=non_clashing_parsed_atoms)
    
    # Mock sys.argv
    test_args = ["stupid_pdb", "--length", "1", "--output", "test_refine.pdb", "--refine-clashes", "1"]
    mocker.patch("sys.argv", test_args)
    mocker.patch("sys.exit")

    main.main()

    # The initial clashing_pdb_content will result in 2 steric clash violations.
    # After the tweak, it should be 0.
    assert "Refinement iteration 1/1. Violations: 2" in caplog.text
    assert "Refinement iteration 1: Reduced violations from 2 to 0." in caplog.text
    assert "Refinement process completed. Reduced total violations from 2 to 0." in caplog.text
    assert "No violations found in the final PDB for" in caplog.text
    sys.exit.assert_not_called()

def test_refine_clashes_no_improvement(mocker, caplog):
    caplog.set_level(logging.INFO) # Set to INFO to capture relevant messages

    # PDB content that causes steric clashes (2 violations: min_distance and VdW overlap)
    initial_clashing_pdb_content = (
        "HEADER    clashing_peptide\n" +
        create_atom_line(1, "CA", "ALA", "A", 1, 0.0, 0.0, 0.0, "C", alt_loc="", insertion_code="") + "\n" +
        create_atom_line(2, "CA", "ALA", "A", 2, 0.500, 0.0, 0.0, "C", alt_loc="", insertion_code="")
    )
    mocker.patch("stupid_pdb.main.generate_pdb_content", return_value=initial_clashing_pdb_content)

    # Mock _apply_steric_clash_tweak to return modified atoms (but still clashing)
    # Ensure coords are numpy arrays as the validator expects them for calculation
    still_clashing_parsed_atoms = [
        {"atom_number": 1, "atom_name": "CA", "alt_loc": "", "residue_name": "ALA", "chain_id": "A", "residue_number": 1, "insertion_code": "", "coords": np.array([0.0, 0.0, 0.0]), "occupancy": 1.0, "temp_factor": 0.0, "element": "C", "charge": ""},
        {"atom_number": 2, "atom_name": "CA", "alt_loc": "", "residue_name": "ALA", "chain_id": "A", "residue_number": 2, "insertion_code": "", "coords": np.array([0.6, 0.0, 0.0]), "occupancy": 1.0, "temp_factor": 0.0, "element": "C", "charge": ""}, # Still clashing
    ]
    mocker.patch.object(PDBValidator, "_apply_steric_clash_tweak", return_value=still_clashing_parsed_atoms)
    
    # Mock sys.argv
    test_args = ["stupid_pdb", "--length", "1", "--output", "test_refine_no_improve.pdb", "--refine-clashes", "2"]
    mocker.patch("sys.argv", test_args)
    mocker.patch("sys.exit")

    main.main()

    # The initial clashing_pdb_content will result in 2 steric clash violations.
    # After the tweak, it should still be 2.
    assert "Refinement iteration 1/2. Violations: 2" in caplog.text
    assert "Refinement iteration 1: No further reduction in violations (2). Stopping refinement." in caplog.text
    assert "Refinement process completed. No change in total violations (2)." in caplog.text
    assert "Final PDB has 2 violations." in caplog.text
    sys.exit.assert_not_called()
