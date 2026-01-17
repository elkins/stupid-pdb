import unittest
import logging
import re
from stupid_pdb.generator import _resolve_sequence, generate_pdb_content, CA_DISTANCE
from stupid_pdb.data import STANDARD_AMINO_ACIDS, AMINO_ACID_ATOMS, ONE_TO_THREE_LETTER_CODE

# Suppress logging during tests to keep output clean
logging.getLogger().setLevel(logging.CRITICAL)


class TestGenerator(unittest.TestCase):

    # --- Tests for _get_sequence ---
    def test_get_sequence_random_length(self):
        """Test if random sequence generation has the correct length."""
        for length in [1, 5, 10, 100]:
            sequence = _resolve_sequence(length=length, user_sequence_str=None)
            self.assertEqual(len(sequence), length)

    def test_get_sequence_random_empty(self):
        """Test random empty sequence request."""
        sequence = _resolve_sequence(length=0, user_sequence_str=None)
        self.assertEqual(len(sequence), 0)
        sequence = _resolve_sequence(length=-5, user_sequence_str=None)
        self.assertEqual(len(sequence), 0)

    def test_get_sequence_random_amino_acids(self):
        """Test if all elements in random sequence are valid amino acids."""
        sequence = _resolve_sequence(length=20, user_sequence_str=None)
        for aa in sequence:
            self.assertIn(aa, STANDARD_AMINO_ACIDS)

    def test_get_sequence_from_1_letter_code(self):
        """Test parsing of a valid 1-letter code sequence."""
        sequence_str = "AGV"
        expected_sequence = ["ALA", "GLY", "VAL"]
        sequence = _resolve_sequence(length=0, user_sequence_str=sequence_str) # length should be ignored
        self.assertEqual(sequence, expected_sequence)

    def test_get_sequence_from_3_letter_code(self):
        """Test parsing of a valid 3-letter code sequence."""
        sequence_str = "ALA-GLY-VAL"
        expected_sequence = ["ALA", "GLY", "VAL"]
        sequence = _resolve_sequence(length=0, user_sequence_str=sequence_str)
        self.assertEqual(sequence, expected_sequence)
    
    def test_get_sequence_from_mixed_case(self):
        """Test parsing of mixed-case sequence strings."""
        sequence_str_1 = "aGv"
        expected_sequence_1 = ["ALA", "GLY", "VAL"]
        self.assertEqual(_resolve_sequence(length=0, user_sequence_str=sequence_str_1), expected_sequence_1)

        sequence_str_2 = "Ala-GlY-vAl"
        expected_sequence_2 = ["ALA", "GLY", "VAL"]
        self.assertEqual(_resolve_sequence(length=0, user_sequence_str=sequence_str_2), expected_sequence_2)

    def test_get_sequence_invalid_1_letter_code(self):
        """Test handling of invalid 1-letter code sequence."""
        sequence_str = "AXG"
        with self.assertRaisesRegex(ValueError, "Invalid 1-letter amino acid code: X"):
            _resolve_sequence(length=0, user_sequence_str=sequence_str)

    def test_get_sequence_invalid_3_letter_code(self):
        """Test handling of invalid 3-letter code sequence."""
        sequence_str = "ALA-XYZ-VAL"
        with self.assertRaisesRegex(ValueError, "Invalid 3-letter amino acid code: XYZ"):
            _resolve_sequence(length=0, user_sequence_str=sequence_str)

    def test_get_sequence_plausible_frequencies(self):
        """
        Test if random sequence generation with plausible frequencies
        adheres to the expected distribution within a tolerance.
        """
        from stupid_pdb.data import AMINO_ACID_FREQUENCIES
        test_length = 10000
        tolerance = 0.02 # 2% deviation allowed

        sequence = _resolve_sequence(length=test_length, use_plausible_frequencies=True)
        self.assertEqual(len(sequence), test_length)

        # Calculate observed frequencies
        observed_counts = {aa: sequence.count(aa) for aa in AMINO_ACID_FREQUENCIES.keys()}
        observed_frequencies = {aa: count / test_length for aa, count in observed_counts.items()}

        # Compare observed with expected frequencies
        for aa, expected_freq in AMINO_ACID_FREQUENCIES.items():
            observed_freq = observed_frequencies.get(aa, 0.0)
            self.assertAlmostEqual(observed_freq, expected_freq, delta=tolerance,
                                   msg=f"Frequency for {aa} (Observed: {observed_freq:.4f}, Expected: {expected_freq:.4f}) out of tolerance.")

    # --- Tests for generate_pdb_content (general) ---
    def test_generate_pdb_content_empty_length(self):
        """Test PDB content generation for zero or negative length when no sequence is provided."""
        with self.assertRaisesRegex(ValueError, "Length must be a positive integer when no sequence is provided."):
            generate_pdb_content(length=0, sequence_str=None)
        with self.assertRaisesRegex(ValueError, "Length must be a positive integer when no sequence is provided."):
            generate_pdb_content(length=-5, sequence_str=None)
    
    def test_generate_pdb_content_empty_sequence_str(self):
        """Test PDB content generation with an empty sequence string."""
        content = generate_pdb_content(length=0, sequence_str="")
        self.assertEqual(content, "")


    # --- Tests for generate_pdb_content (CA only) ---
    def test_generate_pdb_content_num_lines_ca_only(self):
        """Test if the generated PDB content (CA only) has the correct number of ATOM lines."""
        for length in [1, 5, 10, 50]:
            content = generate_pdb_content(length=length, full_atom=False)
            lines = content.strip().split("\n")
            
            atom_lines = [line for line in lines if line.startswith("ATOM")]
            self.assertEqual(len(atom_lines), length) # Verify actual ATOM lines count
            
            for line in atom_lines:
                self.assertTrue(line.startswith("ATOM"))
                self.assertIn("CA", line[12:16]) # Check for CA atom name

    def test_generate_pdb_content_coordinates_ca_only(self):
        """Test if atom coordinates are correctly generated for a linear CA-only chain."""
        length = 5
        content = generate_pdb_content(length=length, full_atom=False)
        lines = content.strip().split("\n")
        
        atom_lines = [line for line in lines if line.startswith("ATOM")]
        expected_x_coords = [i * CA_DISTANCE for i in range(length)]

        for i, line in enumerate(atom_lines):
            # Regex to extract x_coord (8 characters, float)
            match = re.search(r"(\s+\d+\.\d+)\s+\d+\.\d+\s+\d+\.\d+", line)
            if match:
                # Extract the x_coord string and convert to float
                x_coord = float(match.group(1).strip())
                self.assertAlmostEqual(x_coord, expected_x_coords[i], places=3)
            else:
                self.fail(f"Could not parse coordinates from line: {line}")

    def test_generate_pdb_content_atom_residue_numbers_ca_only(self):
        """Test if atom and residue numbers are sequential for CA-only."""
        length = 3
        content = generate_pdb_content(length=length, full_atom=False)
        lines = content.strip().split("\n")
        atom_lines = [line for line in lines if line.startswith("ATOM")]

        for i, line in enumerate(atom_lines):
            atom_num_str = line[6:11].strip()  # ATOM number is chars 7-11 (0-indexed)
            res_num_str = line[22:26].strip()  # Residue number is chars 23-26

            self.assertEqual(int(atom_num_str), i + 1)
            self.assertEqual(int(res_num_str), i + 1)

    def test_generate_pdb_content_residue_names_ca_only(self):
        """Test if residue names are valid for CA-only."""
        length = 5
        content = generate_pdb_content(length=length, full_atom=False)
        lines = content.strip().split("\n")
        atom_lines = [line for line in lines if line.startswith("ATOM")]


        for line in atom_lines:
            res_name = line[17:20].strip()  # Residue name is chars 18-20
            self.assertIn(res_name, STANDARD_AMINO_ACIDS)

    # --- Tests for generate_pdb_content (full_atom=True) ---
    def test_generate_pdb_content_full_atom_more_atoms(self):
        """Test that full_atom generates more atoms than CA-only."""
        length = 1
        ca_only_content = generate_pdb_content(length=length, full_atom=False)
        full_atom_content = generate_pdb_content(length=length, full_atom=True)
        self.assertGreater(len([line for line in full_atom_content.strip().split("\n") if line.startswith("ATOM")]), 
                           len([line for line in ca_only_content.strip().split("\n") if line.startswith("ATOM")]))

    def test_generate_pdb_content_full_atom_backbone_atoms(self):
        """Test for the presence of N, C, O backbone atoms in full_atom mode."""
        length = 1
        content = generate_pdb_content(length=length, full_atom=True)
        lines = [line for line in content.strip().split("\n") if line.startswith("ATOM")]

        atom_names = {line[12:16].strip() for line in lines} # Extract atom names
        
        self.assertIn("N", atom_names)
        self.assertIn("CA", atom_names)
        self.assertIn("C", atom_names)
        self.assertIn("O", atom_names)

    def test_generate_pdb_content_full_atom_side_chain_atoms(self):
        """Test for the presence of side-chain atoms (e.g., CB for ALA) in full_atom mode."""
        length = 10
        content = generate_pdb_content(length=length, full_atom=True)
        lines = [line for line in content.strip().split("\n") if line.startswith("ATOM")]

        atom_names = {line[12:16].strip() for line in lines} # Extract all atom names
        residue_names = {line[17:20].strip() for line in lines} # Extract all residue names

        has_cb_amino_acid = False
        for res in STANDARD_AMINO_ACIDS:
            if res != "GLY" and AMINO_ACID_ATOMS.get(res):
                for atom_def in AMINO_ACID_ATOMS[res]:
                    if atom_def['name'] == 'CB':
                        if res in residue_names:
                            has_cb_amino_acid = True
                            break
            if has_cb_amino_acid:
                break
        
        if has_cb_amino_acid:
            self.assertIn("CB", atom_names, "Expected CB atom in full_atom output for residues like ALA")
        else:
            logging.warning("Test `test_generate_pdb_content_full_atom_side_chain_atoms` could not find an amino acid with a CB atom in the random sequence of length %d. Test passed conditionally.", length)
            
    # --- Tests for PDB Header, TER, END records ---
    def test_generate_pdb_content_no_unintended_blank_lines(self):
        """Test that there are no unintended blank lines in the PDB content."""
        content = generate_pdb_content(length=5)
        lines = content.split("\n")
        
        non_trailing_blank_lines_count = 0
        for i, line in enumerate(lines):
            # Only count blank lines that are not the very last element (potential trailing newline from .join)
            if not line.strip() and i < len(lines) - 1:
                non_trailing_blank_lines_count += 1
        
        # The test should FAIL if it finds any unintended blank lines.
        # We expect 0 unintended blank lines.
        self.assertEqual(non_trailing_blank_lines_count, 0, 
                         f"Found {non_trailing_blank_lines_count} unintended blank lines. Content:\n{content}")

        # Also keep the check for total non-empty lines for overall content structure validation
        non_empty_lines = [line for line in lines if line.strip()]
        expected_content_lines = 19 
        self.assertEqual(len(non_empty_lines), expected_content_lines, 
                         f"Expected {expected_content_lines} non-empty lines, but found {len(non_empty_lines)}. Content:\n{content}")

    def test_generate_pdb_content_header_present(self):
        """Test if the PDB header is present at the beginning."""
        content = generate_pdb_content(length=1)
        lines = content.split("\n")
        self.assertTrue(lines[0].startswith("HEADER"))
        self.assertTrue(lines[1].startswith("TITLE"))

    def test_generate_pdb_content_ter_present(self):
        """Test if the TER record is present and correctly formatted."""
        length = 3
        content = generate_pdb_content(length=length)
        lines = content.strip().split("\n")
        
        ter_line = [line for line in lines if line.startswith("TER")][-1]
        self.assertIsNotNone(ter_line)
        self.assertTrue(ter_line.startswith("TER"))
        self.assertIn("A", ter_line[21]) # Check chain ID

        # Extract last atom number from the preceding ATOM line
        atom_lines = [line for line in lines if line.startswith("ATOM")]
        last_atom_line = atom_lines[-1]
        expected_atom_num = int(last_atom_line[6:11].strip())

        # The TER record atom number should be one greater than the last ATOM record
        ter_atom_num = int(ter_line[6:11].strip())
        self.assertEqual(ter_atom_num, expected_atom_num + 1)
        
        # Check residue name and number of the TER record matches the last ATOM record
        expected_res_name = last_atom_line[17:20].strip()
        expected_res_num = int(last_atom_line[22:26].strip())
        self.assertEqual(ter_line[17:20].strip(), expected_res_name)
        self.assertEqual(int(ter_line[22:26].strip()), expected_res_num)


    def test_generate_pdb_content_end_present(self):
        """Test if the END record is present at the very end."""
        content = generate_pdb_content(length=1)
        lines = content.strip().split("\n")
        self.assertEqual(lines[-1], "END")

    # --- New tests for generate_pdb_content with sequence_str ---
    def test_generate_pdb_content_with_sequence_1_letter(self):
        """Test PDB content generation with a user-provided 1-letter sequence."""
        sequence_str = "AGV"
        expected_sequence = ["ALA", "GLY", "VAL"]
        content = generate_pdb_content(sequence_str=sequence_str)
        lines = content.strip().split("\n")
        atom_lines = [line for line in lines if line.startswith("ATOM")]
        
        self.assertEqual(len(atom_lines), len(expected_sequence))
        for i, line in enumerate(atom_lines):
            res_name = line[17:20].strip()
            self.assertEqual(res_name, expected_sequence[i])

    def test_generate_pdb_content_with_sequence_3_letter(self):
        """Test PDB content generation with a user-provided 3-letter sequence."""
        sequence_str = "ALA-GLY-VAL"
        expected_sequence = ["ALA", "GLY", "VAL"]
        content = generate_pdb_content(sequence_str=sequence_str)
        lines = content.strip().split("\n")
        atom_lines = [line for line in lines if line.startswith("ATOM")]
        
        self.assertEqual(len(atom_lines), len(expected_sequence))
        for i, line in enumerate(atom_lines):
            res_name = line[17:20].strip()
            self.assertEqual(res_name, expected_sequence[i])

    def test_generate_pdb_content_sequence_overrides_length(self):
        """Test that provided sequence's length overrides the 'length' parameter."""
        sequence_str = "AG" # Length 2
        length_param = 5   # Should be ignored
        content = generate_pdb_content(length=length_param, sequence_str=sequence_str)
        lines = content.strip().split("\n")
        atom_lines = [line for line in lines if line.startswith("ATOM")]
        self.assertEqual(len(atom_lines), 2) # Should be 2, not 5

    def test_generate_pdb_content_invalid_sequence_str_raises_error(self):
        """Test that invalid sequence string raises ValueError during PDB generation."""
        invalid_sequence_str = "AXG"
        with self.assertRaises(ValueError):
            generate_pdb_content(sequence_str=invalid_sequence_str)

if __name__ == '__main__':
    unittest.main()
