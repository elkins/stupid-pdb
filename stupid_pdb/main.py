"""
CLI entry point for the stupid-pdb tool.
"""

import argparse
import logging
import datetime
import os
import sys

from .generator import generate_pdb_content
from .validator import PDBValidator

# Get a logger for this module
logger = logging.getLogger(__name__)


def main():
    """
    Main function to parse arguments and generate the PDB file.
    """
    parser = argparse.ArgumentParser(
        description="Generate a PDB file with a random linear amino acid sequence."
    )
    parser.add_argument(
        "--length",
        type=int,
        default=10,
        help="Length of the amino acid sequence (number of residues).",
    )
    parser.add_argument(
        "--output",
        type=str,
        help="Optional: Output PDB filename. If not provided, a default name will be generated.",
    )
    parser.add_argument(
        "--log-level",
        type=str,
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        help="Set the logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL).",
    )
    parser.add_argument(
        "--full-atom",
        action="store_true",
        help="Generate a full atomic representation (N, CA, C, O, and side chain atoms) instead of only C-alpha atoms.",
    )
    parser.add_argument(
        "--sequence",
        type=str,
        help="Specify an amino acid sequence (e.g., 'AGV' or 'ALA-GLY-VAL'). Overrides random generation.",
    )
    parser.add_argument(
        "--plausible-frequencies",
        action="store_true",
        help="Use biologically plausible amino acid frequencies for random sequence generation (ignored if --sequence is provided).",
    )
    parser.add_argument(
        "--validate",
        action="store_true",
        help="Run validation checks (bond lengths and angles, Ramachandran) on the generated PDB.",
    )
    parser.add_argument(
        "--guarantee-valid",
        action="store_true",
        help="If set, repeatedly generate PDB files until a valid one (no violations) is produced. Implies --validate. Will stop after --max-attempts if no valid PDB is found.",
    )
    parser.add_argument(
        "--max-attempts",
        type=int,
        default=100,
        help="Maximum number of regeneration attempts when --guarantee-valid is set.",
    )
    parser.add_argument(
        "--best-of-N",
        type=int,
        default=1,
        help="Generate N PDBs, validate each, and select the one with the fewest violations. Implies --validate. Overrides --guarantee-valid.",
    )

    args = parser.parse_args()

    # Set the logging level based on user input
    log_level = getattr(logging, args.log_level.upper(), None)
    if not isinstance(log_level, int):
        raise ValueError(f"Invalid log level: {args.log_level}")
    logging.getLogger().setLevel(log_level)
    logger.debug("Logging level set to %s.", args.log_level.upper())

    logger.info("Starting PDB file generation process.")
    logger.debug(
        "Parsed arguments: length=%s, output='%s', full_atom=%s, sequence='%s', plausible_frequencies=%s, validate=%s",
        args.length,
        args.output,
        args.full_atom,
        args.sequence,
        args.plausible_frequencies,
        args.validate,
    )

    # If --best-of-N is set, it overrides --guarantee-valid and implies --validate.
    if args.best_of_N > 1:
        args.validate = True
        args.guarantee_valid = False # Disable guarantee-valid if best-of-N is used
        logger.info(f"--best-of-N is set to {args.best_of_N}. Generating multiple PDBs to find the one with fewest violations.")
    elif args.guarantee_valid: # Only apply if best-of-N is not active
        args.validate = True
        logger.info("--guarantee-valid is set. Will attempt to generate a valid PDB.")

    # Validate length only if no sequence is provided
    if args.sequence is None and (args.length is None or args.length <= 0):
        logger.error("Length must be a positive integer when no sequence is provided.")
        sys.exit(1)

    length_for_generator = args.length if args.sequence is None else None

    final_pdb_content = None
    final_violations = None
    min_violations_count = float('inf')

    generation_attempts = 1 if not args.guarantee_valid and args.best_of_N <= 1 else args.max_attempts
    if args.best_of_N > 1:
        generation_attempts = args.best_of_N

    for attempt_num in range(1, generation_attempts + 1):
        logger.info(f"Generation attempt {attempt_num}/{generation_attempts}.")
        current_pdb_content = None
        current_violations = []

        try:
            current_pdb_content = generate_pdb_content(
                length=length_for_generator,
                full_atom=args.full_atom,
                sequence_str=args.sequence,
                use_plausible_frequencies=args.plausible_frequencies,
            )

            if not current_pdb_content:
                logger.warning(f"Failed to generate PDB content in attempt {attempt_num}. Skipping.")
                continue

            if args.validate:
                logger.info("Performing PDB validation checks for current generation...")
                validator = PDBValidator(current_pdb_content)
                validator.validate_bond_lengths()
                validator.validate_bond_angles()
                validator.validate_ramachandran()
                validator.validate_steric_clashes()
                validator.validate_peptide_plane()
                validator.validate_sequence_improbabilities()
                current_violations = validator.get_violations()
            
            if args.guarantee_valid:
                if not current_violations:
                    logger.info(f"Successfully generated a valid PDB file after {attempt_num} attempts.")
                    final_pdb_content = current_pdb_content
                    final_violations = current_violations
                    break # Exit loop, valid PDB found
                else:
                    logger.warning(f"PDB generated in attempt {attempt_num} has {len(current_violations)} violations. Retrying...")
                    if logger.isEnabledFor(logging.DEBUG):
                        logger.debug("--- PDB Validation Report for failed attempt ---")
                        for violation in current_violations:
                            logger.debug(violation)
                        logger.debug("--- End Validation Report ---")
            elif args.best_of_N > 1:
                if len(current_violations) < min_violations_count:
                    min_violations_count = len(current_violations)
                    final_pdb_content = current_pdb_content
                    final_violations = current_violations
                    logger.info(f"Attempt {attempt_num} yielded {len(current_violations)} violations (new minimum).")
                else:
                    logger.info(f"Attempt {attempt_num} yielded {len(current_violations)} violations. Current minimum is {min_violations_count}.")
            else: # No guarantee-valid or best-of-N, just take the first one
                final_pdb_content = current_pdb_content
                final_violations = current_violations
                break

        except ValueError as ve:
            logger.error(f"Error processing sequence during generation (attempt {attempt_num}): {ve}. Skipping.")
        except Exception as e:
            logger.error(f"An unexpected error occurred during generation (attempt {attempt_num}): {e}. Skipping.")

    if final_pdb_content is None:
        logger.error(f"Failed to generate a suitable PDB file after {generation_attempts} attempts.")
        sys.exit(1)

    # After successful generation (and optional validation)
    if args.output:
        output_filename = args.output
        logger.debug("Using user-provided output filename: %s", output_filename)
    else:
        timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        if args.sequence:
            # Use a simplified sequence string for filename to avoid very long names
            sequence_tag = args.sequence.replace("-", "")[
                :10
            ]  # Take first 10 chars, remove hyphens
            output_filename = f"custom_peptide_{sequence_tag}_{timestamp}.pdb"
        else:
            output_filename = f"random_linear_peptide_{args.length}_{timestamp}.pdb"
        logger.debug("Generated default output filename: %s", output_filename)

    try:
        with open(output_filename, "w") as f:
            f.write(final_pdb_content)
        logger.info(
            "Successfully generated PDB file: %s", os.path.abspath(output_filename)
        )

        # Print final validation report
        if final_violations:
            logger.warning(f"--- PDB Validation Report for {os.path.abspath(output_filename)} ---")
            logger.warning(f"Final PDB has {len(final_violations)} violations.")
            for violation in final_violations:
                logger.warning(violation)
            logger.warning("--- End Validation Report ---")
        elif args.validate:
            logger.info(f"No violations found in the final PDB for {os.path.abspath(output_filename)}.")

    except Exception as e:
        logger.error("An unexpected error occurred during file writing: %s", e)
        sys.exit(1)


if __name__ == "__main__":
    main()
