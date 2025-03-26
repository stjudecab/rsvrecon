#!/usr/bin/env python3
"""
RSV Reference Selector

This script analyzes KMA mapping results to select the best matching RSV reference.
It extracts the corresponding FASTA sequence and GFF annotations based on the
top matching reference ID.

Author: Haidong Yi (hyi@stjude.org)
        Lei Li (lei.li@stjude.org)
Date: March 26, 2025
"""

import os
import re
import sys
import argparse
import logging
from pathlib import Path
from typing import Tuple, Dict, Optional, List

import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S"
)
logger = logging.getLogger(__name__)


class ReferenceNotFoundError(Exception):
    """Exception raised when a reference sequence cannot be found."""
    pass


class KMAResultError(Exception):
    """Exception raised when KMA result file is invalid or empty."""
    pass


def parse_args(args=None) -> argparse.Namespace:
    """
    Parse command line arguments.

    Args:
        args: Command line arguments (default: None)

    Returns:
        Parsed arguments as a Namespace object
    """
    parser = argparse.ArgumentParser(
        description="Select the best RSV reference based on KMA mapping results",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
        "--input",
        required=True,
        type=Path,
        help="Path to the KMA results file (.res)"
    )
    parser.add_argument(
        "--ref-fasta",
        required=True,
        type=Path,
        help="Path to the reference FASTA file"
    )
    parser.add_argument(
        "--ref-gff",
        required=True,
        type=Path,
        help="Path to the reference GFF file"
    )
    parser.add_argument(
        "--ref-metadata",
        required=True,
        type=Path,
        help="Path to the reference metadata CSV file"
    )
    parser.add_argument(
        "--output",
        required=True,
        type=Path,
        help="Directory to save the selected reference"
    )
    parser.add_argument(
        "--identity-cutoff",
        type=float,
        default=90.0,
        help="Minimum template identity percentage to classify as RSV"
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Enable verbose logging"
    )

    parsed_args = parser.parse_args(args)

    # Set logging level based on verbosity
    if parsed_args.verbose:
        logger.setLevel(logging.DEBUG)

    return parsed_args


def read_kma_results(kma_file: Path) -> Tuple[str, float]:
    """
    Read KMA results file and return the top match.

    Args:
        kma_file: Path to KMA results file

    Returns:
        Tuple of (reference_id, template_identity)

    Raises:
        KMAResultError: If the KMA file is empty or cannot be parsed
    """
    try:
        if not kma_file.exists():
            raise FileNotFoundError(f"KMA results file not found: {kma_file}")

        kma_df = pd.read_csv(kma_file, sep='\t')

        if kma_df.empty:
            raise KMAResultError("KMA results file is empty")

        # Sort by Score in descending order
        kma_df = kma_df.sort_values('Score', ascending=False)

        # Get top match details
        best_match_id = kma_df.iloc[0, 0].strip()
        template_identity = float(kma_df.iloc[0, 4])

        logger.debug(f"Best match: {best_match_id} (Identity: {template_identity}%)")
        return best_match_id, template_identity

    except pd.errors.EmptyDataError:
        raise KMAResultError("KMA results file is empty or malformed")
    except (IndexError, KeyError) as e:
        raise KMAResultError(f"Invalid KMA results format: {e}")
    except Exception as e:
        raise KMAResultError(f"Error reading KMA results: {e}")


def save_reference_files(
    output_dir: Path,
    ref_id: str,
    fasta_record: SeqRecord,
    gff_records: List[str]
) -> Tuple[Path, Path]:
    """
    Save reference FASTA and GFF files to the output directory.

    Args:
        output_dir: Directory to save files
        ref_id: Reference ID
        fasta_record: BioPython SeqRecord object
        gff_records: List of GFF record lines

    Returns:
        Tuple of (fasta_path, gff_path)
    """
    # Create output directory if it doesn't exist
    output_dir.mkdir(parents=True, exist_ok=True)

    # Define output file paths
    fasta_path = output_dir / f"{ref_id}.fasta"
    gff_path = output_dir / f"{ref_id}.gff"

    # Write FASTA file
    SeqIO.write(fasta_record, fasta_path, "fasta")

    # Write GFF file
    with open(gff_path, 'w') as f:
        f.writelines(gff_records)

    logger.info(f"Saved reference files: {fasta_path}, {gff_path}")
    return fasta_path, gff_path


def get_subtype_from_metadata(ref_id: str, metadata_file: Path) -> str:
    """
    Determine the subtype of a reference from metadata file.

    Args:
        ref_id: Reference ID to look up
        metadata_file: Path to metadata CSV file

    Returns:
        Subtype string

    Raises:
        FileNotFoundError: If metadata file doesn't exist
        KeyError: If reference ID not found in metadata
    """
    try:
        if not metadata_file.exists():
            raise FileNotFoundError(f"Metadata file not found: {metadata_file}")

        info_df = pd.read_csv(metadata_file, index_col=0)

        if ref_id not in info_df.index:
            raise KeyError(f"Reference ID {ref_id} not found in metadata")

        return info_df.loc[ref_id, 'Subtype']

    except pd.errors.EmptyDataError:
        raise ValueError("Metadata file is empty or malformed")


def extract_fasta_record(ref_id: str, fasta_file: Path) -> SeqRecord:
    """
    Extract a sequence record from a FASTA file by ID.

    Args:
        ref_id: Reference ID to extract
        fasta_file: Path to FASTA file

    Returns:
        BioPython SeqRecord object

    Raises:
        ReferenceNotFoundError: If reference ID not found in FASTA file
    """
    if not fasta_file.exists():
        raise FileNotFoundError(f"FASTA file not found: {fasta_file}")

    for record in SeqIO.parse(fasta_file, "fasta"):
        # Extract the ID from the FASTA header
        record_id = record.id.split()[0]
        if record_id == ref_id:
            return record

    raise ReferenceNotFoundError(f"Reference ID {ref_id} not found in FASTA file")


def extract_gff_records(ref_id: str, gff_file: Path) -> List[str]:
    """
    Extract GFF records for a specific reference ID.

    Args:
        ref_id: Reference ID to extract
        gff_file: Path to GFF file

    Returns:
        List of GFF lines for the reference

    Raises:
        ReferenceNotFoundError: If reference ID not found in GFF file
    """
    if not gff_file.exists():
        raise FileNotFoundError(f"GFF file not found: {gff_file}")

    gff_lines = []
    found = 0

    with open(gff_file, 'r') as f:
        for line in f:
            # Skip comments except sequence-region for our reference
            if line.startswith('#'):
                if line.startswith(f"##sequence-region\t{ref_id}") \
                    or line.startswith(f"##sequence-region {ref_id}"):
                    gff_lines.append(line)
                    found = True
                continue

            # Check if line contains our reference ID
            if line.split('\t')[0] == ref_id:
                gff_lines.append(line)
                found = True

    if not found:
        raise ReferenceNotFoundError(f"Reference ID {ref_id} not found in GFF file")

    # Ensure the GFF file has proper header and footer
    if not any(line.startswith("##gff-version") for line in gff_lines):
        gff_lines.insert(0, "##gff-version 3\n")

    if not any(line.startswith("###") for line in gff_lines):
        gff_lines.append("###\n")

    return gff_lines


def determine_rsv_status(
    reference_id: str,
    template_identity: float,
    subtype: str,
    identity_cutoff: float
) -> str:
    """
    Determine RSV status based on identity match.

    Args:
        reference_id: The best match reference ID
        template_identity: The identity percentage from KMA
        subtype: The RSV subtype from metadata
        identity_cutoff: Minimum identity percentage to classify as RSV

    Returns:
        Comma-separated string with classification, reference ID, and subtype
    """
    if template_identity > identity_cutoff:
        return f"{subtype},{reference_id},{subtype}"
    else:
        return f"Not RSV,{reference_id},{subtype}"


def main(args: argparse.Namespace) -> int:
    """
    Main function to process KMA results and select the best reference.

    Args:
        args: Parsed command line arguments

    Returns:
        Exit code (0 for success, non-zero for failure)
    """
    try:
        # Extract values from args
        kma_file = args.input
        fasta_file = args.ref_fasta
        gff_file = args.ref_gff
        metadata_file = args.ref_metadata
        output_dir = args.output
        identity_cutoff = args.identity_cutoff

        # Check that all required files exist
        for file_path, file_desc in [
            (kma_file, "KMA results"),
            (fasta_file, "FASTA reference"),
            (gff_file, "GFF reference"),
            (metadata_file, "Metadata CSV")
        ]:
            if not file_path.exists():
                logger.error(f"Required {file_desc} file not found: {file_path}")
                return 1

        # Read KMA results to get best matching reference
        best_ref_id, template_identity = read_kma_results(kma_file)

        # Clean up reference ID (remove any whitespace)
        best_ref_id = re.sub(r'\s', '', best_ref_id)

        # Get subtype information
        try:
            subtype = get_subtype_from_metadata(best_ref_id, metadata_file)
        except (FileNotFoundError, KeyError) as e:
            logger.warning(f"Could not determine subtype: {e}")
            subtype = "Unknown"

        # Extract reference sequences
        try:
            fasta_record = extract_fasta_record(best_ref_id, fasta_file)
            gff_records = extract_gff_records(best_ref_id, gff_file)
        except ReferenceNotFoundError as e:
            logger.error(str(e))
            return 1

        # Save to output directory
        save_reference_files(output_dir, best_ref_id, fasta_record, gff_records)

        # Determine and print RSV status
        rsv_status = determine_rsv_status(
            best_ref_id,
            template_identity,
            subtype,
            identity_cutoff
        )
        print(f"{best_ref_id},{rsv_status}")

        return 0

    except KMAResultError as e:
        logger.error(f"KMA result error: {e}")
        return 1
    except FileNotFoundError as e:
        logger.error(f"File not found: {e}")
        return 1
    except Exception as e:
        logger.error(f"Unexpected error: {e}")
        logger.debug("Exception details:", exc_info=True)
        return 1


# main entrypoint
if __name__ == "__main__":
    sys.exit(main(parse_args()))
