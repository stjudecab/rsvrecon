#!/usr/bin/env python3

"""
Gene Sequence Extraction Tool

This script extracts coding sequences (CDS) from assembled genome sequences
using coordinates from GFF annotation files.

Author: Haidong Yi (hyi@stjude.org)
        Lei Li (lei.li@stjude.org)
Date: March 31, 2025
"""

import os
import re
import argparse
import logging
import sys
from typing import Tuple, Dict, Optional, List

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def setup_logger(log_file=None):
    """
    Configure and return a logger for the application.

    Args:
        log_file: Optional path to write logs to a file in addition to stdout

    Returns:
        Configured logger instance
    """
    handlers = [logging.StreamHandler(sys.stdout)]
    if log_file:
        handlers.append(logging.FileHandler(log_file))

    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
        handlers=handlers
    )
    return logging.getLogger('extract_gene_sequence')


def parse_gff_for_cds_coordinates(gff_file: str) -> Dict[str, Tuple[str, int, int]]:
    """
    Parse a GFF file to extract CDS coordinates.

    Args:
        gff_file: Path to the GFF annotation file

    Returns:
        Dictionary mapping CDS IDs to their coordinates (sequence_id, start, end)

    Raises:
        FileNotFoundError: If GFF file doesn't exist
        ValueError: If GFF file is improperly formatted
    """
    if not os.path.exists(gff_file):
        raise FileNotFoundError(f"GFF file not found: {gff_file}")

    gene_positions = {}
    with open(gff_file, 'r') as file:
        for line_num, line in enumerate(file, 1):
            if line.startswith('#'):
                continue

            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue

            if fields[2] == 'CDS':
                try:
                    seq_id = fields[0]
                    start = int(fields[3])
                    end = int(fields[4])

                    # Extract CDS ID from attributes field
                    attributes = fields[8]
                    id_match = re.search(r'ID=([^;]+)', attributes)

                    if id_match:
                        gene_name = id_match.group(1)
                    else:
                        # Fallback to last attribute if ID pattern not found
                        attr_parts = attributes.split(';')
                        if attr_parts and '=' in attr_parts[-1]:
                            gene_name = attr_parts[-1].split('=')[1]
                        else:
                            logger.warning(f"Could not extract CDS ID at line {line_num}, skipping")
                            continue

                    gene_positions[gene_name] = (seq_id, start, end)

                except (IndexError, ValueError) as e:
                    logger.warning(f"Error parsing line {line_num}: {e}")

    if not gene_positions:
        raise ValueError(f"No CDS features found in {gff_file}")

    return gene_positions


def extract_gene_sequence(
    assembled_sequence_file: str,
    reference_sequence_file: str,
    gff_file: str,
    cds_name: str
) -> Tuple[str, str]:
    """
    Extract a specific gene sequence from an assembled genome using coordinates from a GFF file.

    Args:
        assembled_sequence_file: Path to the assembled genome FASTA file
        reference_sequence_file: Path to the reference genome FASTA file
        gff_file: Path to the GFF annotation file
        cds_name: Name of the CDS to extract

    Returns:
        Tuple containing (sequence_id, extracted_sequence)

    Raises:
        FileNotFoundError: If input files don't exist
        ValueError: If CDS not found or sequences don't match expected format
    """
    # Validate input files
    for file_path, file_type in [
        (assembled_sequence_file, "assembled genome"),
        (reference_sequence_file, "reference sequence"),
        (gff_file, "GFF annotation")
    ]:
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"{file_type.capitalize()} file not found: {file_path}")

    # Load the assembled genome sequence
    try:
        assembled_genome_record = next(SeqIO.parse(assembled_sequence_file, "fasta"))
        assembled_genome_sequence = str(assembled_genome_record.seq)
        assembled_genome_name = str(assembled_genome_record.id)
    except StopIteration:
        raise ValueError(f"No sequences found in {assembled_sequence_file}")
    except Exception as e:
        raise ValueError(f"Error reading assembled genome: {e}")

    # Load the reference sequence
    try:
        reference = next(SeqIO.parse(reference_sequence_file, "fasta"))
        reference_sequence = str(reference.seq)
        reference_name = str(reference.id)
    except StopIteration:
        raise ValueError(f"No sequences found in {reference_sequence_file}")
    except Exception as e:
        raise ValueError(f"Error reading reference sequence: {e}")

    logger.info(f"Loaded assembled genome: {assembled_genome_name} ({len(assembled_genome_sequence)} bp)")
    logger.info(f"Loaded reference sequence: {reference_name} ({len(reference_sequence)} bp)")

    # Parse GFF to get gene coordinates
    gene_coordinates = parse_gff_for_cds_coordinates(gff_file)

    if cds_name not in gene_coordinates:
        raise ValueError(f"CDS '{cds_name}' not found in GFF file")

    seq_id, start, end = gene_coordinates[cds_name]
    logger.info(f"Found {cds_name} coordinates: {start}-{end} on {seq_id}")

    # Validate the sequence ID matches our assembled genome
    if seq_id != assembled_genome_name:
        logger.warning(
            f"Sequence ID in GFF ({seq_id}) doesn't match assembled genome ID ({assembled_genome_name})"
        )

    # Extract the sequence
    if 1 <= start <= end <= len(assembled_genome_sequence):
        extracted_sequence = assembled_genome_sequence[start-1:end]
        logger.info(f"Extracted sequence of length {len(extracted_sequence)} bp")
        return assembled_genome_name, extracted_sequence
    else:
        raise ValueError(
            f"Invalid coordinates for {cds_name}: {start}-{end} (genome length: {len(assembled_genome_sequence)})"
        )


def write_fasta(sequence_id: str, sequence: str, output_file: str):
    """
    Write a sequence to a FASTA file.

    Args:
        sequence_id: ID for the sequence
        sequence: The nucleotide sequence
        output_file: Path to the output file
    """
    # Create directory if it doesn't exist
    output_dir = os.path.dirname(output_file)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)

    with open(output_file, 'w') as fasta:
        fasta.write(f">{sequence_id}\n{sequence}\n")
    logger.info(f"Sequence written to {output_file}")


def process_multiple_cds(
    assembled_sequence_file: str,
    reference_sequence_file: str,
    gff_file: str,
    cds_list: List[str],
    output_prefix: str,
    output_dir: str = "."
) -> List[str]:
    """
    Process multiple CDS features and write each to a separate FASTA file.

    Args:
        assembled_sequence_file: Path to the assembled genome FASTA file
        reference_sequence_file: Path to the reference genome FASTA file
        gff_file: Path to the GFF annotation file
        cds_list: List of CDS names to extract
        output_prefix: Prefix for output filenames
        output_dir: Directory for output files

    Returns:
        List of paths to created output files
    """
    output_files = []

    # Create output directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for cds_name in cds_list:
        try:
            sequence_id, gene_sequence = extract_gene_sequence(
                assembled_sequence_file,
                reference_sequence_file,
                gff_file,
                cds_name
            )

            output_file = os.path.join(output_dir, f"{output_prefix}_{cds_name}.fasta")
            write_fasta(sequence_id, gene_sequence, output_file)
            output_files.append(output_file)

        except Exception as e:
            logger.error(f"Error processing {cds_name}: {e}")

    return output_files


def main():
    """Main function to parse arguments and run the sequence extraction."""
    parser = argparse.ArgumentParser(
        description="Extract coding sequences from assembled genomes using GFF annotations",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
        "-a", "--assembled",
        required=True,
        help="Path to the assembled genome FASTA file"
    )
    parser.add_argument(
        "-r", "--reference",
        required=True,
        help="Path to the reference genome FASTA file"
    )
    parser.add_argument(
        "-g", "--gff",
        required=True,
        help="Path to the GFF annotation file"
    )
    parser.add_argument(
        "-c", "--cds",
        default="CDS_7",
        help="Name of the CDS feature to extract (comma-separated for multiple)"
    )
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Path to the output FASTA file or directory prefix for multiple CDS"
    )
    parser.add_argument(
        "--outdir",
        default=".",
        help="Output directory for files (will be created if it doesn't exist)"
    )
    parser.add_argument(
        "--log",
        help="Log file path"
    )
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Enable verbose logging"
    )

    # Parse arguments
    args = parser.parse_args()

    # Setup logging
    logger = setup_logger(args.log)

    # Set logging level based on verbosity
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    logger.info("Starting gene sequence extraction")
    logger.info(f"Input assembled genome: {args.assembled}")
    logger.info(f"Input reference: {args.reference}")
    logger.info(f"Input GFF: {args.gff}")

    try:
        # Check if multiple CDS values were provided
        cds_list = [x.strip() for x in args.cds.split(',')]

        if len(cds_list) > 1:
            # Process multiple CDS
            logger.info(f"Processing {len(cds_list)} CDS features: {', '.join(cds_list)}")
            output_files = process_multiple_cds(
                args.assembled,
                args.reference,
                args.gff,
                cds_list,
                os.path.basename(args.output).replace('.fasta', ''),
                args.outdir
            )

            logger.info(f"Successfully extracted {len(output_files)} of {len(cds_list)} sequences")

            # Generate a summary file
            summary_file = os.path.join(args.outdir, f"{os.path.basename(args.output)}_summary.tsv")
            with open(summary_file, 'w') as f:
                f.write("cds_name\toutput_file\tstatus\n")
                for cds in cds_list:
                    expected_file = os.path.join(args.outdir, f"{os.path.basename(args.output).replace('.fasta', '')}_{cds}.fasta")
                    status = "SUCCESS" if os.path.exists(expected_file) else "FAILED"
                    f.write(f"{cds}\t{expected_file}\t{status}\n")

            logger.info(f"Summary written to {summary_file}")

        else:
            # Single CDS extraction
            cds_name = cds_list[0]
            logger.info(f"Processing single CDS: {cds_name}")

            sequence_id, gene_sequence = extract_gene_sequence(
                args.assembled,
                args.reference,
                args.gff,
                cds_name
            )

            # Write to FASTA
            write_fasta(sequence_id, gene_sequence, args.output)

            logger.info(f"Successfully extracted {cds_name} sequence to {args.output}")

        return 0

    except Exception as e:
        logger.error(f"Error: {e}")
        return 1


if __name__ == '__main__':
    # Run the main function
    exit_code = main()

    # Log completion status
    if exit_code == 0:
        logger.info("Gene extraction completed successfully")
    else:
        logger.error("Gene extraction failed")

    sys.exit(exit_code)

