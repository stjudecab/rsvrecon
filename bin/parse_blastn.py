#!/usr/bin/env python3

"""
Whole Genome RSV Genotyper

This script analyzes blastn results to determine RSV (Respiratory Syncytial Virus) genotypes.
It takes a blastn output file and a metadata file as inputs, identifies the best match,
and outputs the genotype information in a CSV format.

Author: Haidong Yi (hyi@stjude.org)
        Lei Li (lei.li@stjude.org)
Date: March 28, 2025
"""

import argparse
import pandas as pd
import os
import sys
import logging


def setup_logging(log_file=None):
    """
    Configure logging format and level.

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
        format='%(asctime)s [%(levelname)s] %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    return logging.getLogger('rsv_genotyper')


def determine_genotype(blast_out_file, meta_file_path, logger):
    """
    Determine genotype from blastn results.

    Args:
        blast_out_file (str): Path to the blastn results file.
        meta_file_path (str): Path to the metadata file.
        logger (logging.Logger): Logger instance.

    Returns:
        tuple: (clade, blast_pct_identity, blast_alignment_length, strain)
    """
    logger.info(f"Reading blast results from: {blast_out_file}")
    logger.info(f"Using metadata from: {meta_file_path}")

    try:
        # Check if blast file exists and is not empty
        if not os.path.exists(blast_out_file) or os.path.getsize(blast_out_file) == 0:
            logger.warning("Blast file is empty or doesn't exist")
            return "No RSV", 0, 0, "No RSV"

        # Parse blast results
        blast_df = pd.read_csv(blast_out_file, sep='\t', header=None)
        if blast_df.empty:
            logger.warning("No blast results found in file")
            return "No RSV", 0, 0, "No RSV"

        blast_df.columns = ['Query', 'Subject', 'Pct_identity', 'Alignment_length',
                           'Number_of_mismatches', 'Number_of_gap', 'Start_in_query',
                           'End_in_query', 'Start_in_subject', 'End_in_subject',
                           'E_value', 'Bit_score']
        blast_df = blast_df.sort_values(by=['Bit_score'], ascending=False)
        logger.info(f"Found {len(blast_df)} blast matches")

        # Check if metadata file exists
        if not os.path.exists(meta_file_path):
            logger.error(f"Metadata file not found: {meta_file_path}")
            sys.exit(1)

        # Load metadata
        meta_df = pd.read_csv(meta_file_path, sep='\t', index_col=0)
        logger.info(f"Loaded metadata with {len(meta_df)} reference entries")

        # Find the best match
        clade = ''
        strain = ''
        blast_pct_identity = 0
        blast_alignment_length = 0
        sel_row_index = 0

        while clade == '' and sel_row_index < len(blast_df):
            selected_ref_name = blast_df.loc[sel_row_index, 'Subject']
            blast_pct_identity = blast_df.loc[sel_row_index, 'Pct_identity']
            blast_alignment_length = blast_df.loc[sel_row_index, 'Alignment_length']

            if selected_ref_name in meta_df.index:
                clade = meta_df.loc[selected_ref_name, 'clade']
                strain = meta_df.loc[selected_ref_name, 'strain']
                logger.info(f"Found match: {selected_ref_name} (clade: {clade}, strain: {strain})")
            else:
                logger.warning(f"Reference {selected_ref_name} not found in metadata")

            sel_row_index += 1

        if clade == '':
            logger.warning("No matching clade found in metadata")
            return "No RSV", 0, 0, "No RSV"

        return clade, blast_pct_identity, blast_alignment_length, strain

    except Exception as e:
        logger.error(f"Error processing blast results: {e}")
        return "No RSV", 0, 0, "No RSV"


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Genotype RSV samples using blastn results')
    parser.add_argument('--blast_out', required=True, help='Path to blastn results file')
    parser.add_argument('--meta_file', required=True, help='Path to metadata file')
    parser.add_argument('--output', default='Genotype.txt', help='Output file path (default: Genotype.txt)')
    parser.add_argument('--log', help="Log file path")

    # Optional parameters (for compatibility with pipeline)
    parser.add_argument('--query', help='Path to query file (not used for processing)')
    parser.add_argument('--reference', help='Path to reference folder (not used for processing)')
    parser.add_argument('--subtype', help='Virus subtype (not used for processing)')
    parser.add_argument('--seq_id', help='Sequence ID (not used for processing)')

    return parser.parse_args()


def main():
    """Main function to determine genotype and write results."""
    args = parse_args()
    logger = setup_logging(args.log)

    logger.info("=== RSV Genotyper ===")
    logger.info(f"Blast results: {args.blast_out}")
    logger.info(f"Metadata file: {args.meta_file}")
    logger.info(f"Output file: {args.output}")

    # Determine genotype
    clade, blast_pct_identity, blast_alignment_length, strain = determine_genotype(
        args.blast_out, args.meta_file, logger
    )

    # Write results to output file
    with open(args.output, 'w') as geno_file:
        result_line = f"{clade},{blast_pct_identity},{blast_alignment_length},{strain}"
        geno_file.write(result_line)
        logger.info(f"Results written to {args.output}: {result_line}")

    if clade == "No RSV":
        logger.warning("No valid RSV match found, most likely not RSV sample or negative control")
    else:
        logger.info(f"Successfully genotyped as clade: {clade}, strain: {strain}")


if __name__ == "__main__":
    sys.exit(main())

