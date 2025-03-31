#!/usr/bin/env python3
"""
G Gene RSV Genotyper

This script determines the genotype of Respiratory Syncytial Virus (RSV) based on
G gene sequence analysis using BLAST results. It takes pre-generated BLAST output
as input and identifies the closest matching reference sequence to determine the
G gene clade.

Author: Haidong Yi (hyi@stjude.org)
        Lei Li (lei.li@stjude.org)
Date: March 31, 2025
"""


import argparse
import pandas as pd
import os
import sys
import logging


def setup_logging():
    """
    Configure logging format and level.

    Args:
        log_file: Optional path to write logs to a file in addition to stdout

    Returns:
        Configured logger instance
    """
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s [%(levelname)s] %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    return logging.getLogger('g_rsv_genotyper')


def determine_g_genotype(blast_out_file, logger):
    """
    Determine G gene genotype from blastn results.

    Args:
        blast_out_file (str): Path to the blastn results file.
        logger (logging.Logger): Logger instance.

    Returns:
        tuple: (clade, blast_pct_identity, blast_alignment_length, subject_name)
    """
    logger.info(f"Reading G gene blast results from: {blast_out_file}")

    try:
        # Check if blast file exists and is not empty
        if not os.path.exists(blast_out_file) or os.path.getsize(blast_out_file) == 0:
            logger.warning("G gene blast file is empty or doesn't exist")
            return "Low G coverage", 0, 0, "unassigned"

        # Parse blast results
        blast_df = pd.read_csv(blast_out_file, sep='\t', header=None)
        if blast_df.empty:
            logger.warning("No G gene blast results found in file")
            return "Low G coverage", 0, 0, "unassigned"

        blast_df.columns = ['Query', 'Subject', 'Pct_identity', 'Alignment_length',
                           'Number_of_mismatches', 'Number_of_gap', 'Start_in_query',
                           'End_in_query', 'Start_in_subject', 'End_in_subject',
                           'E_value', 'Bit_score']
        blast_df = blast_df.sort_values(by=['Bit_score'], ascending=False)
        logger.info(f"Found {len(blast_df)} G gene blast matches")

        # Find the best match
        if len(blast_df) > 0:
            subject_name = blast_df.loc[0, 'Subject']
            blast_pct_identity = blast_df.loc[0, 'Pct_identity']
            blast_alignment_length = blast_df.loc[0, 'Alignment_length']
            clade = subject_name.split('_')[0]
            logger.info(f"Best G gene match: {subject_name} (clade: {clade})")
            return clade, blast_pct_identity, blast_alignment_length, subject_name
        else:
            logger.warning("No G gene matches found")
            return "Low G coverage", 0, 0, "unassigned"

    except Exception as e:
        logger.error(f"Error processing G gene blast results: {e}")
        return "Low G coverage", 0, 0, "unassigned"


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Genotype RSV G gene using blastn results')
    parser.add_argument('--blast_out', required=True, help='Path to G gene blastn results file')
    parser.add_argument('--output', default='Genotype_G.txt', help='Output file path (default: Genotype_G.txt)')
    parser.add_argument('--log', help="Log file path")


    # Optional parameters (for compatibility with pipeline)
    parser.add_argument('--query', help='Path to query file (not used for processing)')
    parser.add_argument('--subtype', help='Virus subtype (not used for processing)')
    parser.add_argument('--render_script', help='Path to render script (not used for processing)')
    parser.add_argument('--seq_id', help='Sequence ID (not used for processing)')

    return parser.parse_args()


def main():
    """Main function to determine G gene genotype and write results."""
    args = parse_args()
    logger = setup_logging(args.log)

    logger.info("=== RSV G Gene Genotyper ===")
    logger.info(f"G gene blast results: {args.blast_out}")
    logger.info(f"Output file: {args.output}")

    # Determine G gene genotype
    clade, blast_pct_identity, blast_alignment_length, subject_name = determine_g_genotype(
        args.blast_out, logger
    )

    # Write results to output file
    with open(args.output, 'w') as geno_file:
        result_line = f"{clade},{blast_pct_identity},{blast_alignment_length},{subject_name}"
        geno_file.write(result_line)
        logger.info(f"Results written to {args.output}: {result_line}")

    if clade == "Low G coverage":
        logger.warning("No valid G gene match found, possibly low coverage or negative control")
    else:
        logger.info(f"Successfully genotyped G gene as clade: {clade}")


if __name__ == "__main__":
    sys.exit(main())
