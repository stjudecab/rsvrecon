#!/usr/bin/env python3
"""
IGV Processing Script for Nextflow

This script processes IGV WIG files to generate consensus sequences
based on nucleotide coverage data and a reference genome.

Author: Haidong Yi (hyi@stjude.org)
        Lei Li (lei.li@stjude.org)
Date: March 27, 2025
"""

import os
import re
import sys
import logging
import argparse
from pathlib import Path
from typing import Dict, List, Tuple
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger("IGVProcessor")


class IGVProcessor:
    """A class for processing IGV WIG files and generating consensus sequences."""

    def __init__(self, wig_file: str, ref_file: str, cutoff: int):
        """
        Initialize the IGV processor.

        Args:
            wig_file: Path to the WIG file
            ref_file: Path to the reference FASTA file
            cutoff: Coverage cutoff for base calling
        """
        self.wig_file = Path(wig_file)
        self.ref_file = Path(ref_file)
        self.cutoff = cutoff

        # Validate input files
        self._validate_inputs()

        # Load reference genome
        self.ref_genome_record, self.ref_genome_sequence = self._load_reference()
        self.ref_length = len(self.ref_genome_sequence)
        logger.info(f"Loaded reference genome: {self.ref_genome_record.id} ({self.ref_length} bp)")

    def _validate_inputs(self) -> None:
        """Validate that input files exist and are readable."""
        if not self.wig_file.exists():
            logger.error(f"WIG file not found: {self.wig_file}")
            sys.exit(1)

        if not self.ref_file.exists():
            logger.error(f"Reference file not found: {self.ref_file}")
            sys.exit(1)

        if not os.access(self.wig_file, os.R_OK):
            logger.error(f"No read permission for WIG file: {self.wig_file}")
            sys.exit(1)

        if not os.access(self.ref_file, os.R_OK):
            logger.error(f"No read permission for reference file: {self.ref_file}")
            sys.exit(1)

    def _load_reference(self) -> Tuple[SeqRecord, str]:
        """
        Load the reference genome sequence.

        Returns:
            Tuple containing the reference genome record and sequence as a string
        """
        try:
            ref_genome_record = next(SeqIO.parse(self.ref_file, "fasta"))
            ref_genome_sequence = str(ref_genome_record.seq)
            return ref_genome_record, ref_genome_sequence
        except StopIteration:
            logger.error(f"No sequence found in reference file: {self.ref_file}")
            sys.exit(1)
        except Exception as e:
            logger.error(f"Error reading reference file: {str(e)}")
            sys.exit(1)

    def process(self) -> str:
        """
        Process the WIG file and generate a consensus sequence.

        Returns:
            The consensus sequence as a string
        """
        # Initialize sequence array with Ns
        sequence_array = ["N"] * self.ref_length

        try:
            with open(self.wig_file, 'r') as file:
                line_count = 0
                base_called_count = 0

                for line_num, line in enumerate(file, 1):
                    line = line.strip()

                    # Skip header lines
                    if line.startswith('track') or line.startswith('#') or line.startswith('variableStep chrom='):
                        continue

                    # Process data lines
                    elif line[0].isdigit():
                        line_count += 1
                        match = re.match(r'^(\d+)\t(\d+)\.0\t(\d+)\.0\t(\d+)\.0\t(\d+)\.0.+', line)

                        if match:
                            pos, a, c, g, t = map(int, match.groups())
                            pos = pos - 1  # Convert to 0-based indexing

                            if pos >= self.ref_length:
                                logger.warning(f"Position {pos+1} exceeds reference length at line {line_num}, skipping")
                                continue

                            count_dict = {'A': a, 'C': c, 'G': g, 'T': t}
                            cov = a + c + g + t

                            if cov > self.cutoff:
                                max_base = max(count_dict, key=count_dict.get)
                                sequence_array[pos] = max_base
                                base_called_count += 1
                        else:
                            logger.warning(f"Malformed data line at line {line_num}: {line}")
                    else:
                        logger.warning(f"Unrecognized line format at line {line_num}: {line}")

            # Create the consensus sequence
            sequence = ''.join(sequence_array)

            # Log processing results
            logger.info(f"Processed {line_count} data lines from {self.wig_file}")
            logger.info(f"Called {base_called_count} bases ({base_called_count/self.ref_length:.2%} of reference)")

            return sequence

        except Exception as e:
            logger.error(f"Error processing WIG file: {str(e)}")
            sys.exit(1)

def write_fasta(sequence: str, output_file: str, sample_id: str) -> None:
    """
    Write the consensus sequence to a FASTA file.

    Args:
        sequence: The consensus sequence
        output_file: Path to the output file
        sample_id: Sample identifier for the FASTA header
    """
    try:
        output_path = Path(output_file)

        # Ensure the output directory exists
        output_path.parent.mkdir(parents=True, exist_ok=True)

        with open(output_path, 'w') as fasta:
            fasta.write(f'>{sample_id}\n{sequence}\n')

        logger.info(f"Consensus sequence written to {output_file}")
    except Exception as e:
        logger.error(f"Error writing output file: {str(e)}")
        sys.exit(1)

def parse_arguments() -> argparse.Namespace:
    """
    Parse command-line arguments.

    Returns:
        Parsed arguments
    """
    parser = argparse.ArgumentParser(
        description='Process IGV WIG files to generate consensus sequences',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
        '--wig',
        required=True,
        help='Path to the WIG file'
    )
    parser.add_argument(
        '--ref-fasta',
        required=True,
        help='Path to reference FASTA file'
    )
    parser.add_argument(
        '--igv-cutoff',
        type=int,
        default=10,
        help='Coverage cutoff for base calling'
    )
    parser.add_argument(
        '--out',
        required=True,
        help='Path to output FASTA file'
    )
    parser.add_argument(
        '--id',
        required=True,
        help='Sample identifier for FASTA header'
    )
    parser.add_argument(
        '--log-level',
        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'],
        default='INFO',
        help='Set the logging level'
    )

    return parser.parse_args()

def main() -> None:
    """Main entry point for the script."""
    # Parse command-line arguments
    args = parse_arguments()

    # Set log level
    logger.setLevel(getattr(logging, args.log_level))

    try:
        # Log script start
        logger.info(f"Starting IGV processing with cutoff: {args.igv_cutoff}")

        # Process IGV file
        processor = IGVProcessor(args.wig, args.ref_fasta, args.igv_cutoff)
        consensus_sequence = processor.process()

        # Write output
        write_fasta(consensus_sequence, args.out, args.id)

        # Log script completion
        logger.info("IGV processing completed successfully")

    except Exception as e:
        logger.error(f"Unexpected error: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    main()
