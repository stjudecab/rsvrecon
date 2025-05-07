#!/usr/bin/env python3
"""
RSV report files generation script

This script takes the output files from the `stjudecab/rsvrecon` pipeline for
generating the report files for analyzing Respiratory Syncytial Virus (RSV)
genomic data, including sequence assembly, mutation detection, and phylogenetic analysis.

Author: Haidong Yi (hyi@stjude.org)
        Lei Li (lei.li@stjude.org)
Date: May 6th, 2025
"""

import os
import re
import json
import glob
import shutil
import logging
import argparse
import pandas as pd
from io import StringIO
from typing import Dict, List, Tuple, Optional, Any
from Bio import SeqIO
from Bio.Seq import Seq


class RSVAnalysisPipeline:
    """Main class for the RSV analysis pipeline."""

    def __init__(self, meta_file: str, version_file: str, output_dir: str, reference_dir: str,
                 coverage_threshold: int, next_rsv_pipe_res: Optional[str] = None):
        """
        Initialize the RSV Analysis Pipeline.

        Args:
            meta_file: Path to the metadata file containing sample information
            version_file: Path to the version file containing pipeline version information
            output_dir: Directory to store output files
            reference_dir: Directory containing reference files
            coverage_threshold: Minimum coverage threshold for gene analysis
            next_rsv_pipe_res: Optional path to NEXT-RSV-PIPE results for comparison
        """
        self.meta_file = meta_file
        self.version_file = version_file
        self.output_dir = output_dir
        self.reference_dir = reference_dir
        self.coverage_threshold = coverage_threshold
        self.next_rsv_pipe_res = next_rsv_pipe_res
        self.version = self._get_version()

        # Create output directory if it doesn't exist
        os.makedirs(self.output_dir, exist_ok=True)

        # Configure logging
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s [%(levelname)s] %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S',
            handlers=[
                logging.FileHandler(os.path.join(self.output_dir, "pipeline.log")),
                logging.StreamHandler()
            ]
        )
        self.logger = logging.getLogger(__name__)

    def _get_version(self) -> str:
        """
        Get the current version of the pipeline from version.txt.

        Returns:
            The version string
        """
        try:
            with open(self.version_file, "r") as file:
                return file.readline().strip()
        except FileNotFoundError:
            self.logger.warning("Version file not found. Using 'unknown' as version.")
            return "unknown"

    def run(self) -> None:
        """Execute the full analysis pipeline."""
        self.logger.info(f"Starting RSV Report of pipeline version: {self.version}")
        self.logger.info(f"Processing metadata file: {self.meta_file}")

        # Initialize output files and data structures
        output_files = self._initialize_output_files()
        subtype_a_samples = []
        subtype_b_samples = []
        f_protein_mutations_a = {}
        f_protein_mutations_b = {}

        # Read metadata file
        meta_df = pd.read_csv(self.meta_file, sep='\t', index_col=0)

        # Write CSV header
        with open(output_files['report'], 'w') as report_file:
            report_file.write(self._generate_report_header())

        # Process each sample
        for sample_name in sorted(meta_df.index.values):
            self.logger.info(f"Processing sample: {sample_name}")

            try:
                # Process the sample and get results
                results = self._process_sample(sample_name, meta_df, output_files)

                if results:
                    # Update subtype lists and F protein mutations
                    if results['subtype'].startswith('A'):
                        subtype_a_samples.append(sample_name)
                        if 'f_mutations' in results:
                            f_protein_mutations_a[sample_name] = results['f_mutations']
                    elif results['subtype'].startswith('B'):
                        subtype_b_samples.append(sample_name)
                        if 'f_mutations' in results:
                            f_protein_mutations_b[sample_name] = results['f_mutations']
            except Exception as e:
                self.logger.error(f"Error processing sample {sample_name}: {str(e)}")

        # Generate mutation reports and phylogenetic tree files
        self._generate_f_mutation_reports(f_protein_mutations_a, f_protein_mutations_b, output_files)

        if subtype_a_samples:
            self._prepare_phylogenetic_tree_files('A', subtype_a_samples, output_files)

        if subtype_b_samples:
            self._prepare_phylogenetic_tree_files('B', subtype_b_samples, output_files)

        self.logger.info("Analysis complete! Report and sequences have been generated.")

    def _initialize_output_files(self) -> Dict[str, str]:
        """
        Initialize output file paths.

        Returns:
            Dictionary containing all output file paths
        """
        output_files = {
            'report': os.path.join(self.output_dir, "Report.csv"),
            'sequence': os.path.join(self.output_dir, "Sequence.fasta"),
            'sequence_error': os.path.join(self.output_dir, "Sequence_err.fasta"),
            'sequence_a': os.path.join(self.output_dir, "Sequence_A.fasta"),
            'sequence_b': os.path.join(self.output_dir, "Sequence_B.fasta"),
            'f_protein_a': os.path.join(self.output_dir, 'F_protein_seq_A.fasta'),
            'f_protein_b': os.path.join(self.output_dir, 'F_protein_seq_B.fasta'),
            'f_report_a': os.path.join(self.output_dir, 'F_mutation_A_report.csv'),
            'f_report_b': os.path.join(self.output_dir, 'F_mutation_B_report.csv'),
            'f_key_position_a': os.path.join(self.output_dir, 'F_protein_key_position_A.csv'),
            'f_key_position_b': os.path.join(self.output_dir, 'F_protein_key_position_B.csv'),
            'tree_seq_a': os.path.join(self.output_dir, 'tree_sequences_A.fasta'),
            'tree_anno_a': os.path.join(self.output_dir, 'tree_sequences_A.csv'),
            'tree_seq_b': os.path.join(self.output_dir, 'tree_sequences_B.fasta'),
            'tree_anno_b': os.path.join(self.output_dir, 'tree_sequences_B.csv')
        }

        return output_files

    def _generate_report_header(self) -> str:
        """
        Generate the report CSV header.

        Returns:
            Header string for the CSV report
        """
        header = [
            f"Pipeline version: {self.version}",
            "Sample name,before_filtering_total_reads,before_filtering_q20_rate,before_filtering_q30_rate," +
            "after_filtering_total_reads,after_filtering_q20_rate,after_filtering_q30_rate,QC rate," +
            "Uniquely mapped reads %,MULTI-MAPPING READS %,UNMAPPED READS%,CHIMERIC READS%," +
            "Subtype,reference_accession,ref_subtype," +
            "F protein mutations," +
            "NS1_cov,NS2_cov,N_cov,P_cov,M_cov,SH_cov,G_cov,F_cov,M2-1_cov,M2-2_cov,L_cov," +
            "Whole Genome Clade(NextClade),Whole Genome Clade(Blast)"
        ]

        return "\n".join(header) + "\n"

    def _process_sample(self, sample_name: str, meta_df: pd.DataFrame,
                        output_files: Dict[str, str]) -> Optional[Dict[str, Any]]:
        """
        Process a single sample through the entire analysis pipeline.

        Args:
            sample_name: Name of the sample to process
            meta_df: Metadata DataFrame containing sample information
            output_files: Dictionary of output file paths

        Returns:
            Dictionary containing sample analysis results or None if processing failed
        """
        # Get QC information
        qc_info = self._parse_qc_data(meta_df.loc[sample_name, 'qc_fastp'])
        if not qc_info:
            return None

        # Get mapping information
        flagstat_file = meta_df.loc[sample_name, 'flagstat']
        if not os.path.exists(flagstat_file):
            self._write_error_line(sample_name, qc_info, output_files['report'])
            return None

        # Parse flagstat data
        mapping_stats = self._parse_flagstat(flagstat_file)

        # Determine RSV subtype
        genotype_file = meta_df.loc[sample_name, 'whg_genotype']
        subtype_info = self._determine_subtype(genotype_file, meta_df.loc[sample_name, 'kma_out'])
        assert subtype_info['subtype'].startswith('A') or subtype_info['subtype'].startswith('B'), \
            f"{sample_name}'s subtype is not 'A' or 'B'."

        # Process genome sequence
        genome_sequence = self._extract_genome_sequence(meta_df.loc[sample_name, 'assembly'])
        if genome_sequence:
            self._write_genome_sequences(sample_name, genome_sequence, subtype_info['subtype'][0], output_files)

        # Detect F protein mutations if sample is RSV
        f_mutation_results = {}
        if subtype_info['subtype'] not in ['Not RSV', '']:
            # detect the F mutation per sample
            f_mutation_results = self._detect_f_mutations(
                sample_name,
                subtype_info['subtype'][0],
                meta_df.loc[sample_name, 'assembly'],
                meta_df.loc[sample_name, 'ref_fasta'],
                meta_df.loc[sample_name, 'ref_gff'],
                output_files
            )

        # Calculate gene coverage
        wig_file = meta_df.loc[sample_name, 'igv_out']
        gff_file = meta_df.loc[sample_name, 'ref_gff']
        gene_coverage = self._calculate_gene_coverage(wig_file, gff_file)

        # Get NextClade and BLAST genotype information
        nextclade_info = self._get_nextclade_info(meta_df.loc[sample_name, 'nextclade_out']
                                                  if subtype_info['subtype'] != 'Not RSV' else None)
        blast_info = self._get_blast_info(genotype_file if subtype_info['subtype'] != 'Not RSV' else None)

        # Write results to report file
        self._write_report_line(
            sample_name, qc_info, mapping_stats, subtype_info,
            f_mutation_results, gene_coverage, nextclade_info, blast_info,
            output_files['report']
        )

        # Return results for further processing
        return {
            'subtype': subtype_info['subtype'],
            'f_mutations': f_mutation_results.get('mutation_list', [])
        }

    def _parse_qc_data(self, qc_json_file: str) -> Optional[Dict[str, Any]]:
        """
        Parse quality control data from FastP JSON output.

        Args:
            qc_json_file: Path to the FastP JSON file

        Returns:
            Dictionary containing QC metrics or None if file cannot be read
        """
        try:
            with open(qc_json_file, 'r') as file:
                qc_data = json.load(file)

            before = qc_data['summary']['before_filtering']
            after = qc_data['summary']['after_filtering']
            qc_rate = (after['total_reads'] / before['total_reads']) * 100 if before['total_reads'] > 0 else 0

            return {
                'before_total': before['total_reads'],
                'before_q20': before['q20_rate'],
                'before_q30': before['q30_rate'],
                'after_total': after['total_reads'],
                'after_q20': after['q20_rate'],
                'after_q30': after['q30_rate'],
                'qc_rate': qc_rate
            }
        except (FileNotFoundError, json.JSONDecodeError, KeyError) as e:
            self.logger.error(f"Error parsing QC data: {str(e)}")
            return None

    def _parse_flagstat(self, flagstat_file: str) -> Dict[str, float]:
        """
        Parse alignment statistics from samtools flagstat output.

        Args:
            flagstat_file: Path to the flagstat file

        Returns:
            Dictionary containing mapping statistics
        """
        try:
            with open(flagstat_file, 'r') as file:
                lines = [line.strip() for line in file.readlines()]

            # Extract key metrics
            metrics = {'total_reads': 0, 'primary_mapped': 0, 'total_mapped': 0, 'supplementary': 0}

            for line in lines:
                if 'in total' in line:
                    metrics['total_reads'] = int(line.split('+')[0].strip())
                elif 'primary mapped' in line:
                    metrics['primary_mapped'] = int(line.split('+')[0].strip())
                elif 'mapped (' in line and 'primary' not in line:
                    metrics['total_mapped'] = int(line.split('+')[0].strip())
                elif 'supplementary' in line:
                    metrics['supplementary'] = int(line.split('+')[0].strip())

            # Calculate percentages
            total = metrics['total_reads'] if metrics['total_reads'] > 0 else 1  # Avoid division by zero

            return {
                'uniquely_mapped': (metrics['primary_mapped'] / total) * 100,
                'multi_mapped': ((metrics['total_mapped'] - metrics['primary_mapped']) / total) * 100,
                'unmapped': ((total - metrics['total_mapped']) / total) * 100,
                'chimeric': (metrics['supplementary'] / total) * 100
            }
        except (FileNotFoundError, ValueError, IndexError) as e:
            self.logger.error(f"Error parsing flagstat data: {str(e)}")
            return {'uniquely_mapped': 0, 'multi_mapped': 0, 'unmapped': 100, 'chimeric': 0}

    def _determine_subtype(self, genotype_file: str, kma_out_file: str) -> Dict[str, str]:
        """
        Determine the RSV subtype from genotyping data.

        Args:
            genotype_file: Path to the genotype file
            kma_out_file: Path to the KMA output file

        Returns:
            Dictionary containing subtype information
        """
        # Check BLAST genotype results
        blast_subtype = 'Not RSV'
        blast_alignment_length = 0

        if os.path.isfile(genotype_file):
            try:
                with open(genotype_file, 'r') as file:
                    genotype_data = file.readline().strip().split(',')
                    blast_subtype = genotype_data[0]
                    blast_alignment_length = int(genotype_data[2])

                # Require minimum alignment length
                if blast_alignment_length < 1500:
                    blast_subtype = 'Not RSV'
            except (FileNotFoundError, IndexError, ValueError) as e:
                self.logger.error(f"Error reading genotype file: {str(e)}")

        # Check KMA reference
        ref_subtype = 'Not RSV'
        ref_accession = ''

        try:
            kma_df = pd.read_csv(kma_out_file, sep='\t')
            if not kma_df.empty:
                kma_df = kma_df.sort_values('Score', ascending=False)
                ref_accession = kma_df.iloc[0, 0]
                template_identity = float(kma_df.iloc[0, 4])

                # Get reference subtype information
                ref_info_file = os.path.join(self.reference_dir, 'ReferenceCandidate', 'ReferenceCandidate_subtype.csv')
                ref_info = pd.read_csv(ref_info_file, delimiter=',', index_col=0)
                ref_subtype = ref_info.loc[ref_accession, 'Subtype']

                # Apply identity threshold
                if template_identity <= 10:  # Using the threshold from the original code
                    ref_subtype = 'Not RSV'
        except (FileNotFoundError, pd.errors.EmptyDataError, KeyError, IndexError) as e:
            self.logger.error(f"Error determining subtype from KMA: {str(e)}")

        # Final subtype determination (prioritize BLAST result)
        final_subtype = blast_subtype if blast_subtype != 'Not RSV' else ref_subtype

        return {
            'subtype': final_subtype,
            'ref_accession': ref_accession,
            'ref_subtype': ref_subtype
        }

    def _extract_genome_sequence(self, genome_file: str) -> Optional[str]:
        """
        Extract genome sequence from FASTA file.

        Args:
            genome_file: Path to the genome FASTA file

        Returns:
            Genome sequence string or None if file cannot be read
        """
        try:
            with open(genome_file, 'r') as file:
                sequence = ''
                for line in file:
                    if not line.startswith('>'):
                        sequence += line.strip()
                return sequence
        except FileNotFoundError:
            self.logger.error(f"Genome file not found: {genome_file}")
            return None

    @staticmethod
    def _write_genome_sequences(sample_name: str, sequence: str, subtype: str,
                                output_files: Dict[str, str]) -> None:
        """
        Write genome sequences to appropriate output files.

        Args:
            sample_name: Sample name
            sequence: Genome sequence
            subtype: RSV subtype
            output_files: Dictionary of output file paths
        """
        fasta_entry = f">{sample_name}\n{sequence}\n"

        # Write to appropriate files based on subtype
        if subtype == 'Not RSV':
            with open(output_files['sequence_error'], 'a') as file:
                file.write(fasta_entry)
        elif subtype == 'A':
            with open(output_files['sequence'], 'a') as file:
                file.write(fasta_entry)
            with open(output_files['sequence_a'], 'a') as file:
                file.write(fasta_entry)
        elif subtype == 'B':
            with open(output_files['sequence'], 'a') as file:
                file.write(fasta_entry)
            with open(output_files['sequence_b'], 'a') as file:
                file.write(fasta_entry)

    def _detect_f_mutations(self, sample_name: str, subtype: str,
                            genome_file: str, ref_file: str, gff_file: str,
                            output_files: Dict[str, str]) -> Dict[str, Any]:
        """
        Detect mutations in the F protein.

        Args:
            sample_name: Sample name
            subtype: RSV subtype
            genome_file: Path to the genome file
            ref_file: Path to the reference genome file
            gff_file: Path to the GFF file
            output_files: Dictionary of output file paths

        Returns:
            Dictionary containing F protein mutation information
        """
        results = {'mutation_str': '', 'mutation_list': []}

        # Skip if not RSV
        if subtype == 'Not RSV':
            return results

        # Select appropriate mutation file based on subtype
        mutation_file = os.path.join(
            self.reference_dir,
            'F_mutation',
            f"RSV_{subtype}_F_Mutation.csv"
        )

        # Extract F gene sequence
        sequence_data = self._extract_gene_sequence(
            genome_file, ref_file, gff_file, 'CDS_8'
        )

        if not sequence_data:
            return results

        # Translate to protein sequence
        aa_sequence = self._translate_to_protein(sequence_data[1])

        # Write protein sequence to file
        f_protein_file = output_files['f_protein_a'] if subtype == 'A' else output_files['f_protein_b']
        with open(f_protein_file, 'a') as file:
            file.write(f">{sample_name}\n{aa_sequence}\n")

        # Detect mutations
        mutations = self._detect_protein_mutations(aa_sequence, mutation_file)

        # Format mutation string
        mutation_strings = [f"{m[0]}({m[1]})" for m in mutations]

        return {
            'mutation_str': "|".join(mutation_strings),
            'mutation_list': mutation_strings
        }

    def _extract_gene_sequence(self, fasta_file: str, ref_file: str, gff_file: str,
                               gene_id: str) -> Optional[Tuple[str, str]]:
        """
        Extract a specific gene sequence from a genome.

        Args:
            fasta_file: Path to the genome FASTA file
            ref_file: Path to the reference genome FASTA file
            gff_file: Path to the GFF file
            gene_id: ID of the gene to extract

        Returns:
            Tuple of (sequence_name, gene_sequence) or None if extraction fails
        """
        try:
            # Load the assembled genome sequence
            assembled_record = next(SeqIO.parse(fasta_file, "fasta"))
            assembled_sequence = str(assembled_record.seq)
            assembled_name = str(assembled_record.id)

            # Load the reference
            reference = next(SeqIO.parse(ref_file, "fasta"))
            reference_sequence = str(reference.seq)

            # Make sure sequences are compatible
            if len(reference_sequence) != len(assembled_sequence):
                self.logger.warning(
                    f"Reference and assembled sequence length mismatch: {len(reference_sequence)} vs {len(assembled_sequence)}"
                )
                return None

            # Parse GFF to get gene coordinates
            gene_coordinates = self._parse_gff_for_cds(gff_file)

            if gene_id not in gene_coordinates:
                self.logger.error(f"Gene ID {gene_id} not found in GFF file")
                return None

            start = gene_coordinates[gene_id][1]
            end = gene_coordinates[gene_id][2]

            # Extract the gene sequence
            extracted_sequence = assembled_sequence[start - 1:end]

            return assembled_name, extracted_sequence
        except (FileNotFoundError, StopIteration, KeyError, IndexError) as e:
            self.logger.error(f"Error extracting gene sequence: {str(e)}")
            return None

    def _parse_gff_for_cds(self, gff_file: str) -> Dict[str, Tuple[str, int, int]]:
        """
        Parse a GFF file to extract CDS coordinates.

        Args:
            gff_file: Path to the GFF file

        Returns:
            Dictionary mapping gene IDs to (sequence_id, start, end) tuples
        """
        gene_positions = {}

        try:
            with open(gff_file, 'r') as file:
                for line in file:
                    if not line.startswith('#'):
                        fields = line.strip().split('\t')
                        if len(fields) >= 9 and fields[2] == 'CDS':
                            seq_id = fields[0]
                            start = int(fields[3])
                            end = int(fields[4])

                            # Extract gene name from attributes
                            attributes = fields[8]
                            gene_name = None

                            # Try to find ID in attributes
                            match = re.search(r'ID=([^;]+)', attributes)
                            if match:
                                gene_name = match.group(1)

                            # Fall back to last attribute
                            if not gene_name:
                                gene_name = attributes.split(';')[-1].split('=')[1]

                            gene_positions[gene_name] = (seq_id, start, end)

            return gene_positions
        except (FileNotFoundError, IndexError) as e:
            self.logger.error(f"Error parsing GFF file: {str(e)}")
            return {}

    @staticmethod
    def _translate_to_protein(nucleotide_sequence: str) -> str:
        """
        Translate a nucleotide sequence to amino acid sequence.

        Args:
            nucleotide_sequence: Nucleotide sequence string

        Returns:
            Amino acid sequence string
        """
        seq = Seq(nucleotide_sequence)
        return str(seq.translate())

    def _detect_protein_mutations(self, protein_sequence: str, mutation_file: str) -> List[List[str]]:
        """
        Detect mutations in a protein sequence compared to a reference.

        Args:
            protein_sequence: Protein sequence to analyze
            mutation_file: Path to the mutation definition file

        Returns:
            List of [mutation_name, mutation_type] lists
        """
        # Load mutation definitions
        single_mutations = []
        single_mutation_indices = []
        co_mutations = []

        try:
            with open(mutation_file, 'r') as file:
                for line in file:
                    mutation_def = line.strip()

                    if '+' in mutation_def:
                        # Co-mutations
                        co_mutations.append(mutation_def.split('+'))
                    else:
                        # Single mutations
                        parts = mutation_def.split(',')
                        if len(parts) == 3:
                            single_mutations.append([parts[0], parts[2]])
                            single_mutation_indices.append(parts[1])

            # Create DataFrame for single mutations
            single_mutation_df = pd.DataFrame(
                single_mutations,
                index=single_mutation_indices,
                columns=['original', 'alternative']
            )

            # Detect mutations
            detected_mutations = []

            # Check single mutations
            for position in single_mutation_df.index:
                pos_idx = int(position) - 1
                if pos_idx >= len(protein_sequence):
                    continue

                current_aa = protein_sequence[pos_idx]
                original_aa = single_mutation_df.loc[position]['original']
                alternative_aa = single_mutation_df.loc[position]['alternative']

                # Skip if matches reference
                if current_aa in original_aa:
                    continue

                # Skip ambiguous positions
                if current_aa == 'X':
                    continue

                # Check if it's a known mutation
                mutation_type = 'Novel'
                if current_aa in alternative_aa:
                    mutation_type = 'Reported'

                mutation_name = f"{original_aa}{position}{current_aa}"
                detected_mutations.append([mutation_name, mutation_type])

            # Check co-mutations
            for co_mutation in co_mutations:
                co_mutation_names = []
                all_present = True

                for mutation in co_mutation:
                    parts = mutation.split(',')
                    if len(parts) == 3:
                        position = int(parts[1])
                        expected_aa = parts[2]

                        if position - 1 >= len(protein_sequence):
                            all_present = False
                            break

                        current_aa = protein_sequence[position - 1]
                        if current_aa not in expected_aa:
                            all_present = False
                            break

                        co_mutation_names.append(f"{parts[0]}{parts[1]}{current_aa}")

                if all_present:
                    detected_mutations.append(['+'.join(co_mutation_names), 'Reported'])

            return detected_mutations

        except (FileNotFoundError, pd.errors.EmptyDataError, KeyError, IndexError) as e:
            self.logger.error(f"Error detecting protein mutations: {str(e)}")
            return []

    def _calculate_gene_coverage(self, wig_file: str, gff_file: str) -> Dict[str, float]:
        """
        Calculate coverage for each gene from an alignment WIG file.

        Args:
            wig_file: Path to the coverage WIG file
            gff_file: Path to the GFF file

        Returns:
            Dictionary mapping gene IDs to coverage values
        """
        try:
            # Parse gene coordinates from GFF
            gene_coordinates = self._parse_gff_for_cds(gff_file)

            # Read coverage data from WIG file
            coverage = self._read_wig_coverage(wig_file)

            if not coverage.any():
                return {gene: 0 for gene in gene_coordinates}

            index_pool = list(coverage.index)

            # Calculate coverage for each gene
            gene_coverage = {}
            for gene in gene_coordinates:
                start = gene_coordinates[gene][1]
                end = gene_coordinates[gene][2]
                gene_len = end - start + 1

                # Get positions within gene boundaries
                gene_positions = [x for x in index_pool if start <= x <= end]

                if gene_positions:
                    # Calculate average coverage
                    gene_coverage[gene] = sum(coverage[gene_positions].values) / gene_len
                else:
                    gene_coverage[gene] = 0

            return gene_coverage

        except Exception as e:
            self.logger.error(f"Error calculating gene coverage: {str(e)}")
            return {}

    def _read_wig_coverage(self, wig_file: str) -> Optional[pd.Series]:
        """
        Read coverage data from a WIG file.

        Args:
            wig_file: Path to the WIG file

        Returns:
            Pandas Series with position-indexed coverage data
        """
        try:
            # Skip header lines and read the data
            coverage_df = pd.read_csv(
                wig_file,
                skiprows=3,
                delimiter='\t',
                index_col=0,
                header=None,
                usecols=range(6)
            )

            # Set column names
            coverage_df.columns = ['A', 'C', 'G', 'T', 'N']

            # Sum across all nucleotides to get total coverage at each position
            total_coverage = coverage_df.sum(axis=1)

            return total_coverage

        except (FileNotFoundError, pd.errors.EmptyDataError) as e:
            self.logger.error(f"Error reading WIG file: {str(e)}")
            return None

    def _get_nextclade_info(self, nextclade_file: Optional[str]) -> Dict[str, str]:
        """
        Extract genotype information from NextClade output.

        Args:
            nextclade_file: Path to the NextClade TSV file or None

        Returns:
            Dictionary with NextClade genotype information
        """
        if not nextclade_file or not os.path.isfile(nextclade_file):
            return {'whole_genome_clade': 'Not RSV', 'g_clade': 'Not RSV'}

        try:
            nextclade_df = pd.read_csv(nextclade_file, index_col=0, header=0, sep=';')
            return {
                'whole_genome_clade': nextclade_df.iloc[0, 1],
                'g_clade': nextclade_df.iloc[0, 2]
            }
        except (FileNotFoundError, pd.errors.EmptyDataError, IndexError) as e:
            self.logger.error(f"Error reading NextClade file: {str(e)}")
            return {'whole_genome_clade': 'Not RSV', 'g_clade': 'Not RSV'}

    def _get_blast_info(self, genotype_file: Optional[str]) -> Dict[str, str]:
        """
        Extract genotype information from BLAST output.

        Args:
            genotype_file: Path to the BLAST output file or None

        Returns:
            Dictionary with BLAST genotype information
        """
        if not genotype_file or not os.path.isfile(genotype_file):
            return {'whole_genome_clade': 'Not RSV'}

        try:
            with open(genotype_file, 'r') as file:
                genotype_data = file.readline().strip().split(',')

            if len(genotype_data) >= 3 and int(genotype_data[2]) >= 1500:
                return {'whole_genome_clade': genotype_data[0]}
            else:
                return {'whole_genome_clade': 'Not RSV'}

        except (FileNotFoundError, IndexError, ValueError) as e:
            self.logger.error(f"Error reading BLAST genotype file: {str(e)}")
            return {'whole_genome_clade': 'Not RSV'}

    def _write_report_line(self, sample_name: str, qc_info: Dict[str, Any],
                           mapping_stats: Dict[str, float], subtype_info: Dict[str, str],
                           f_mutation_results: Dict[str, Any], gene_coverage: Dict[str, float],
                           nextclade_info: Dict[str, str], blast_info: Dict[str, str],
                           report_file: str) -> None:
        """
        Write a single sample's results to the report file.

        Args:
            sample_name: Sample name
            qc_info: Quality control information
            mapping_stats: Mapping statistics
            subtype_info: Subtype information
            f_mutation_results: F protein mutation results
            gene_coverage: Gene coverage data
            nextclade_info: NextClade genotype information
            blast_info: BLAST genotype information
            report_file: Path to the report file
        """
        # Format QC information
        qc_str = (f"{qc_info['before_total']},{qc_info['before_q20']},{qc_info['before_q30']},"
                  f"{qc_info['after_total']},{qc_info['after_q20']},{qc_info['after_q30']},{qc_info['qc_rate']:.2f}")

        # Format mapping statistics
        mapping_str = (f"{mapping_stats['uniquely_mapped']:.2f},{mapping_stats['multi_mapped']:.2f},"
                       f"{mapping_stats['unmapped']:.2f},{mapping_stats['chimeric']:.2f}")

        # Format F protein mutations
        mutations_str = f_mutation_results.get('mutation_str', '')

        # Check if F gene has low coverage
        if not mutations_str and 'CDS_8' in gene_coverage and gene_coverage['CDS_8'] < self.coverage_threshold:
            mutations_str = 'Low coverage of F gene'

        # Format gene coverage data
        gene_ids = ['CDS_1', 'CDS_2', 'CDS_3', 'CDS_4', 'CDS_5',
                    'CDS_6', 'CDS_7', 'CDS_8', 'CDS_9', 'CDS_10', 'CDS_11']
        coverage_values = [gene_coverage.get(gene, 0) for gene in gene_ids]
        coverage_str = ','.join(f"{value:.2f}" for value in coverage_values)

        # Write to report file
        with open(report_file, 'a') as file:
            file.write(
                f"{sample_name},{qc_str},{mapping_str},{subtype_info['subtype']},"
                f"{subtype_info['ref_accession']},{subtype_info['ref_subtype']},"
                f"{mutations_str},{coverage_str},{nextclade_info['whole_genome_clade']},"
                f"{blast_info['whole_genome_clade']}\n"
            )

    @staticmethod
    def _write_error_line(sample_name: str, qc_info: Dict[str, Any], report_file: str) -> None:
        """
        Write an error line for a sample that couldn't be processed.

        Args:
            sample_name: Sample name
            qc_info: Quality control information
            report_file: Path to the report file
        """
        qc_str = (f"{qc_info['before_total']},{qc_info['before_q20']},{qc_info['before_q30']},"
                  f"{qc_info['after_total']},{qc_info['after_q20']},{qc_info['after_q30']},{qc_info['qc_rate']:.2f}")

        with open(report_file, 'a') as file:
            file.write(
                f"{sample_name},{qc_str},0,0,100,0,Not RSV,NA,Not RSV,,"
                f"0,0,0,0,0,0,0,0,0,0,0,Not RSV,Not RSV\n"
            )

    def _generate_f_mutation_reports(self, f_mutations_a: Dict[str, List[str]],
                                     f_mutations_b: Dict[str, List[str]],
                                     output_files: Dict[str, str]) -> None:
        """
        Generate reports for F protein mutations.

        Args:
            f_mutations_a: Dictionary of F protein mutations for RSV-A samples
            f_mutations_b: Dictionary of F protein mutations for RSV-B samples
            output_files: Dictionary of output file paths
        """
        # Extract key residue positions
        self._extract_key_residues(
            output_files['f_protein_a'],
            os.path.join(self.reference_dir, 'F_mutation', 'RSV_A_F_Mutation.csv'),
            output_files['f_key_position_a']
        )

        self._extract_key_residues(
            output_files['f_protein_b'],
            os.path.join(self.reference_dir, 'F_mutation', 'RSV_B_F_Mutation.csv'),
            output_files['f_key_position_b']
        )

        # Generate F mutation report for RSV-A
        self._write_f_mutation_report(
            'A',
            f_mutations_a,
            output_files['f_report_a']
        )

        # Generate F mutation report for RSV-B
        self._write_f_mutation_report(
            'B',
            f_mutations_b,
            output_files['f_report_b']
        )

    def _extract_key_residues(self, fasta_file: str, mutation_file: str, output_csv: str) -> None:
        """
        Extract key residue positions from a protein sequence alignment.

        Args:
            fasta_file: Path to the protein sequence FASTA file
            mutation_file: Path to the mutation definition file
            output_csv: Path to the output CSV file
        """
        try:
            # Read the FASTA file
            sequence_names = []
            sequences = []

            for record in SeqIO.parse(fasta_file, "fasta"):
                sequence_names.append(record.id)
                sequences.append(list(str(record.seq)))

            if not sequence_names:
                return

            # Create a DataFrame from sequences
            seq_df = pd.DataFrame(sequences, index=sequence_names)
            seq_df.columns = range(1, seq_df.shape[1] + 1)  # 1-based positions

            # Read mutation file and extract positions
            with open(mutation_file, 'r') as file:
                content = file.read().replace('+', '\n')

            # Parse positions from mutation definitions
            mutation_data = pd.read_csv(
                StringIO(content),
                header=None,
                names=["Residue", "Position", "Substitutions"]
            )
            positions = mutation_data["Position"].astype(int).tolist()

            # Keep only unique positions, sorted
            unique_positions = sorted(set(positions))

            # Subset the DataFrame to keep only these positions
            result_df = seq_df[unique_positions]

            # Write to CSV
            result_df.to_csv(output_csv, index=True, header=True)

        except (FileNotFoundError, pd.errors.EmptyDataError) as e:
            self.logger.error(f"Error extracting key residues: {str(e)}")

    def _write_f_mutation_report(self, subtype: str, f_mutations: Dict[str, List[str]],
                                 output_file: str) -> None:
        """
        Write F protein mutation report to a CSV file.

        Args:
            subtype: RSV subtype (A or B)
            f_mutations: Dictionary mapping sample names to mutation lists
            output_file: Path to the output CSV file
        """
        if not f_mutations:
            return

        try:
            # Get mutation definitions
            mutation_file = os.path.join(
                self.reference_dir,
                'F_mutation',
                f"RSV_{subtype}_F_Mutation.csv"
            )

            # Read positions from mutation file
            positions = []
            with open(mutation_file, 'r') as file:
                for line in file:
                    parts = line.strip().split(',')
                    if len(parts) >= 2:
                        positions.append(parts[1])

            # Write header
            with open(output_file, 'w') as file:
                # Version information
                file.write(f"Pipeline version: {self.version}\n")

                # CSV header
                headers = ['Sample']
                for line in open(mutation_file, 'r'):
                    headers.append(re.sub(',', '', line.strip()))

                file.write(','.join(headers) + '\n')

                # Write sample data
                for sample, mutations in f_mutations.items():
                    # Initialize empty result list
                    result_list = [''] * len(positions)

                    # Fill in mutations
                    for mutation in mutations:
                        # Extract position number
                        match = re.search(r'\d+', mutation)
                        if match:
                            position = match.group()
                            if position in positions:
                                index = positions.index(position)
                                # Remove '(Reported)' from mutation string
                                mutation_str = mutation.replace('(Reported)', '')
                                result_list[index] = mutation_str

                    # Write sample line
                    file.write(f"{sample},{','.join(result_list)}\n")

        except (FileNotFoundError, IndexError) as e:
            self.logger.error(f"Error writing F mutation report: {str(e)}")

    def _prepare_phylogenetic_tree_files(self, subtype: str, sample_names: List[str],
                                         output_files: Dict[str, str]) -> None:
        """
        Prepare files for phylogenetic tree analysis.

        Args:
            subtype: RSV subtype (A or B)
            sample_names: List of sample names
            output_files: Dictionary of output file paths
        """
        if not sample_names:
            return

        try:
            # Define reference files
            reference_fasta = os.path.join(
                self.reference_dir,
                'TreeReference',
                f"representative_ref_{subtype}.fasta"
            )
            reference_anno = os.path.join(
                self.reference_dir,
                'TreeReference',
                f"representative_ref_{subtype}.csv"
            )

            # Define output files
            tree_seq_file = output_files[f'tree_seq_{subtype.lower()}']
            tree_anno_file = output_files[f'tree_anno_{subtype.lower()}']

            # Copy reference annotation file
            shutil.copyfile(reference_anno, tree_anno_file)

            # Copy sequences for this subtype
            sequence_file = output_files[f'sequence_{subtype.lower()}']
            shutil.copyfile(sequence_file, tree_seq_file)

            # Append reference sequences
            with open(tree_seq_file, 'a') as dst_file:
                with open(reference_fasta, 'r') as src_file:
                    shutil.copyfileobj(src_file, dst_file)

            # Add sample annotations
            with open(tree_anno_file, 'a') as anno_file:
                for sample in sample_names:
                    anno_file.write(f"{sample},RSVrecon,Query\n")

            # Add NEXT-RSV-PIPE sequences if available
            if self.next_rsv_pipe_res:
                self._add_next_rsv_sequences(sample_names, tree_seq_file, tree_anno_file)

        except FileNotFoundError as e:
            self.logger.error(f"Error preparing phylogenetic tree files: {str(e)}")

    def _add_next_rsv_sequences(self, sample_names: List[str],
                                tree_seq_file: str, tree_anno_file: str) -> None:
        """
        Add NEXT-RSV-PIPE sequences to the phylogenetic tree files.

        Args:
            sample_names: List of sample names
            tree_seq_file: Path to the tree sequence file
            tree_anno_file: Path to the tree annotation file
        """
        # Find FASTA files in NEXT-RSV-PIPE results
        fa_files = glob.glob(os.path.join(self.next_rsv_pipe_res, '*.fa'))
        if not fa_files:
            return

        with open(tree_seq_file, 'a') as seq_file, open(tree_anno_file, 'a') as anno_file:
            for filename in fa_files:
                with open(filename, 'r') as file:
                    lines = file.readlines()
                    current_sample = None
                    sequence_lines = []

                    for line in lines:
                        if line.startswith('>'):
                            # Process previous sequence if we have one
                            if current_sample and current_sample in sample_names:
                                seq_file.write(f">{current_sample}|NEXT_RSV\n")
                                anno_file.write(f"{current_sample}|NEXT_RSV,NEXT_RSV,Query\n")
                                for seq_line in sequence_lines:
                                    seq_file.write(seq_line)

                            # Extract new sample name
                            header = line.lstrip('>').rstrip()
                            current_sample = re.sub('_R_.+', '', header)
                            sequence_lines = []
                        else:
                            sequence_lines.append(line)

                    # Process the last sequence
                    if current_sample and current_sample in sample_names:
                        seq_file.write(f">{current_sample}|NEXT_RSV\n")
                        anno_file.write(f"{current_sample}|NEXT_RSV,NEXT_RSV,Query\n")
                        for seq_line in sequence_lines:
                            seq_file.write(seq_line)


def main():
    """Main entry point for the script."""
    parser = argparse.ArgumentParser(description='RSV Analysis Pipeline')
    parser.add_argument('--meta', required=True, help='Path to metadata file')
    parser.add_argument('--version', required=True, help='Path to version file')
    parser.add_argument('--output', required=True, help='Output directory')
    parser.add_argument('--reference', required=True, help='Reference directory')
    parser.add_argument('--coverage-threshold', type=int, default=50,
                        help='Minimum coverage threshold for gene analysis')
    parser.add_argument('--next-rsv', help='Path to NEXT-RSV-PIPELINE results')

    args = parser.parse_args()

    # Run the pipeline
    pipeline = RSVAnalysisPipeline(
        args.meta,
        args.version,
        args.output,
        args.reference,
        args.coverage_threshold,
        args.next_rsv
    )

    pipeline.run()


if __name__ == "__main__":
    main()
