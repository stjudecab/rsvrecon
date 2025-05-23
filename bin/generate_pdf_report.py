#!/usr/bin/env python3
"""
RSV Report Generator

This module generates comprehensive PDF reports for RSV (Respiratory Syncytial Virus)
detection and analysis from clinical samples. It includes summary information,
phylogenetic analysis, and detailed coverage statistics for each sample.

Author: Haidong Yi (hyi@stjude.org)
        Lei Li (lei.li@stjude.org)
Date: May 22, 2025
"""

import os
import logging
from pathlib import Path
from datetime import datetime
from typing import List, Tuple, Dict, Union
from dataclasses import dataclass

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib import gridspec
from matplotlib.colors import LinearSegmentedColormap, BoundaryNorm
from PIL import Image as PILImage

from reportlab.lib.colors import black
from reportlab.lib.styles import getSampleStyleSheet, TA_LEFT, ParagraphStyle
from reportlab.platypus import (
    SimpleDocTemplate, Paragraph, Spacer, Flowable, Table,
    Image, TableStyle, PageBreak, ActionFlowable
)
from reportlab.lib import colors
from reportlab.lib.units import inch


# Configure the logger
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    handlers=[
        logging.FileHandler("generate_report.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)


@dataclass
class GenomeInfo:
    """Store genome interval information for visualization and analysis."""
    intervals: List[Tuple[int, int]]
    xticks: List[int]
    gene_names: List[str]


# Constants
COVERAGE_COLORMAP = {
    'colors': ['#D3D3D3', '#d1fbd8', '#a3ee5b', '#40ba0f'],
    'bins': [0, 50, 500, 1000, 50000]
}

GENOME_INFO = {
    'SubtypeA': GenomeInfo(
        intervals = [(70, 420), (599, 375), (1111, 1176), (2318, 726), (3226, 771), (4266, 195), (4652, 966),
            (5697, 1725), (7640, 585), (8158, 273), (8532, 6498)],
        xticks = [70, 800, 1600, 2700, 3550, 4366, 5000, 6300, 7650, 8400, 11500],
        gene_names = ['NS1', 'NS2', 'N', 'P', 'M', 'SH', 'G', 'F', 'M2-1', 'M2-2', 'L']
    ),
    'SubtypeB': GenomeInfo(
        intervals = [(57, 420), (584, 375), (1097, 1176), (2305, 726), (3220, 771), (4259, 198), (4646, 933),
            (5676, 1725), (7627, 588), (8180, 273), (8518, 6501)],
        xticks = [70, 800, 1600, 2700, 3550, 4366, 5000, 6300, 7650, 8400, 11500],
        gene_names = ['NS1', 'NS2', 'N', 'P', 'M', 'SH', 'G', 'F', 'M2-1', 'M2-2', 'L']
    )
}
GENOME_INFO['default'] = GENOME_INFO['SubtypeA']

GENOME_COLORS = [
    'red', 'yellow', 'green', 'cyan', 'blue', 'magenta',
    '#CCCC00', '#800080', '#D2B48C', '#8B4513', '#ADD8E6'
]


def color_gradient(value: float) -> colors.Color:
    """Return a color from green to red based on the input value (0 to 1)."""
    return colors.Color(1 - value, value, 0)


def color_gradient_matplotlib(value: float) -> np.ndarray:
    """Create a matplotlib color gradient from red to green."""
    cmap = mcolors.LinearSegmentedColormap.from_list("", ["#ff4d4d", "#4dff4d"])
    return cmap(value / 100)


def percentage_to_number(percentage: str) -> float:
    """Convert percentage string to float value."""
    return float(percentage.strip('%')) / 100


def make_coverage_heatmap(df: pd.DataFrame, png_file: str) -> None:
    """
    Create a heatmap visualization of gene coverage.

    Args:
        df: DataFrame with coverage data
        png_file: Path to save the output PNG file
    """
    number_of_rows = df.shape[0]
    height = 2 + number_of_rows * 0.15
    df.columns = df.columns.str.replace('_cov', '')

    # Create custom colormap
    colors = COVERAGE_COLORMAP['colors']
    n_bins = COVERAGE_COLORMAP['bins']
    cmap = LinearSegmentedColormap.from_list('custom_cmap', colors, N=len(n_bins) - 1)
    norm = BoundaryNorm(boundaries=n_bins, ncolors=len(n_bins) - 1)

    # Create figure
    plt.figure(figsize=(9, height))
    plt.rcParams['font.size'] = 7
    ax = sns.heatmap(
        df, annot=False, cmap=cmap, norm=norm,
        cbar_kws={"aspect": 20, "shrink": 0.5},
        linewidths=0.2, linecolor='black'
    )

    # Style the figure
    plt.yticks(rotation=0)
    plt.xticks(rotation=0)

    # Highlight specific columns
    xticks = ax.get_xticklabels()
    xticks[6].set_fontweight('bold')
    xticks[7].set_fontweight('bold')
    ax.set_xticklabels(xticks)

    # Draw custom grid lines
    for y in range(len(df)):
        plt.plot([6, 6], [y - 1, y + 1], color='black', linewidth=1)
        plt.plot([8, 8], [y - 1, y + 1], color='black', linewidth=1)

    # Save and close
    plt.savefig(png_file, dpi=300, bbox_inches='tight')
    plt.close()


def get_best_blast_hit(blast_res: str) -> Tuple[str, float, float]:
    """
    Extract the best BLAST hit from results file.

    Args:
        blast_res: Path to BLAST results file

    Returns:
        Tuple containing (reference_name, alignment_length, percent_identity)
    """
    blast_df = pd.read_csv(blast_res, sep='\t', header=None)
    blast_df.columns = [
        'Query', 'Subject', 'Pct_identity', 'Alignment_length',
        'Number_of_mismatches', 'Number_of_gap', 'Start_in_query',
        'End_in_query', 'Start_in_subject', 'End_in_subject',
        'E_value', 'Bit_score'
    ]
    blast_df = blast_df.sort_values(by=['Bit_score'], ascending=False)

    selected_ref_name = blast_df.loc[0, 'Subject']
    return (
        selected_ref_name,
        blast_df.loc[0, 'Alignment_length'],
        blast_df.loc[0, 'Pct_identity']
    )


def float_to_percentage(value: float) -> str:
    """Format float as percentage string with 2 decimal places."""
    return f"{value:.2%}"


def count_greater_than_n(array: Union[np.ndarray, List], n: int) -> int:
    """Count values in array greater than n."""
    return len([x for x in array if x > n])


def SNP_calling(wig_file: str, cutoff: float, genotype_text: str,
                gff_path: str, out_path: str, prefix_name: str,
                cov_cutoff: int) -> List[List]:
    """
    Generate SNP calling data and visualization.

    Args:
        wig_file: Path to coverage WIG file
        cutoff: SNP calling cutoff value
        genotype_text: Subtype text (SubtypeA or SubtypeB)
        gff_path: Path to GFF file
        out_path: Output directory
        prefix_name: Prefix for output files
        cov_cutoff: Coverage cutoff

    Returns:
        List of lists with SNP data for tabular display
    """
    # Get genome info for the genotype
    genome_info = GENOME_INFO.get(genotype_text, GENOME_INFO['default'])

    # Read coverage data
    cov_df = pd.read_csv(wig_file, skiprows=3, delimiter='\t', index_col=0, header=None, usecols=range(6))
    column_names = ['A', 'C', 'G', 'T', 'N']
    cov_df.columns = column_names
    cov = cov_df.sum(axis=1)
    percentage_df = cov_df.divide(cov, axis=0)
    percentage_df['cov'] = cov

    # remove columns that are lowly covered
    cov_df = cov_df[percentage_df['cov'] > cov_cutoff]
    percentage_df = percentage_df[percentage_df['cov'] > cov_cutoff]
    del percentage_df['cov']

    # Filter rows with >= 2 columns above cutoff
    num_columns_to_check = 2
    filtered_pct_df = percentage_df[(percentage_df > cutoff).sum(axis=1) >= num_columns_to_check]
    filtered_df = cov_df[(percentage_df > cutoff).sum(axis=1) >= num_columns_to_check]

    #  Get cutoff from GFF and apply
    end_pos_cds = last_pos_gff(gff_path)
    filtered_pct_df = filtered_pct_df[filtered_pct_df.index < end_pos_cds]
    filtered_df = filtered_df[filtered_df.index < end_pos_cds]
    cov = filtered_df.sum(axis=1)
    filtered_df['cov'] = cov


    plt.figure(figsize=(12, 4))
    gs = gridspec.GridSpec(2, 1, height_ratios=[4, 0.5], hspace=0.05)
    font_axis = {'family': 'sans-serif', 'color': 'black', 'weight': 'bold', 'size': 12}

    # Plot 1: Coverage and base distribution
    gs1 = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec=gs[0])
    ax_cov = plt.subplot(gs1[0])
    ax_cov.set_ylabel('Coverage', fontdict=font_axis)

    ax_cov.bar(filtered_pct_df.index, filtered_pct_df['A'], label='A', width=100)
    ax_cov.bar(filtered_pct_df.index, filtered_pct_df['C'], bottom=filtered_pct_df['A'], label='C', width=100)
    ax_cov.bar(
        filtered_pct_df.index, filtered_pct_df['G'],
        bottom=filtered_pct_df['A'] + filtered_pct_df['C'], label='G', width=100
    )
    ax_cov.bar(
        filtered_pct_df.index, filtered_pct_df['T'],
        bottom=filtered_pct_df['A'] + filtered_pct_df['C'] + filtered_pct_df['G'], label='T', width=100
    )
    ax_cov.legend()
    ax_cov.set_xticks([])

    # Plot 2: Gene map
    gs2 = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs[1], hspace=0.05)
    ax_genome = plt.subplot(gs2[0], sharex=ax_cov)

    for interval, color in zip(genome_info.intervals, GENOME_COLORS):
        ax_genome.broken_barh([interval], (10, 9), facecolors=color, edgecolors=color)

    ax_genome.set_xticks(genome_info.xticks)
    ax_genome.set_xticklabels(genome_info.gene_names)
    ax_genome.set_yticks([])
    ax_genome.set_ylabel('Genes', rotation=0, fontdict=font_axis, labelpad=30)
    ax_genome.yaxis.set_label_coords(0, -0.1)

    # Remove borders
    for spine in ax_genome.spines.values():
        spine.set_visible(False)

    # Save figure
    png_file = os.path.join(out_path, f"{prefix_name}_snp_figure.png")
    plt.savefig(png_file, dpi=300)
    plt.close()

    # Prepare data table
    gene_positions = parse_gff(gff_path)
    data = [['Position', 'Coverage', 'A', 'A%', 'C', 'C%', 'G', 'G%', 'T', '%T', 'Gene']]

    for index_label in filtered_df.index:
        sub_list = [
            int(index_label), int(filtered_df.loc[index_label, 'cov']),
            int(filtered_df.loc[index_label, 'A']), float_to_percentage(filtered_pct_df.loc[index_label, 'A']),
            int(filtered_df.loc[index_label, 'C']), float_to_percentage(filtered_pct_df.loc[index_label, 'C']),
            int(filtered_df.loc[index_label, 'G']), float_to_percentage(filtered_pct_df.loc[index_label, 'G']),
            int(filtered_df.loc[index_label, 'T']), float_to_percentage(filtered_pct_df.loc[index_label, 'T']),
            find_gene_at_position(gene_positions, index_label)
        ]

        data.append(sub_list)

    return data

def last_pos_gff(gff_file: str) -> int:
    """
    Get the end position of the last CDS in GFF file.

    Args:
        gff_file: Path to GFF file

    Returns:
        End position of the last CDS or 15000 if file cannot be read
    """
    try:
        gff_df = pd.read_csv(gff_file, sep='\t', comment='#', header=None)

        # Assign column names
        gff_df.columns = [
            'seqid', 'source', 'type', 'start', 'end',
            'score', 'strand', 'phase', 'attributes'
        ]

        # Filter for CDS entries and get last position
        cds_df = gff_df[gff_df['type'] == 'CDS']
        last_cds_end = cds_df['end'].iloc[-1]
        return int(last_cds_end)
    except Exception as e:
        logger.warning(f"Error reading GFF file: {e}. Using default value 15000.")
        return 15000


def parse_gff(gff_file: str) -> Dict[str, Tuple[str, int, int]]:
    """
    Parse GFF file to extract gene positions.

    Args:
        gff_file: Path to GFF file

    Returns:
        Dictionary mapping gene names to (seq_id, start, end) tuples
    """
    gene_positions = {}
    with open(gff_file, 'r') as file:
        for line in file:
            if not line.startswith('#'):
                fields = line.strip().split('\t')
                if len(fields) >= 9 and fields[2] == 'gene':
                    seq_id = fields[0]
                    start = int(fields[3])
                    end = int(fields[4])
                    gene_name = fields[8].split(';')[0].split('=')[1]
                    gene_positions[gene_name] = (seq_id, start, end)
    return gene_positions

def find_gene_at_position(gene_positions: Dict[str, Tuple[str, int, int]],
                          position: int) -> str:
    """
    Find which gene contains a given position.

    Args:
        gene_positions: Dictionary of gene positions from parse_gff
        position: Position to query

    Returns:
        Gene name or "Not in CDS" if not found
    """
    for gene_name, (seq_id, start, end) in gene_positions.items():
        if start <= position <= end:
            return gene_name
    return 'Not in CDS'


@dataclass
class ReportConfig:
    """Configuration class for report generation parameters."""
    csv_file: str
    working_folder: str
    temp_folder: str
    asset_folder: str
    meta_file: str
    version_file: str
    logo_file: str
    tree_a_file: str
    tree_b_file: str
    igv_cutoff: int = 50
    line_width: int = 440
    snp_cutoff: float = 0.2


@dataclass
class GenomeIntervals:
    """Genome interval definitions for different RSV subtypes."""
    SUBTYPE_A = [
        (70, 420), (599, 375), (1111, 1176), (2318, 726), (3226, 771),
        (4266, 195), (4652, 966), (5697, 1725), (7640, 585), (8158, 273), (8532, 6498)
    ]
    SUBTYPE_B = [
        (57, 420), (584, 375), (1097, 1176), (2305, 726), (3220, 771),
        (4259, 198), (4646, 933), (5676, 1725), (7627, 588), (8180, 273), (8518, 6501)
    ]
    GENE_NAMES = ['NS1', 'NS2', 'N', 'P', 'M', 'SH', 'G', 'F', 'M2-1', 'M2-2', 'L']
    GENOME_XTICKS = [70, 800, 1600, 2700, 3550, 4366, 5000, 6300, 7650, 8400, 11500]
    GENOME_COLORS = [
        'red', 'yellow', 'green', 'cyan', 'blue', 'magenta',
        '#CCCC00', '#800080', '#D2B48C', '#8B4513', '#ADD8E6'
    ]


class LineFlowable(Flowable):
    """Custom flowable for drawing horizontal lines in PDF reports."""

    def __init__(self, width: int, height: int = 0):
        """
        Initialize line flowable.

        Args:
            width: Width of the line
            height: Height position of the line
        """
        Flowable.__init__(self)
        self.width = width
        self.height = height

    def draw(self):
        """Draw the line on the canvas."""
        self.canv.line(0, self.height, self.width, self.height)


class RSVReportGenerator:
    """Main class for generating RSV detection PDF reports."""

    def __init__(self, config: ReportConfig):
        """
        Initialize the report generator.

        Args:
            config: Configuration object containing all necessary parameters
        """
        self.config = config
        self.resource_path = Path(self.config.asset_folder)
        self.logo_path = Path(self.config.logo_file)

        # Validate input files
        self._validate_inputs()

        # Setup temporary folder
        self._setup_temp_folder()

        # Load data
        self.meta_df = self._load_meta_data()
        self.csv_df = self._load_csv_data()

        # Get version information
        self.version_info = self._get_version_info()

    def _validate_inputs(self) -> None:
        """Validate that all required input files exist."""
        required_files = [
            self.config.csv_file,
            self.config.meta_file,
            self.config.version_file,
            self.config.logo_file
        ]

        for file_path in required_files:
            if not Path(file_path).exists():
                raise FileNotFoundError(f"Required file not found: {file_path}")

    def _setup_temp_folder(self) -> None:
        """Create temporary folder if it doesn't exist."""
        temp_path = Path(self.config.temp_folder)
        temp_path.mkdir(parents=True, exist_ok=True)
        logger.info(f"Using temporary folder: {temp_path}")

    def _load_meta_data(self) -> pd.DataFrame:
        """Load metadata from TSV file."""
        try:
            return pd.read_csv(self.config.meta_file, sep='\t', index_col=0)
        except Exception as e:
            logger.error(f"Failed to load metadata: {e}")
            raise

    def _load_csv_data(self) -> pd.DataFrame:
        """Load CSV data with proper formatting."""
        try:
            df = pd.read_csv(self.config.csv_file, skiprows=1)
            # Format QC rate as decimal
            df['QC rate'] = df['QC rate'] / 100
            return df
        except Exception as e:
            logger.error(f"Failed to load CSV data: {e}")
            raise

    @staticmethod
    def _get_version(version_file: str) -> str:
        """Read version from the version file."""
        with open(version_file, "r") as file:
            return file.readline().strip()

    def _get_version_info(self) -> Dict[str, str]:
        """Get version information from version file."""
        try:
            return {
                'version': self._get_version(self.config.version_file),
                'date': datetime.now().strftime('%Y-%m-%d')
            }
        except Exception as e:
            logger.warning(f"Could not read version info: {e}")
            return {'version': 'Unknown', 'date': 'Unknown'}

    def _create_document_header(self, elements: List) -> None:
        """Create the document header with title and logo."""
        styles = getSampleStyleSheet()

        # Title
        title_style = styles['Title']
        title_style.fontSize = 20
        title_style.alignment = TA_LEFT
        title = Paragraph("Detection of RSV from clinical samples", title_style)
        elements.append(title)

        # Logo
        if self.logo_path.exists():
            try:
                pil_img = PILImage.open(self.logo_path)
                original_width, original_height = pil_img.size
                new_height = 60
                new_width = original_width * (new_height / original_height)
                img_logo = Image(str(self.logo_path), width=new_width, height=new_height)
                img_logo.hAlign = 'LEFT'
                elements.append(img_logo)
            except Exception as e:
                logger.warning(f"Could not load logo: {e}")

        # Divider line
        line = LineFlowable(self.config.line_width)
        elements.append(line)

    def _create_summary_section(self, elements: List) -> None:
        """Create the summary section with sample overview table."""
        styles = getSampleStyleSheet()

        # Section title
        subtitle = Paragraph("Summary", styles['Heading2'])
        elements.append(subtitle)

        # Information text
        info_text = (
            f"Pipeline Version: {self.version_info['version']}<br/>"
            "Subtypes of each reference are highlighted in different colors: "
            "<font color='red'>Subtype A</font> and <font color='blue'>Subtype B</font><br/>"
            f"Genotype calling is based on "
            f"<a href='https://nextstrain.org/rsv/a/genome'><b>Nextstrain</b> and "
            f"<b>Nextclade3</b>, Data updated {self.version_info['date']}</a><br/>"
            "A clade label with * indicates the clade of best blast hit strain "
            "due to negative results from NextClade3<br/>"
            "Click the sample name to jump to the detail section<br/>"
        )

        paragraph = Paragraph(info_text, styles['BodyText'])
        elements.append(paragraph)
        elements.append(Spacer(1, 6))

        # Create summary table
        self._create_summary_table(elements)

    def _create_summary_table(self, elements: List) -> None:
        """Create the main summary table."""
        # Prepare data for table
        subset_df = self.csv_df[[
            'Sample name', 'QC rate', 'Uniquely mapped reads(%)', 'Subtype',
            'reference_accession', 'ref_subtype', 'Whole Genome Clade(NextClade)',
            'Whole Genome Clade(Blast)', 'G_cov'
        ]].copy()

        data = []
        custom_style = ParagraphStyle(
            name='CustomStyle',
            parent=getSampleStyleSheet()['BodyText'],
            fontSize=9
        )

        # Process each row
        for _, row in subset_df.iterrows():
            processed_row = self._process_table_row(row, custom_style)
            data.append(processed_row)

        # Add header
        headers = ['Sample name', 'Pass QC', 'Mapping rate', 'Clade', 'Sign']
        data.insert(0, headers)

        # Create and style table
        col_widths = [1.5 * inch, 0.7 * inch, 4 * inch, 0.8 * inch, 0.4 * inch]
        table = Table(data, colWidths=col_widths)

        # Apply styling
        self._style_summary_table(table, data)
        elements.append(table)
        elements.append(Spacer(1, 12))

    def _process_table_row(self, row: pd.Series, style: ParagraphStyle) -> List:
        """Process a single table row with images and links."""
        sample_name = row['Sample name']
        mapping_rate = row['Uniquely mapped reads(%)']
        subtype = row['Subtype']
        ref_subtype = row['ref_subtype']

        # Create mapping rate visualization
        mapping_image = self._create_mapping_visualization(
            sample_name, mapping_rate, subtype, ref_subtype
        )

        # Create status icon
        status_icon = self._get_status_icon(subtype, mapping_rate)

        # Create clickable sample name link
        sample_link = Paragraph(
            f"<a href='#{sample_name}'>{sample_name}</a>",
            style
        )

        # Determine genotype text
        if row['Whole Genome Clade(NextClade)'] in ['A', 'B']:
            genotype_text = row['Whole Genome Clade(Blast)'] + '*'
        else:
            genotype_text = row['Whole Genome Clade(NextClade)']

        return [
            sample_link,
            f'{row["QC rate"]:.2%}',
            mapping_image,
            genotype_text,
            status_icon
        ]

    def _create_mapping_visualization(self, sample_name: str, mapping_rate: float,
                                      subtype: str, ref_subtype: str) -> Image:
        """Create a horizontal bar chart for mapping rate visualization."""
        color_dict = {'SubtypeA': 'red', 'SubtypeB': 'blue', 'Not RSV': 'gray'}

        # Create plot
        plt.figure(figsize=(7, 0.3))
        x = [f"{sample_name}"]
        y = [mapping_rate]
        my_color = color_gradient_matplotlib(mapping_rate)

        plt.barh(x, y, color=[my_color])
        plt.xlim(0, 100)

        # Style plot
        plt.xlabel('')
        plt.ylabel('')
        plt.title('')
        plt.xticks([])
        plt.yticks(color=color_dict.get(ref_subtype, 'gray'))

        # Remove borders
        for spine in plt.gca().spines.values():
            spine.set_visible(False)

        # Add text label
        plt.text(mapping_rate + 1, 0, f'{mapping_rate}%',
                 color='black', va='center')

        # Save and return image
        png_file = Path(self.config.temp_folder) / f'{sample_name}_mapping_figure.png'
        plt.savefig(png_file, dpi=300)
        plt.close()

        return Image(str(png_file), width=4 * inch, height=0.2 * inch)

    def _get_status_icon(self, subtype: str, mapping_rate: float) -> Image:
        """Get appropriate status icon based on subtype and mapping rate."""
        if subtype == 'Not RSV':
            icon_name = 'error.png'
        elif mapping_rate > 80:
            icon_name = 'correct.png'
        else:
            icon_name = 'warning.png'

        icon_path = self.resource_path / icon_name
        return Image(str(icon_path), width=0.2 * inch, height=0.2 * inch)

    @staticmethod
    def _style_summary_table(table: Table, data: List) -> None:
        """Apply styling to the summary table."""
        style = TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), colors.grey),
            ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
            ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
            ('FONTSIZE', (0, 0), (-1, -1), 9),
            ('BOTTOMPADDING', (0, 0), (-1, 0), 12),
            ('BACKGROUND', (0, 1), (-1, -1), colors.white),
            ('GRID', (0, 0), (-1, -1), 1, colors.black)
        ])

        # Add background colors based on values
        for row in range(1, len(data)):
            # QC rate column
            qc_value = percentage_to_number(data[row][1])
            style.add('BACKGROUND', (1, row), (1, row), color_gradient(qc_value))

            # Clade columns
            clade_value = data[row][3][0] if data[row][3] else 'N'
            color_dict = {'A': '#ffcccc', 'B': '#ccccff', 'N': '#d3d3d3'}
            clade_color = color_dict.get(clade_value, '#d3d3d3')
            style.add('BACKGROUND', (3, row), (4, row), clade_color)

        table.setStyle(style)

    def _create_coverage_section(self, elements: List) -> None:
        """Create the coverage heatmap section."""
        styles = getSampleStyleSheet()

        elements.append(PageBreak())
        subtitle = Paragraph("Coverage by genes", styles['Heading2'])
        elements.append(subtitle)

        # Prepare coverage data
        coverage_columns = [
            "NS1_cov", "NS2_cov", "N_cov", "P_cov", "M_cov",
            "SH_cov", "G_cov", "F_cov", "M2-1_cov", "M2-2_cov", "L_cov"
        ]

        plot_df = self.csv_df[coverage_columns].copy()
        plot_df.index = self.csv_df['Sample name']

        # Create and add heatmap
        png_file = Path(self.config.temp_folder) / 'coverage_heatmap.png'
        make_coverage_heatmap(plot_df, str(png_file))

        # Resize image appropriately
        img_heatmap = self._resize_image(png_file, max_width=480, max_height=650)
        img_heatmap.hAlign = 'LEFT'
        elements.append(img_heatmap)

    def _resize_image(self, image_path: Path, max_width: int, max_height: int) -> Image:
        """Resize image while maintaining aspect ratio."""
        pil_img = PILImage.open(image_path)
        original_width, original_height = pil_img.size

        # Calculate new dimensions
        width_ratio = max_width / original_width
        height_ratio = max_height / original_height
        scale_ratio = min(width_ratio, height_ratio)

        new_width = original_width * scale_ratio
        new_height = original_height * scale_ratio

        return Image(str(image_path), width=new_width, height=new_height)

    def _create_phylogenetic_section(self, elements: List) -> None:
        """Create phylogenetic analysis section if tree files exist."""
        tree_a_file = self.config.tree_a_file
        tree_b_file = self.config.tree_b_file
        if not tree_a_file and not tree_b_file:
            return

        styles = getSampleStyleSheet()

        elements.append(PageBreak())
        subtitle = Paragraph("Phylogenetic Analysis", styles['Heading2'])
        elements.append(subtitle)

        info_text = "Phylogenetic trees are generated using whole genome sequences."
        paragraph = Paragraph(info_text, styles['BodyText'])
        elements.append(paragraph)

        # Add subtype A tree if exists
        if tree_a_file and Path(tree_a_file).exists():
            subtitle_a = Paragraph("Subtype A", styles['Heading3'])
            elements.append(subtitle_a)

            img_tree_a = self._resize_image(
                Path(self.config.tree_a_file), max_width=480, max_height=600
            )
            img_tree_a.hAlign = 'LEFT'
            elements.append(img_tree_a)

        # Add subtype B tree if exists
        if tree_b_file and Path(tree_b_file).exists():
            if tree_a_file and Path(tree_a_file).exists():
                elements.append(PageBreak())

            subtitle_b = Paragraph("Subtype B", styles['Heading3'])
            elements.append(subtitle_b)

            img_tree_b = self._resize_image(
                Path(self.config.tree_b_file), max_width=480, max_height=600
            )
            img_tree_b.hAlign = 'LEFT'
            elements.append(img_tree_b)

    def _create_details_section(self, elements: List) -> None:
        """Create detailed analysis section for each sample."""
        styles = getSampleStyleSheet()

        elements.append(PageBreak())
        subtitle = Paragraph("Details", styles['Heading2'])
        elements.append(subtitle)

        csv_df_indexed = pd.read_csv(self.config.csv_file, skiprows=1, index_col=0)
        samples = self.meta_df.index.values

        for i, sample in enumerate(sorted(samples)):
            if i != 0:
                elements.append(PageBreak())
            self._create_sample_detail(elements, sample, csv_df_indexed, styles)

    def _create_sample_detail(self, elements: List, sample: str,
                              csv_df: pd.DataFrame, styles) -> None:
        """Create detailed analysis for a single sample."""
        # Sample header
        sample_title = Paragraph(
            f"<a name='{sample}'/>Sample: {sample}",
            styles['Heading3']
        )
        elements.append(sample_title)

        # Genotype information
        self._add_genotype_section(elements, sample, csv_df, styles)

        # QC details
        self._add_qc_section(elements, sample, csv_df, styles)

        # Coverage analysis (if not "Not RSV")
        ref_subtype = csv_df.loc[sample, 'ref_subtype']
        if ref_subtype != "Not RSV":
            self._add_coverage_section(elements, sample, csv_df, styles)
            self._add_snp_section(elements, sample, csv_df, styles)

    def _add_genotype_section(self, elements: List, sample: str,
                              csv_df: pd.DataFrame, styles) -> None:
        """Add genotype calling information."""
        subtitle = Paragraph('Genotype calls', styles['Heading4'])
        elements.append(subtitle)

        # Determine genotype text
        nextclade_result = csv_df.loc[sample, 'Whole Genome Clade(NextClade)']
        if nextclade_result in ['A', 'B']:
            genotype_text = csv_df.loc[sample, 'Whole Genome Clade(Blast)'] + '*'
        else:
            genotype_text = nextclade_result

        # Format F protein mutations
        f_mutations = csv_df.loc[sample, 'F protein mutations']
        if isinstance(f_mutations, str):
            f_mutations = f_mutations.replace("|", ", ").replace('(Reported)', '')
        else:
            f_mutations = ''

        # Create genotype paragraph with appropriate icon
        icon_size = '20'
        if genotype_text == "Not RSV":
            icon_path = self.resource_path / 'error.png'
            genotype_para = (f'<img src="{icon_path}" valign="middle" '
                             f'width="{icon_size}" height="{icon_size}"/>  {genotype_text}')
        else:
            mapping_rate = int(csv_df.loc[sample, 'Uniquely mapped reads(%)'])
            if mapping_rate > 80:
                icon_path = self.resource_path / 'correct.png'
            else:
                icon_path = self.resource_path / 'warning.png'

            genotype_para = (f'<img src="{icon_path}" valign="middle" '
                             f'width="{icon_size}" height="{icon_size}"/>  '
                             f'<b>{genotype_text}</b> (based on whole genome)<br/>'
                             f'F protein mutations: <b>{f_mutations}</b><br/><br/>')

        paragraph = Paragraph(genotype_para)
        elements.append(Spacer(1, 4))
        elements.append(paragraph)

        # Add BLAST hit information if not "Not RSV"
        if genotype_text != "Not RSV":
            self._add_blast_hit_table(elements, sample)

    def _add_blast_hit_table(self, elements: List, sample: str) -> None:
        """Add BLAST hit comparison table."""
        # Get BLAST results
        gisaid_blast = self.meta_df.loc[sample, 'blast_gisaid']
        nextstrain_blast = self.meta_df.loc[sample, 'whg_blastout']

        gisaid_name, gisaid_length, gisaid_identity = get_best_blast_hit(gisaid_blast)
        nextstrain_name, nextstrain_length, nextstrain_identity = get_best_blast_hit(nextstrain_blast)

        gisaid_name = gisaid_name.split('|')[0] if '|' in gisaid_name else gisaid_name

        data = [
            ['', 'Best hit', 'Identity(%)', 'Alignment Length'],
            ['GISAID', gisaid_name, gisaid_identity, gisaid_length],
            ['NextStrain', nextstrain_name, nextstrain_identity, nextstrain_length],
        ]

        table = Table(data)
        table.setStyle(self._get_standard_table_style())
        elements.append(table)

    def _add_qc_section(self, elements: List, sample: str,
                        csv_df: pd.DataFrame, styles) -> None:
        """Add QC details section."""
        subtitle = Paragraph('QC Details', styles['Heading4'])
        elements.append(subtitle)

        data = [
            ['', 'Total reads', 'Q20 pct', 'Q30 pct'],
            ['Raw data',
             csv_df.loc[sample, 'before_filtering_total_reads'],
             float_to_percentage(csv_df.loc[sample, 'before_filtering_q20_rate']),
             float_to_percentage(csv_df.loc[sample, 'before_filtering_q30_rate'])],
            ['Filtered data',
             csv_df.loc[sample, 'after_filtering_total_reads'],
             float_to_percentage(csv_df.loc[sample, 'after_filtering_q20_rate']),
             float_to_percentage(csv_df.loc[sample, 'after_filtering_q30_rate'])]
        ]

        table = Table(data, hAlign='LEFT')
        table.setStyle(self._get_standard_table_style())
        elements.append(table)

    def _add_coverage_section(self, elements: List, sample: str,
                              csv_df: pd.DataFrame, styles) -> None:
        """Add coverage analysis section."""
        subtitle = Paragraph('Coverage summary', styles['Heading4'])
        elements.append(subtitle)

        # Reference genome information
        ref_accession = csv_df.loc[sample, 'reference_accession']
        ref_text = (f"Reference genome: <b>"
                    f"<a href='https://www.ncbi.nlm.nih.gov/nuccore/{ref_accession}'>"
                    f"{ref_accession}, click for details</a></b>")

        paragraph = Paragraph(ref_text, styles['BodyText'])
        elements.append(paragraph)

        # Get genome intervals based on subtype
        subtype = csv_df.loc[sample, 'Subtype']
        if subtype == 'SubtypeB':
            intervals = GenomeIntervals.SUBTYPE_B
        else:
            intervals = GenomeIntervals.SUBTYPE_A

        # Read and process coverage data
        wig_file = self.meta_df.loc[sample, 'igv_out']
        coverage_data = self._process_coverage_data(wig_file)

        # Add coverage statistics
        self._add_coverage_statistics(elements, coverage_data, styles)

        # Create and add coverage plot
        coverage_image = self._create_coverage_plot(
            sample, coverage_data, intervals
        )

        info_text = ("The coverage plot is under <b>Log scale</b>, "
                     "the red, green and blue lines indicate coverage = 10, 50 and 500<br/>")
        paragraph = Paragraph(info_text, styles['BodyText'])
        elements.append(paragraph)
        elements.append(coverage_image)

    @staticmethod
    def _process_coverage_data(wig_file: str) -> Tuple[List[int], List[int]]:
        """Process coverage data from WIG file."""
        try:
            cov_df = pd.read_csv(wig_file, skiprows=3, delimiter='\t', index_col=0)
            coverage = cov_df.sum(axis=1)
            positions = cov_df.index.tolist()
            return positions, coverage.tolist()
        except Exception as e:
            logger.error(f"Failed to process coverage data: {e}")
            return [], []

    @staticmethod
    def _add_coverage_statistics(elements: List,
                                 coverage_data: Tuple[List[int], List[int]],
                                 styles) -> None:
        """Add coverage statistics to the report."""
        positions, coverage = coverage_data

        if not coverage:
            return

        max_pos = max(positions) if positions else 0

        # Add divider line
        line = LineFlowable(300)
        elements.append(line)

        # Coverage statistics
        stats = [
            ('Genome size (BP)', str(max_pos)),
            ('20x coverage (BP)', f"{count_greater_than_n(coverage, 19)} "
                                  f"({float_to_percentage(count_greater_than_n(coverage, 19) / max_pos) if max_pos > 0 else '0%'})"),
            ('50x coverage (BP)', f"{count_greater_than_n(coverage, 49)} "
                                  f"({float_to_percentage(count_greater_than_n(coverage, 49) / max_pos) if max_pos > 0 else '0%'})"),
            ('100x coverage (BP)', f"{count_greater_than_n(coverage, 99)} "
                                   f"({float_to_percentage(count_greater_than_n(coverage, 99) / max_pos) if max_pos > 0 else '0%'})"),
            ('500x coverage (BP)', f"{count_greater_than_n(coverage, 499)} "
                                   f"({float_to_percentage(count_greater_than_n(coverage, 499) / max_pos) if max_pos > 0 else '0%'})"),
            ('Average sequencing depth', str(int(sum(coverage) / len(coverage))) if coverage else '0')
        ]

        for label, value in stats:
            text = f'<b>{label}:</b>   {value}'
            paragraph = Paragraph(text, styles['BodyText'])
            elements.append(paragraph)

        # Add closing divider
        elements.append(Spacer(1, 2))
        elements.append(line)
        elements.append(Spacer(1, 4))

    def _create_coverage_plot(self, sample: str,
                              coverage_data: Tuple[List[int], List[int]],
                              intervals: List[Tuple[int, int]]) -> Image:
        """Create coverage visualization plot."""
        positions, coverage = coverage_data

        if not coverage:
            # Return empty plot if no data
            plt.figure(figsize=(12, 4))
            plt.text(0.5, 0.5, 'No coverage data available',
                     ha='center', va='center', transform=plt.gca().transAxes)
            png_file = Path(self.config.temp_folder) / f'{sample}_coverage_figure.png'
            plt.savefig(png_file, dpi=300)
            plt.close()
            return Image(str(png_file), width=600, height=200)

        # Create figure with subplots
        fig = plt.figure(figsize=(12, 4))
        gs = gridspec.GridSpec(2, 1, height_ratios=[4, 0.5], hspace=0.05)
        font_axis = {'family': 'sans-serif', 'color': 'black', 'weight': 'bold', 'size': 12}

        # Coverage plot
        gs1 = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec=gs[0])
        ax_cov = plt.subplot(gs1[0])
        ax_cov.set_ylabel('Coverage', fontdict=font_axis)
        ax_cov.fill_between(positions, coverage, color='gray')
        ax_cov.set_xticks([])
        ax_cov.set_yscale('log')

        # Add reference lines
        ax_cov.axhline(y=10, color='red', linestyle='--')
        ax_cov.axhline(y=50, color='green', linestyle='--')
        ax_cov.axhline(y=500, color='blue', linestyle='--')

        # Gene annotation plot
        gs2 = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs[1], hspace=0.05)
        ax_genome = plt.subplot(gs2[0], sharex=ax_cov)

        for interval, color in zip(intervals, GenomeIntervals.GENOME_COLORS):
            ax_genome.broken_barh([interval], (10, 9),
                                  facecolors=color, edgecolors=color)

        ax_genome.set_xticks(GenomeIntervals.GENOME_XTICKS)
        ax_genome.set_xticklabels(GenomeIntervals.GENE_NAMES)
        ax_genome.set_yticks([])
        ax_genome.set_ylabel('Genes', rotation=0, fontdict=font_axis, labelpad=30)
        ax_genome.yaxis.set_label_coords(0, -0.1)

        # Remove borders
        for spine in ax_genome.spines.values():
            spine.set_visible(False)

        # Save plot
        png_file = Path(self.config.temp_folder) / f'{sample}_coverage_figure.png'
        plt.savefig(png_file, dpi=300)
        plt.close()

        return Image(str(png_file), width=600, height=200)

    def _add_snp_section(self, elements: List, sample: str,
                         csv_df: pd.DataFrame, styles) -> None:
        """Add SNP analysis section."""
        ref_subtype = csv_df.loc[sample, 'ref_subtype']
        subtype = csv_df.loc[sample, 'Subtype']

        if "Not RSV" in (subtype, ref_subtype):
            return

        # Get required files
        wig_file = self.meta_df.loc[sample, 'igv_out']
        gff_file = self.meta_df.loc[sample, 'ref_gff']

        try:
            # Perform SNP calling
            table_data = SNP_calling(
                wig_file,
                self.config.snp_cutoff,
                subtype,
                gff_file,
                self.config.temp_folder,
                sample,
                self.config.igv_cutoff
            )

            # Add section header
            subtitle = Paragraph('SNP details', styles['Heading4'])
            elements.append(subtitle)

            # Add SNP figure
            png_file = Path(self.config.temp_folder) / f'{sample}_snp_figure.png'
            if png_file.exists():
                img_snp = Image(str(png_file), width=600, height=200)
                elements.append(img_snp)

            # Add SNP table
            if table_data:
                table = Table(table_data)
                table.setStyle(self._get_snp_table_style())
                elements.append(table)

        except Exception as e:
            logger.error(f"Failed to process SNP data for sample {sample}: {e}")
            error_text = "SNP analysis could not be completed for this sample."
            paragraph = Paragraph(error_text, styles['BodyText'])
            elements.append(paragraph)

    def _get_standard_table_style(self) -> TableStyle:
        """Get standard table styling."""
        return TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), colors.grey),
            ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
            ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
            ('LEFTMARGIN', (0, 0), (-1, -1), 0),
            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
            ('FONTSIZE', (0, 0), (-1, 0), 9),
            ('BOTTOMPADDING', (0, 0), (-1, 0), 12),
            ('BACKGROUND', (0, 1), (-1, -1), colors.beige),
            ('GRID', (0, 0), (-1, -1), 1, colors.black)
        ])

    def _get_snp_table_style(self) -> TableStyle:
        """Get SNP table specific styling."""
        return TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), colors.lightblue),
            ('TEXTCOLOR', (0, 0), (-1, 0), colors.white),
            ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
            ('LINEABOVE', (0, 0), (-1, 0), 2, colors.black),
            ('LINEBELOW', (0, 0), (-1, 0), 1, colors.black),
            ('LINEBELOW', (0, -1), (-1, -1), 2, colors.black),
        ])

    def generate_report(self) -> str:
        """
        Generate the complete PDF report.

        Returns:
            str: Path to the generated PDF file
        """
        try:
            # Setup output path
            pdf_path = Path(self.config.working_folder) / "Report.pdf"
            doc = SimpleDocTemplate(str(pdf_path))
            elements = []

            logger.info("Starting PDF report generation...")

            # Create document sections
            self._create_document_header(elements)
            self._create_summary_section(elements)
            self._create_coverage_section(elements)
            self._create_phylogenetic_section(elements)
            self._create_details_section(elements)

            # Build the document
            logger.info("Building PDF document...")
            doc.build(elements)

            return str(pdf_path)

        except Exception as e:
            logger.error(f"Failed to generate report: {e}")
            raise


def main():
    """Main entry point for the report generator."""
    import argparse

    parser = argparse.ArgumentParser(description='Generate comprehensive RSV analysis pdf reports.')
    parser.add_argument('--csv', required=True, help='Path to CSV file with analysis results')
    parser.add_argument('--output-dir', required=True, help='Working directory for outputs')
    parser.add_argument('--temp-dir', required=True, help='Directory for temporary files')
    parser.add_argument('--asset-dir', required=True, help='Directory for asset files')
    parser.add_argument('--meta', required=True, help='Path to manifest file with metadata')
    parser.add_argument('--version', required=True, help='Path to version file')
    parser.add_argument('--logo', required=True, help='Path to logo file')
    parser.add_argument('--tree-a', required=False, default='', help='Path to phylogenetic tree image for subtype A')
    parser.add_argument('--tree-b', required=False, default='', help='Path to phylogenetic tree image for subtype B')
    parser.add_argument('--igv-cutoff', type=int, default=50, help='Coverage cutoff for IGV visualization')

    args = parser.parse_args()

    logger.info("Initializing RSV report generator")
    config = ReportConfig(
        csv_file=args.csv,
        meta_file=args.meta,
        working_folder=args.output_dir,
        temp_folder=args.temp_dir,
        asset_folder=args.asset_dir,
        logo_file=args.logo,
        version_file=args.version,
        tree_a_file=args.tree_a,
        tree_b_file=args.tree_b,
        igv_cutoff=args.igv_cutoff
    )

    try:
        # Create report generator and generate report
        generator = RSVReportGenerator(config)
        output_path = generator.generate_report()
        logger.info(f"Report generated successfully: {output_path}")

    except Exception as e:
        logger.error(f"Report generation failed: {e}")
        raise


if __name__ == "__main__":
    main()
