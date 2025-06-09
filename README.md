# stjudecab/rsvrecon

<h3 align="center">
    <img src="https://raw.githubusercontent.com/stjudecab/rsvrecon/main/assets/rsvrecon_logo.png" alt="rsvrecon_logo" width="250"/>
</h3>
<br>

<p align="center">
<a href="https://github.com/stjudecab/rsvrecon/actions/workflows/ci.yml"><img alt="GitHub Actions Workflow Status" src="https://img.shields.io/github/actions/workflow/status/stjudecab/rsvrecon/ci.yml?branch=main&style=for-the-badge&logo=github&label=Test&labelColor=363a4f&color=f2cdcd"></a>
<a href="https://github.com/stjudecab/rsvrecon/actions/workflows/linting.yml"><img alt="GitHub Actions Workflow Status" src="https://img.shields.io/github/actions/workflow/status/stjudecab/rsvrecon/linting.yml?branch=main&style=for-the-badge&logo=github&label=Lint&labelColor=363a4f&color=a6d189"></a>
<a href="https://github.com/stjudecab/rsvrecon/stargazers"><img alt="GitHub Repo stars" src="https://img.shields.io/github/stars/stjudecab/rsvrecon?style=for-the-badge&logo=starship&labelColor=363a4f&color=b7bdf8"></a>
<a href="https://github.com/stjudecab/rsvrecon/releases/latest"><img alt="GitHub Release" src="https://img.shields.io/github/v/release/stjudecab/rsvrecon?style=for-the-badge&logo=github&labelColor=363a4f&color=89dceb"></a>
<a href="https://github.com/stjudecab/rsvrecon/issues"><img alt="GitHub Issues or Pull Requests" src="https://img.shields.io/github/issues/stjudecab/rsvrecon?style=for-the-badge&logo=data%3Aimage%2Fsvg%2Bxml%3Bbase64%2CPHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHZpZXdCb3g9IjAgMCAyNTYgMjU2Ij4KPHBhdGggZD0iTTIxNiwzMlYxOTJhOCw4LDAsMCwxLTgsOEg3MmExNiwxNiwwLDAsMC0xNiwxNkgxOTJhOCw4LDAsMCwxLDAsMTZINDhhOCw4LDAsMCwxLTgtOFY1NkEzMiwzMiwwLDAsMSw3MiwyNEgyMDhBOCw4LDAsMCwxLDIxNiwzMloiIHN0eWxlPSJmaWxsOiAjQ0FEM0Y1OyIvPgo8L3N2Zz4%3D&labelColor=363a4f&color=f5a97f"></a>
<br/>
<a href="https://www.nextflow.io/"><img alt="Nextflow" src="https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A524.04.2-23aa62.svg?style=for-the-badge"></a>
<a href="https://www.docker.com/"><img alt="run with docker" src="https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker&style=for-the-badge"></a>
<a href="https://sylabs.io/docs/"><img alt="run with singularity" src="https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000&style=for-the-badge"></a>
<a href="https://github.com/stjudecab/rsvrecon/blob/main/LICENSE"><img alt="GitHub License" src="https://img.shields.io/github/license/stjudecab/rsvrecon?style=for-the-badge&labelColor=363a4f&color=eba0ac"></a>
<a href="https://cloud.seqera.io/launch?pipeline=https://github.com/stjudecab/rsvrecon"><img alt="Launch on Seqera Platform" src="https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7?style=for-the-badge"></a>
</p>

[//]: # "[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)"

## News and Updates

## Introduction

**stjudecab/rsvrecon** is a bioinformatics workflow developed to assemble and analyze genomic sequences of Respiratory
Syncytial Virus (RSV) from Next-Generation Sequencing (NGS) data. It identifies genomic variations within RSV samples
and highlights clinically relevant genomic features. To simplify interpretation, the workflow generates easy-to-understand
HTML and PDF reports summarizing the results.

Built using [Nextflow](https://www.nextflow.io), the pipeline offers scalability, portability, and reproducibility
across diverse computational infrastructures. Dependency management is simplified by employing containerization
technologies such as `Docker`, `Singularity`, and `Conda`.

This pipeline utilizes the [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) framework, featuring
modularized processes with independent software environments, thereby making updates and maintenance straightforward.
Processes are also integrated, whenever feasible, with the [nf-core/modules](https://github.com/nf-core/modules)
repository to enhance usability and foster community contributions.

The schematic overview of the **stjudecab/rsvrecon** workflow is shown below:

<p align="center">
    <img src="https://raw.githubusercontent.com/stjudecab/rsvrecon/main/assets/rsvrecon_workflow.png" alt="rsvrecon_workflow" />
</p>

## Pipeline Overview

Briefly, the `rsvrecon` pipeline performs the following major steps:

1. **Merge Raw Reads (optional)**
   Concatenate re-sequenced FastQ files ([cat](http://www.linfo.org/cat.html)).

2. **Quality Control of Raw Reads**
   Evaluate sequencing quality ([FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)).

3. **Adapter Trimming**
   Remove adapters and low-quality bases ([fastp](https://github.com/OpenGene/fastp)).

4. **RSV Database Mapping**
   Map reads against RSV-specific databases ([KMA](https://github.com/genomicepidemiology/kma)).

5. **Sequence Alignment**

   - Align reads to a reference ([BWA](https://bio-bwa.sourceforge.net/bwa.shtml)(default) or [STAR](https://github.com/alexdobin/STAR)).
   - Sort and index aligned sequences ([SAMtools](https://sourceforge.net/projects/samtools/files/samtools/)).
   - Alignment quality assessment ([SAMtools](https://sourceforge.net/projects/samtools/files/samtools/)).

6. **Genome Assembly**

   - Calculate genome coverage ([IGVtools](https://github.com/igvteam/igv)).
   - Perform reference-guided genome assembly ([custom script](https://github.com/stjudecab/rsvrecon/blob/dev/bin/assemble_sequence.py)).

7. **Variant Identification**
   Determine viral clade assignments, mutations, and sequence quality ([NextClade](https://github.com/nextstrain/nextclade)).

8. **Genotyping (Whole Genome & G-gene)**
   - Perform BLAST search against reference genomes ([BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastSearch)).
   - Execute multiple sequence alignment ([Mafft](https://mafft.cbrc.jp/alignment/server/index.html)).
   - Generate a phylogenetic tree ([FastTree](https://www.microbesonline.org/fasttree/)).

## Pipeline Usage

> [!NOTE]
> If you’re new to Nextflow or nf-core, please refer to the [nf-core installation guide](https://nf-co.re/docs/usage/installation).
> Run the workflow first with `-profile test` to ensure proper functionality before applying it to your actual data.

Prepare your sample metadata in a CSV format (`samplesheet.csv`) with FASTQ files:

```csv
sample,fastq_1,fastq_2
sample_1,sample_1_R1_001.fastq.gz,sample_1_R2_001.fastq.gz
sample_2,sample_2_R1_001.fastq.gz,sample_2_R2_001.fastq.gz
```

> [!NOTE]
> If a sample has multiple sequencing lanes or replicates, list each replicate in a separate row with the
> same `sample` ID. The pipeline will automatically merge these reads before analysis. Spaces in sample IDs will be
> automatically converted to underscores (`_`).

Run the workflow using the following command structure:

```bash
# For St. Jude HPC users, specify the institutional profile, e.g., '-profile stjude'
nextflow run stjudecab/rsvrecon \
    -profile <docker/singularity/.../institution_config> \
    --input samplesheet.csv \
    --outdir <OUTDIR> \
    <args>
```

> [!WARNING]
> Pipeline parameters should only be provided via CLI arguments or Nextflow's `-params-file` option. Custom
> configuration files (`-c`) can specify system and pipeline configurations **except for parameters**.
> For detailed guidance, consult the [configuration documentation](https://nf-co.re/docs/usage/getting_started/configuration#custom-configuration-files).

For further information, please check the complete [usage documentation](./docs/usage.md) and [pipeline output descriptions](./docs/output.md).

## Pipeline Output

Pipeline results are organized per sample ID within the specified output directory (`<OUTDIR>`), structured by analysis stages:

- `<OUTDIR>/<sample_id>/fastq`: Contains trimmed FASTQ files ready for alignment.
- `<OUTDIR>/<sample_id>/bam`: Includes aligned BAM files and index files (BAI).
- `<OUTDIR>/<sample_id>/qc`: Quality control outputs for raw and aligned reads.
- `<OUTDIR>/<sample_id>/reference`: Reference genomes and database files used for the analysis.
- `<OUTDIR>/<sample_id>/variant_calling`: Variant calling and identification results.
- `<OUTDIR>/<sample_id>/phylogeny_tree`: Results related to phylogenetic analyses.
- `<OUTDIR>/<sample_id>/assembly`: Genome assemblies and coverage summaries.
- _(and more)_

In addition to sample-level results, our pipeline delivers batch-level results as well, structured as follows:

- `<OUTDIR>/batch_qcs/`: Contains [MultiQC](https://seqera.io/multiqc/)-based QC reports for both raw and filtered (trimmed) data.
- `<OUTDIR>/batch_reports/`: Include combined results across samples as clinical reports in both `PDF` and `HTML` formats.

## Credits

**stjudecab/rsvrecon** was developed by Haidong Yi ([@HaidYi](https://github.com/HaidYi)) and Lei Li ([@LeiLi-Uchicago](https://github.com/LeiLi-Uchicago)) at the
[Center for Applied Bioinformatics (CAB)](https://www.stjude.org/research/why-st-jude/shared-resources/center-for-applied-bioinformatics-cab.html),
[St. Jude Children's Research Hospital](https://www.stjude.org/). The pipeline design incorporates community-driven best
practices, especially inspired by [nf-core](https://nf-co.re/). We also thank the wider CAB team for their valuable inputs and feedbacks.

![StJude_CAB](assets/report_logo.png)

## Contributions and Support

We encourage community contributions. Please adhere to [nf-core guidelines](https://nf-co.re/developers/guidelines) for
maintaining consistency. Suggestions and improvements are welcome through pull requests or issues on our [GitHub repository](https://github.com/stjudecab/rsvrecon).

## Disclaimer

The pipeline [logo](./assets/rsvrecon_logo.png) is initially generated through [ChatGPT](https://chatgpt.com/)'s
new `4o Image Generation` function using the pipeline introduction as the prompt.

## Citations

@article {Li2025.rsvrecon,
author = {Li, Lei and Yi, Haidong and Brazelton, Jessica N. and Webby, Richard and Hayden, Randall T. and Wu, Gang and Hijano, Diego R.},
title = {Bridging Genomics and Clinical Medicine: RSVrecon Enhances RSV Surveillance with Automated Genotyping and Clinically-important Mutation Reporting},
elocation-id = {2025.06.03.657184},
year = {2025},
doi = {10.1101/2025.06.03.657184},
publisher = {Cold Spring Harbor Laboratory},
URL = {https://www.biorxiv.org/content/early/2025/06/09/2025.06.03.657184},
journal = {bioRxiv}
}
