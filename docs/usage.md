# stjudecab/rsvrecon: Usage

## Preparing Your Sample Sheet

Before running the pipeline, you need to create a sample information sheet in `CSV` format (for example, `samplesheet.csv`).
Specify the location of this file using the `--input` parameter:

```
--input </path/to/samplesheet.csv>
```

### Sample Sheet Structure

Your sample sheet must contain exactly three columns (`sample`, `fastq_1`, and `fastq_2`) and a header row as shown below:

```csv
sample,fastq_1,fastq_2
sample_1,sample_1_R1_001.fastq.gz,sample_1_R2_001.fastq.gz
sample_2,sample_2_R1_001.fastq.gz,sample_2_R2_001.fastq.gz
```

> **Note:** Spaces in sample names will automatically be converted to underscores (`_`) by the pipeline to prevent
> potential downstream issues.

| Column    | Description                                                                                                                                    |
| --------- | ---------------------------------------------------------------------------------------------------------------------------------------------- |
| `sample`  | Unique identifier for each sample. If a sample has multiple sequencing libraries or runs, this ID must remain consistent across multiple rows. |
| `fastq_1` | Full path to the gzipped FastQ file containing Illumina read 1. The filename must end with `.fastq.gz` or `.fq.gz`.                            |
| `fastq_2` | Full path to the gzipped FastQ file containing Illumina read 2. The filename must end with `.fastq.gz` or `.fq.gz`.                            |

> [!NOTE]
> Multiple Sequencing Runs
> If the same sample has been sequenced multiple times (e.g., across different lanes or runs), include each sequencing
> run as a separate row in the samplesheet using the same sample ID. The pipeline will automatically concatenate the
> reads from these runs before proceeding with downstream analysis.

## Run the pipeline

The command used for running the pipeline is as follows:

```bash
# For St. JUDE HPC users, specify the stjude profile: -profile stjude
nextflow run stjudecab/rsvrecon \
    -r <VERSION> \
    -profile <docker/singularity/.../intitution_config> \
    --input samplesheet.csv \
    --outdir <OUTDIR> \
    <args>
```

- `-r <VERSION>`: (Optional but recommended) Specifies the pipeline version for reproducibility.
- `-profile <PROFILE>`: Required. Select the configuration profile (e.g., your institution profile). If you use St. Jude
  HPC, you can use St. Jude's profile with `-profile stjude`.
- `--input <samplesheet.csv>`: Path to your sample sheet. (see [samplesheet](#preparing-your-sample-sheet) for more details).
- `--outdir <output_dir>`: Directory for output files. The pipeline will create this directory if it doesn’t exist.

> [!IMPORTANT]
> If running the pipeline on St. Jude HPC, use the institution-level profile by specifying -profile stjude.
> For more details, see the [profile documentation](https://nf-co.re/configs/stjude/).

When you run the command, the pipeline will create the following items in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

### Use a Parameter File

To avoid retyping command-line options every time, you can define them in a `YAML` or `JSON` file.
Run the pipeline using your parameter file like this:

```bash
nextflow run stjudecab/rsvrecon -params-file params.yaml
```

For example, your `params.yaml` file might look like:

```yaml
input: "./samplesheet.csv"
outdir: "./results"
# Additional advanced parameters...
```

> [!WARNING]
> Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must
> only be used for [tuning process resource specifications](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources),
> other infrastructural tweaks (such as output directories), or module arguments (args).

You can also generate these parameter files via [nf-core/launch](https://nf-co.re/launch).

### Update the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version.
When running the pipeline after this, it will always use the cached version if available -
even if the pipeline has been updated since. To make sure that you’re running the latest version of the pipeline,
make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull stjudecab/rsvrecon
```

You can also run the pipeline from local Git repositories or other Git hosting services (GitLab, Bitbucket) by providing
the appropriate URL:

```bash
nextflow run http://your-git-server.com/username/repository
```

Or with HTTPS:

```bash
nextflow run https://your-git-server.com/username/repository
```

To target a specific branch, tag, or commit, add the `-r` option:

```bash
nextflow run https://your-git-server.com/username/repository -r branch_name
nextflow run https://your-git-server.com/username/repository -r v1.0
nextflow run https://your-git-server.com/username/repository -r a1b2c3d
```

Nextflow supports various Git hosting platforms including GitHub, GitLab, Bitbucket, and any other Git server,
as long as Nextflow can access it via HTTP/HTTPS or SSH protocols.

### Reproducibility

It is a good idea to specify the pipeline version when running the pipeline on your data. This ensures that a specific
version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll
be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [stjudecab/rsvrecon releases page](https://github.com/stjudecab/rsvrecon/releases) and find the latest
pipeline version - numeric only (e.g., `0.1.0`). Then specify this when running the pipeline with `-r` (one hyphen) - e.g.,
`-r 0.1.0`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you
look back in the future.

To further assist in reproducibility, you can use share and reuse [parameter files](#run-the-pipeline) to repeat
pipeline runs with the same settings without having to write out a command with every single parameter.

> [!TIP]
> If you wish to share such profile (such as upload as supplementary material for academic publications),
> make sure to NOT include cluster specific paths to files, nor institutional specific profiles.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs)
when it runs, making multiple config profiles for various institutional clusters available at run time.
For more information and to check if your system is supported, please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the
`PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer environment.

- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://hpc.github.io/charliecloud/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `wave`
  - A generic configuration profile to enable [Wave](https://seqera.io/wave/) containers. Use together with one of the above (requires Nextflow ` 24.03.0-edge` or later).
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same,
continuing from where it got to previously. For input to be considered the same, not only the names must be identical
but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data,
you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has
a default set of requirements for number of CPUs, memory and time. For most of the pipeline steps, if the job exits with
any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18)
it will automatically be resubmitted with higher resources request (2 x original, then 3 x original).
If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/usage/configuration#max-resources)
and [tuning workflow resources](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources) section of the nf-core website.

### Custom Containers

In some cases, you may wish to change the container or conda environment used by a pipeline steps for a particular tool.
By default, nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or
[bioconda](https://bioconda.github.io/) projects. However, in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the
[updating tool versions](https://nf-co.re/docs/usage/configuration#updating-tool-versions) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline.
Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the
[customising tool arguments](https://nf-co.re/docs/usage/configuration#customising-tool-arguments) section of the nf-core website.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation
are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea
to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please
can you test that the config file works with your pipeline of choice using the `-c` parameter.
You can then create a pull request to the `nf-core/configs` repository with the addition of your config file,
associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)),
and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

> [!TIP]
> Fortunately, if you are in St.Jude and intend to use the HPC to run this pipeline, St. Jude profile is available at
> [nf-core/configs](https://nf-co.re/configs/). You can use it via `-profile stjude` when launching jobs. For more
> details and functionality about this profile, please refer to the official [profile documentation](https://nf-co.re/configs/stjude/).

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about
creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on
the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop
if you log out of your session. The logs are saved to a file.

Alternatively, you can use `tmux` / `zellij` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
