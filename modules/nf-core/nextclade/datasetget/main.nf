process NEXTCLADE_DATASETGET {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/93/936786744b34cf016b948026a6b4e9489011424e15c28dfb2f7d03c31bb4afb5/data' :
        'community.wave.seqera.io/library/nextclade:3.11.0--155203da8341cfe6' }"

    input:
    tuple val(meta), path(ref_id)
    val tag

    output:
    tuple val(meta), path("*.db") , emit: dataset
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def version = tag ? "--tag ${tag}" : ''
    """
    # extract the dataset from the read_id file using 'cut'
    strain_info = \$(cat ref_id)
    strain_id=\$(echo \$strain_info | cut -d ',' -f1)
    subtype=\$(echo \$strain_info | cut -d ',' -f2)
    ref_strain_id2=\$(echo \$strain_info | cut -d ',' -f3)
    ref_subtype=\$(echo \$strain_info | cut -d ',' -f4)

    # Determine the reference path based on the subtype
    if [[ "\$ref_subtype" == "SubtypeA" ]]; then
        dataset="rsv_a"
    elif [[ "\$ref_subtype" == "SubtypeB" ]]; then
        dataset="rsv_b"
    else
        echo "Unknown subtype: \$dataset"
        exit 1
    fi

    # download the nextclade rsv database
    nextclade \\
        dataset \\
        get \\
        $args \\
        --name \$dataset \\
        $version \\
        --output-dir \$dataset.nextclade.db

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nextclade: \$(echo \$(nextclade --version 2>&1) | sed 's/^.*nextclade //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${dataset}"
    """
    mkdir -p ${prefix}
    touch ${prefix}/CHANGELOG.md
    touch ${prefix}/README.md
    touch ${prefix}/genome_annotation.gff3
    touch ${prefix}/pathogen.json
    touch ${prefix}/reference.fasta
    touch ${prefix}/sequences.fasta
    touch ${prefix}/tree.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nextclade: \$(echo \$(nextclade --version 2>&1) | sed 's/^.*nextclade //; s/ .*\$//')
    END_VERSIONS
    """

}
