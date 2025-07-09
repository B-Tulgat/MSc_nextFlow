#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Parameters
params.pod5 = "m6A_DRACH.pod5"
params.reference = "drach_context_strands.fa"
params.model_name = "rna004_130bps_sup@v5.2.0"
params.mod_model = "rna004_130bps_sup@v5.2.0_m6A_DRACH@v1"
params.model_path = "$PWD/models/${params.model_name}"
params.mod_model_path = "$PWD/models/${params.mod_model}"

// Process to setup directories and download models
process SETUP_AND_DOWNLOAD_MODELS {
    container 'nanopore:latest'

    input:
    val model_name
    val mod_model

    output:
    path model_name, emit: base_model
    path mod_model, emit: mod_model

    script:
    """
    mkdir -p models

    if [ ! -d "./${model_name}" ]; then
        echo "Downloading base model: ${model_name}"
        micromamba run -n nanopore dorado download --model ${model_name}
    else
        echo "Base model ${model_name} already exists"
    fi

    if [ ! -d "./${mod_model}" ]; then
        echo "Downloading modification model: ${mod_model}"
        micromamba run -n nanopore dorado download --model ${mod_model}
    else
        echo "Modification model ${mod_model} already exists"
    fi
    """
}

// Process for basecalling with modification detection
process BASECALL_WITH_MODS {
    container 'nanopore:latest'

    input:
    path base_model
    path mod_model
    path pod5_file
    path reference

    output:
    path "m6A_DRACH.bam", emit: raw_bam

    script:
    """
    micromamba run -n nanopore dorado basecaller ${base_model} ${pod5_file} \
        --modified-bases-models ${mod_model} \
        --reference ${reference} > m6A_DRACH.bam
    """
}

// Process to sort and index BAM file
process SORT_AND_INDEX_BAM {
    container 'nanopore:latest'

    input:
    path raw_bam

    output:
    path "m6A_DRACH.sorted.bam", emit: sorted_bam
    path "m6A_DRACH.sorted.bam.bai", emit: bam_index

    script:
    """
    micromamba run -n nanopore samtools sort ${raw_bam} -o m6A_DRACH.sorted.bam
    micromamba run -n nanopore samtools index m6A_DRACH.sorted.bam
    """
}

// Process to call modifications using modkit
process CALL_MODIFICATIONS {
    container 'nanopore:latest'

    input:
    path sorted_bam
    path bam_index
    path reference

    output:
    path "m6A.bed", emit: mod_bed

    script:
    """
    micromamba run -n nanopore modkit pileup ${sorted_bam} m6A.bed --ref ${reference}
    """
}

// Workflow
workflow {
    pod5_ch = Channel.fromPath(params.pod5)
    reference_ch = Channel.fromPath(params.reference)

    SETUP_AND_DOWNLOAD_MODELS(params.model_name, params.mod_model)

    BASECALL_WITH_MODS(
        SETUP_AND_DOWNLOAD_MODELS.out.base_model,
        SETUP_AND_DOWNLOAD_MODELS.out.mod_model,
        pod5_ch,
        reference_ch
    )

    SORT_AND_INDEX_BAM(BASECALL_WITH_MODS.out.raw_bam)

    CALL_MODIFICATIONS(
        SORT_AND_INDEX_BAM.out.sorted_bam,
        SORT_AND_INDEX_BAM.out.bam_index,
        reference_ch
    )

    CALL_MODIFICATIONS.out.mod_bed.view { "✅ Pipeline completed. BED file: $it" }
}
