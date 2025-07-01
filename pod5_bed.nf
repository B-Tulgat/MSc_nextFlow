#! /usr/bin/env nextflow
params.pod5 = "m6A_DRACH.pod5"
params.reference = "drach_context_strands.fa"
params.model_name = "rna004_130bps_sup@v5.2.0"
params.mod_model = "rna004_130bps_sup@v5.2.0_m6A_DRACH@v1"

workflow {
    setup_directories()
        download_models()
        basecall()
        sort_bam()
        index_bam()
        call_mods()
}

process setup_directories {
output:
    file "models"
        file "outdir"

        script:
        """
        mkdir -p models outdir
        """
}

process download_models {
input:
    file "models"

        output:
        file "models/${params.model_name}"
        file "models/${params.mod_model}"

        script:
        """
        dorado download ${params.model_name} --directory models
        dorado download ${params.mod_model} --directory models
        """
}

process basecall {
input:
    file pod5 from file(params.pod5)
        file reference from file(params.reference)
        file model from file("models/${params.model_name}")
        file mod_model from file("models/${params.mod_model}")

        output:
        file "outdir/m6A_DRACH.bam"

        script:
        """
        dorado basecaller models/${params.model_name} $pod5 \
        --modified-bases-models models/${params.mod_model} \
        --reference $reference > outdir/m6A_DRACH.bam
        """
}

process sort_bam {
input:
    file bam from file("outdir/m6A_DRACH.bam")

        output:
        file "outdir/m6A_DRACH.sorted.bam"

        script:
        """
        samtools sort -o outdir/m6A_DRACH.sorted.bam $bam
        """
}

process index_bam {
input:
    file sorted_bam from file("outdir/m6A_DRACH.sorted.bam")

        output:
        file "outdir/m6A_DRACH.sorted.bam.bai"

        script:
        """
        samtools index $sorted_bam
        """
}

process call_mods {
input:
    file sorted_bam from file("outdir/m6A_DRACH.sorted.bam")
        file reference from file(params.reference)

        output:
        file "outdir/m6A.bed"

        script:
        """
        modkit pileup $sorted_bam outdir/m6A.bed --ref $reference
        """
}
