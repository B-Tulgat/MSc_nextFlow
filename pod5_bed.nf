#! /usr/bin/env nextflow

nextflow.enable.dsl = 2

// Pipeline parameters
params.pod5_dir = null
params.reference_fa = null
params.bam_file = null
params.outdir = "results"
params.min_depth = 25
params.dorado_model = "rna004_130bps_sup@v3.0.1_m6A_DRACH@v1"
params.minimap2_preset = "map-ont"
params.minimap2_kmer = 14

// Validate required parameters
if (!params.pod5_dir) error "Please specify --pod5_dir"
if (!params.reference_fa) error "Please specify --reference_fa"

log.info """\
        NANOPORE RNA-SEQ ANALYSIS PIPELINE
        ===================================
        pod5_dir        : ${params.pod5_dir}
reference_fa    : ${params.reference_fa}
bam_file        : ${params.bam_file}
outdir          : ${params.outdir}
dorado_model    : ${params.dorado_model}
min_depth       : ${params.min_depth}
"""

// Process 1: Basecalling with Dorado
process DORADO_BASECALLER {
    tag "basecalling"
        publishDir "${params.outdir}/basecalling", mode: 'copy'

        container 'nanopore-rnaseq:latest'

        input:
        path pod5_files

        output:
        path "basecalled.bam", emit: bam
        path "basecalling_summary.txt", emit: summary

        script:
        """
        dorado basecaller \\
        ${params.dorado_model} \\
        ${pod5_files} \\
        --output-dir . \\
        --summary-file basecalling_summary.txt > basecalled.bam
        """
}

// Process 2: Convert BAM to FASTQ with methylation tags
process BAM_TO_FASTQ {
    tag "bam2fastq"
        publishDir "${params.outdir}/fastq", mode: 'copy'

        container 'nanopore-rnaseq:latest'

        input:
        path bam

        output:
        path "converted.fastq", emit: fastq

        script:
        """
        samtools fastq \\
        -F 3840 \\
        -T MM,ML \\
        ${bam} > converted.fastq
        """
}

// Process 3: Alignment with minimap2
process MINIMAP2_ALIGN {
    tag "alignment"
        publishDir "${params.outdir}/alignment", mode: 'copy'

        container 'nanopore-rnaseq:latest'

        input:
        path fastq
        path reference

        output:
        path "aligned.bam", emit: bam
        path "aligned.bam.bai", emit: bai

        script:
        """
        minimap2 \\
        -ax ${params.minimap2_preset} \\
        -k ${params.minimap2_kmer} \\
        ${reference} \\
        ${fastq} | \\
        samtools sort -o aligned.bam -

        samtools index aligned.bam
        """
}

// Process 4: Generate bedMethyl files with modkit
process MODKIT_PILEUP {
    tag "modkit_pileup"
        publishDir "${params.outdir}/methylation", mode: 'copy'

        container 'nanopore-rnaseq:latest'

        input:
        path bam
        path bai
        path reference

        output:
        path "methylation.bed", emit: bedmethyl
        path "modkit_summary.txt", emit: summary

        script:
        """
        modkit pileup \\
        ${bam} \\
        methylation.bed \\
        --ref ${reference} \\
        --log-filepath modkit_summary.txt \\
        --threads ${task.cpus}
    """
}

// Process 5: Filter m6A sites and create BED files
process FILTER_M6A_SITES {
    tag "filter_m6a"
        publishDir "${params.outdir}/m6a_sites", mode: 'copy'

        container 'nanopore-rnaseq:latest'

        input:
        path bedmethyl
        val min_depth

        output:
        path "m6a_sites_filtered.bed", emit: filtered_bed
        path "m6a_statistics.txt", emit: stats

        script:
        """
#!/usr/bin/env python3
        import pandas as pd

# Read bedMethyl file
        df = pd.read_csv('${bedmethyl}', sep='\\t', header=None,
                names=['chrom', 'start', 'end', 'name', 'score', 'strand',
                'start2', 'end2', 'color', 'coverage', 'freq'])

# Filter for m6A sites with sufficient depth and non-zero frequency
        filtered = df[(df['coverage'] >= ${min_depth}) & (df['freq'] > 0)]

# Save filtered BED file
        filtered[['chrom', 'start', 'end', 'name', 'score', 'strand']].to_csv(
                'm6a_sites_filtered.bed', sep='\\t', header=False, index=False)

# Generate statistics
        with open('m6a_statistics.txt', 'w') as f:
        f.write(f"Total sites: {len(df)}\\n")
        f.write(f"Filtered sites (depth>={min_depth}, freq>0): {len(filtered)}\\n")
        f.write(f"Average coverage: {df['coverage'].mean():.2f}\\n")
        f.write(f"Average frequency: {df['freq'].mean():.4f}\\n")
        """
}

// Process 6: Create IGV visualization files
process CREATE_IGV_SESSION {
    tag "igv_session"
        publishDir "${params.outdir}/visualization", mode: 'copy'

        container 'nanopore-rnaseq:latest'

        input:
        path bed_file
        path reference
        path bam
        path bai

        output:
        path "igv_session.xml", emit: session
        path "tracks_info.txt", emit: info

        script:
        """
#!/usr/bin/env python3

# Create IGV session XML
        session_xml = '''<?xml version="1.0" encoding="UTF-8"?>
        <Session genome="mm39" hasGeneTrack="true" hasSequenceTrack="true" version="8">
        <Resources>
        <Resource path="${reference}"/>
        <Resource path="${bam}"/>
        <Resource path="${bed_file}"/>
        </Resources>
        <Panel height="500" name="Panel1" width="1000">
        <Track attributeKey="Reference sequence" clazz="org.broad.igv.track.SequenceTrack" fontSize="10" id="Reference sequence" name="Reference sequence" visible="true"/>
        <Track attributeKey="${bam}" clazz="org.broad.igv.sam.AlignmentTrack" fontSize="10" id="${bam}" name="RNA-seq Alignments" visible="true">
        <RenderOptions/>
        </Track>
        <Track attributeKey="${bed_file}" clazz="org.broad.igv.track.FeatureTrack" fontSize="10" id="${bed_file}" name="m6A Sites" visible="true">
        <RenderOptions/>
        </Track>
        </Panel>
        <PanelLayout dividerFractions="0.6,0.8"/>
        <HiddenAttributes>
        <Attribute name="DATA FILE"/>
        <Attribute name="DATA TYPE"/>
        <Attribute name="NAME"/>
        </HiddenAttributes>
        </Session>'''

        with open('igv_session.xml', 'w') as f:
        f.write(session_xml)

# Create tracks info
        with open('tracks_info.txt', 'w') as f:
        f.write("IGV Visualization Files Created\\n")
        f.write("==============================\\n")
        f.write(f"Reference: {reference}\\n")
        f.write(f"BAM file: {bam}\\n")
        f.write(f"m6A sites BED: {bed_file}\\n")
        f.write("\\nTo view in IGV:\\n")
        f.write("1. Open IGV\\n")
        f.write("2. Load session file: igv_session.xml\\n")
        f.write("3. Navigate to Htt gene region\\n")
        """
}

// Process 7: Generate summary report
process GENERATE_REPORT {
    tag "report"
        publishDir "${params.outdir}/report", mode: 'copy'

        container 'nanopore-rnaseq:latest'

        input:
        path basecalling_summary
        path modkit_summary
        path m6a_stats

        output:
        path "pipeline_report.html", emit: report

        script:
        """
#!/usr/bin/env python3

        html_content = '''<!DOCTYPE html>
        <html>
        <head>
        <title>Nanopore RNA-seq Analysis Report</title>
        <style>
        body { font-family: Arial, sans-serif; margin: 40px; }
    h1, h2 { color: #2c3e50; }
    .section { margin: 20px 0; padding: 15px; border-left: 4px solid #3498db; }
    .stats { background-color: #f8f9fa; padding: 10px; border-radius: 5px; }
    </style>
        </head>
        <body>
        <h1>Nanopore RNA-seq Analysis Report</h1>
        <div class="section">
        <h2>Pipeline Overview</h2>
        <p>This report summarizes the analysis of Nanopore direct RNA sequencing data
        following the protocol for HdhQ7/Q7 and Hdh+/Q111 mouse striatal tissue.</p>
        </div>

        <div class="section">
        <h2>Processing Steps</h2>
        <ol>
        <li>Basecalling with Dorado (${params.dorado_model})</li>
        <li>BAM to FASTQ conversion with methylation tags</li>
        <li>Alignment to mouse genome (mm39) with minimap2</li>
        <li>Methylation analysis with modkit</li>
        <li>m6A site filtering (depth ≥ ${params.min_depth})</li>
        <li>IGV visualization preparation</li>
        </ol>
        </div>

        <div class="section">
        <h2>Output Files</h2>
        <ul>
        <li><strong>alignment/</strong> - BAM files and indices</li>
        <li><strong>methylation/</strong> - bedMethyl files</li>
        <li><strong>m6a_sites/</strong> - Filtered m6A sites in BED format</li>
        <li><strong>visualization/</strong> - IGV session files</li>
        </ul>
        </div>

        <div class="section">
        <h2>Next Steps</h2>
        <p>Load the IGV session file to visualize m6A modifications along the Htt gene.
        The filtered BED files contain high-confidence m6A sites for further analysis.</p>
        </div>
        </body>
        </html>'''

        with open('pipeline_report.html', 'w') as f:
        f.write(html_content)
        """
}

// Workflow definition
workflow {
    // Input channels
    pod5_ch = Channel.fromPath("${params.pod5_dir}/*.pod5").collect()
        reference_ch = Channel.fromPath(params.reference_fa)

        // Main workflow
        if (params.bam_file) {
            // If BAM file is provided, start from alignment step
            bam_ch = Channel.fromPath(params.bam_file)
                bai_ch = Channel.fromPath("${params.bam_file}.bai")
        } else {
            // Full pipeline starting from pod5 files
            DORADO_BASECALLER(pod5_ch)
                BAM_TO_FASTQ(DORADO_BASECALLER.out.bam)
                MINIMAP2_ALIGN(BAM_TO_FASTQ.out.fastq, reference_ch)
                bam_ch = MINIMAP2_ALIGN.out.bam
                bai_ch = MINIMAP2_ALIGN.out.bai
        }

    // Methylation analysis
    MODKIT_PILEUP(bam_ch, bai_ch, reference_ch)
        FILTER_M6A_SITES(MODKIT_PILEUP.out.bedmethyl, params.min_depth)

        // Visualization
        CREATE_IGV_SESSION(
                FILTER_M6A_SITES.out.filtered_bed,
                reference_ch,
                bam_ch,
                bai_ch
                )

        // Generate report
        if (params.bam_file) {
            // Create empty summary files for BAM-only mode
            basecall_summary = Channel.empty()
        } else {
            basecall_summary = DORADO_BASECALLER.out.summary
        }

    GENERATE_REPORT(
            basecall_summary.ifEmpty(file('NO_FILE')),
            MODKIT_PILEUP.out.summary,
            FILTER_M6A_SITES.out.stats
            )
}

workflow.onComplete {
    log.info ( workflow.success ? "\nPipeline completed successfully!" : "Pipeline failed" )
}
