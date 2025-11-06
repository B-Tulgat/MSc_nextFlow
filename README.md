# ONT DRS Remora Pipeline



A concise, reproducible Nextflow pipeline that reproduces the EPI2ME basecalling workflow's m6A detection accuracy in DRACH motifs using **Dorado + Remora** and **Modkit** for validation. This repository contains the pipeline, example data downloads, visualization notebooks, and validation scripts.


<p align="center">
  <img alt="image" src="https://github.com/user-attachments/assets/7ae9d3de-c6a0-49d4-86df-34e7b3cead86" width="60%" />
</p>

---

## Quick summary

This pipeline accepts a `.pod5` file as input and produces a `pod5_bed.nf` result that contains detected m6A calls in BED-like format. The minimal flow is:

```
raw POD5 (input) -> Dorado + Remora -> BAM -> modkit pileup -> BED (output)
```

Paths used in examples throughout this README:

* `/pod5` — place or download your `.pod5` here
* `/outdir` — where `pod5_bed.nf` emits the `pod5_bed.nf` result and other outputs
* `/reference` — FASTA reference(s)

---

## Repository

```
/                       <- repo root
├─ pod5_bed.nf          <- Nextflow pipeline
├─ custom.config        <- nextflow config for running EPI2ME workflow
├─ reference/           <- place reference FASTA here for EPI2ME workflow
├─ pod5/                <- example pod5 files
├─ outdir/              <- pipeline outputs (BAMs, BEDs, logs)
├─ wf-basecalling out-dir/    <- pipeline outputs (BAMs, Execution Timeline, logs)
├─ Visualization/
│  ├─ BED_Visualization.ipynb
│  └─ Pod5_Visualization.ipynb
└─ README.md
```

---

## Requirements

Install these before running the pipeline:

* `Nextflow` (tested with `nextflow` latest stable)
* `Docker` (or modify pipeline for Singularity/Podman)
* `AWS CLI` (for fetching the example data in this README using `--no-sign-request`)
* `modkit` (for validation — see Modkit docs for install instructions)

> Tip: use a single `Dockerfile` with pinned tool versions so results are reproducible.

---

## Quickstart

Build the image used by Nextflow:

```bash
docker build -t nanopore:latest .
```

---

### Example data downloads

Download the example `.pod5` and reference used for the validation in this repo. These commands use the public ONT S3 bucket (no credentials required):

```bash
aws s3 cp s3://ont-open-data/rna-modbase-validation_2025.03/subset/m6A_DRACH.pod5 ./pod5 --no-sign-request
aws s3 cp s3://ont-open-data/rna-modbase-validation_2025.03/references/drach_context_strands.fa ./reference --no-sign-request
```

Run Nextflow (example):

```bash
nextflow run pod5_bed.nf --with-docker
```

<p align="center">
  <img alt="image" src="https://github.com/user-attachments/assets/7ae9d3de-c6a0-49d4-86df-34e7b3cead86" width="60%" />
</p>


Example `epi2me-labs/wf-basecalling` workflow (EPI2ME Nextflow workflow) that produces BAM output:

```bash
nextflow run epi2me-labs/wf-basecalling \
  -c custom.config \
  -w /tmp/work \
  --input ./pod5 \
  --ref ./reference/drach_context_strands.fa \
  --basecaller_cfg rna004_130bps_sup@v5.2.0 \
  --remora_cfg rna004_130bps_sup@v5.2.0_m6A_DRACH@v1 \
  --output_fmt bam
```

> Use `--output_fmt bam` to ensure BAM (not CRAM) is produced.

<p align="center">
  <img src="https://github.com/user-attachments/assets/e203b82f-ef5b-4a49-8501-1f5ce65c8ce2" width="70%">
</p>

## Output locations

* `./wf-basecalling out-dir/SAMPLE.pass.bam` (when using epi2me workfow)
* `outdir/` (when using the custom `pod5_bed.nf` pipeline)



---

## Running validation (Modkit)

Download ground truth BED files used in the validation:

```bash
aws s3 cp s3://ont-open-data/rna-modbase-validation_2025.03/references/drach_context_A_sites.bed ./reference --no-sign-request
aws s3 cp s3://ont-open-data/rna-modbase-validation_2025.03/references/drach_context_m6A_sites.bed ./reference --no-sign-request
```

Assuming you have `m6A_detected_sorted.bam` from the pipeline, validate with Modkit:

```bash
modkit validate \
  --bam-and-bed m6A_detected_sorted.bam drach_context_m6A_sites.bed \
  --log-filepath merged_validation.log \
  -o merged_validation_results.tsv
```

The TSV produced (`merged_validation_results.tsv`) contains the metrics compared to the ground truth.

<p align="center">
  <img src="https://github.com/user-attachments/assets/e8881086-5656-4e99-8f23-f5d3f88fb102" width="45%">
</p>


```bash
modkit validate --bam-and-bed SAMPLE.pass.bam drach_context_m6A_sites.bed
```


## Comparison between the EPI2ME-workflow and Custom-Workflow

<table>
  <tr>
    <td align="center" width="45.8%">
      <img src="https://github.com/user-attachments/assets/e8881086-5656-4e99-8f23-f5d3f88fb102" width="100%">
      <br><b>Custom 99.51% Filtered Accuracy</b>
    </td>
    <td align="center" width="50%">
      <img src="https://github.com/user-attachments/assets/ef129f97-345e-4965-865d-facbd2604db2" width="100%">
      <br><b>EPI2ME 99.55% Filtered Accuracy</b>
    </td>
  </tr>
</table>

---

## Visualization notebooks

We include two Jupyter notebooks that demonstrate how to visualize results and pod5 signals:

* `Visualization/BED_Visualization.ipynb` — shows how to render modkit BED-style outputs into PNG/figures.
* `Visualization/Pod5_Visualization.ipynb` — nanopore picoampere signal visualization + metadata.

<img width="1335" height="409" alt="image" src="https://github.com/user-attachments/assets/d3707d01-8ef8-461e-a644-93878ff40d24" />
<img width="1200" height="400" alt="image" src="https://github.com/user-attachments/assets/923c98ae-59e3-47bd-9751-27418b352780" />




---

## IGV compatibility fix

`modkit pileup` produces a bed-like file that IGV may not read directly. Convert it with this `awk` one-liner:

```bash
awk 'BEGIN{OFS="\t"} {print $1, $2, $3, $4, $11, $6}' out.bed > igv_ready.bed
```

---



These numbers were produced using the dataset in `rna-modbase-validation_2025.03` and validated with Modkit.

---



## EPI2ME Report Summary

<p align="center">
  <img src="https://github.com/user-attachments/assets/692987f1-a81c-4719-82c6-4fea933348eb" width="80%">
</p>


## Nextflow Execution Timeline

<p align="center">
  <img src="https://github.com/user-attachments/assets/ffdd7f79-9a08-4868-b240-36c81215d228" width="80%">
</p>

