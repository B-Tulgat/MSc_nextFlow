# ONT DRS Remora Pipeline Setup and Validation

This guide demonstrates how our **custom Nextflow pipeline** reproduce the accuracy rate reported running the **official EPI2ME basecalling workflow** for N6-Methyladenosine (m6A) detection on DRACH motifs within oligonucleotides using **Modkit**.

---

## ⚙️ Setup

### Requirements

Make sure the following tools are installed:

* [AWS CLI](https://docs.aws.amazon.com/cli/latest/userguide/getting-started-install.html)
* [Docker](https://docs.docker.com/engine/install)
* [Nextflow](https://www.nextflow.io/docs/latest/install.html)

---

### Build the Docker Image

```bash
docker build -t nanopore:latest .
```

---

### Download Example Data

#### Example POD5 file

```bash
aws s3 cp s3://ont-open-data/rna-modbase-validation_2025.03/subset/m6A_DRACH.pod5 . --no-sign-request
```

#### Reference file

```bash
aws s3 cp s3://ont-open-data/rna-modbase-validation_2025.03/references/drach_context_strands.fa . --no-sign-request
```

> **Note:**
> For EPI2ME workflows, place the reference file inside the `reference/` directory.
> For the **custom workflow**, proceed directly with the above download command.

---

## Running the pipeline

```bash
nextflow run pod5_bed.nf --with-docker
```

<p align="center">
  <img src="https://github.com/user-attachments/assets/68d822af-523b-4890-8408-eff29babf3e0" width="42.5%">
  <img src="https://github.com/user-attachments/assets/1a4bd377-4275-4ecf-b9f2-6e2410133310" width="45%">
</p>

---

## Modkit Validation

### Download Ground Truth Files

```bash
aws s3 cp s3://ont-open-data/rna-modbase-validation_2025.03/references/drach_context_A_sites.bed . --no-sign-request
aws s3 cp s3://ont-open-data/rna-modbase-validation_2025.03/references/drach_context_m6A_sites.bed . --no-sign-request
```

### Validate Results

We have the resulting BAM file `m6A_detected_sorted.bam` from `m6A_DRACH.pod5`.

```bash
modkit validate \
  --bam-and-bed m6A_detected_sorted.bam drach_context_m6A_sites.bed \
  --log-filepath merged_validation.log \
  -o merged_validation_results.tsv
```

<p align="center">
  <img src="https://github.com/user-attachments/assets/e8881086-5656-4e99-8f23-f5d3f88fb102" width="45%">
</p>

---

## EPI2ME Nextflow Workflow: `epi2me-labs/wf-basecalling`

In the `--input` argument, specify the directory containing your `.pod5` files.
Each `.pod5` file will be processed separately.

```bash
nextflow run epi2me-labs/wf-basecalling \
  -c custom.config \
  -w /tmp/work \
  --input ./ \
  --ref ./reference/drach_context_strands.fa \
  --basecaller_cfg rna004_130bps_sup@v5.2.0 \
  --remora_cfg rna004_130bps_sup@v5.2.0_m6A_DRACH@v1 \
  --output_fmt bam
```

> `--output_fmt bam` ensures BAM output (instead of the default CRAM).

<p align="center">

  <img src="https://github.com/user-attachments/assets/e203b82f-ef5b-4a49-8501-1f5ce65c8ce2" width="70%">
</p>

After successful execution, the resulting file `SAMPLE.pass.bam` will be located in `./output/`.

Validate the basecall:

```bash
modkit validate --bam-and-bed SAMPLE.pass.bam drach_context_m6A_sites.bed
```


## Comparison between the EPI2ME-workflow and Custom-Workflow

<p align="center">
  <figure style="display:inline-block; text-align:center; margin:10px; width:41%;">
    <img src="https://github.com/user-attachments/assets/e8881086-5656-4e99-8f23-f5d3f88fb102" width="100%">
    <figcaption><b>Custom 99.51% Filtered Accuracy</b></figcaption>
  </figure>
  <figure style="display:inline-block; text-align:center; margin:10px; width:45%;">
    <img src="https://github.com/user-attachments/assets/ef129f97-345e-4965-865d-facbd2604db2" width="100%">
    <figcaption><b>EPI2ME 99.55% Filtered Accuracy</b></figcaption>
  </figure>
</p>
---

## EPI2ME Report Summary

<p align="center">
  <img src="https://github.com/user-attachments/assets/692987f1-a81c-4719-82c6-4fea933348eb" width="90%">
</p>

---

## Nextflow Execution Timeline

<p align="center">
  <img src="https://github.com/user-attachments/assets/ffdd7f79-9a08-4868-b240-36c81215d228" width="90%">
</p>

---

## IGV Compatibility Fix

The BED file from `modkit pileup` is **not directly compatible** with IGV.
Run the following command to convert it into a compatible format:

```bash
awk 'BEGIN{OFS="\t"} {print $1, $2, $3, $4, $11, $6}' out.bed > igv_ready.bed
```

---

## 📈 Results Summary

We successfully **reproduced the official EPI2ME results** using a simplified custom workflow, achieving comparable performance:

| Workflow        | Accuracy (%) |
| --------------- | ------------ |
| Official EPI2ME | **99.55**    |
| Custom Workflow | **99.51**    |

---

