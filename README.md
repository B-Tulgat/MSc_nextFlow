Download example pod5 file from [https://epi2me.nanoporetech.com/rna-mod-validation-data] with the following command:
```bash
aws s3 cp  s3://ont-open-data/rna-modbase-validation_2025.03/subset/m6A_DRACH.pod5 . --no-sign-request
```
![image](https://github.com/user-attachments/assets/68d822af-523b-4890-8408-eff29babf3e0)
![image](https://github.com/user-attachments/assets/1a4bd377-4275-4ecf-b9f2-6e2410133310)

Confusion matrix:

Download the necessary reference ground truth files:
```
aws s3 cp s3://ont-open-data/rna-modbase-validation_2025.03/references/drach_context_A_sites.bed . --no-sign-request
aws s3 cp s3://ont-open-data/rna-modbase-validation_2025.03/references/drach_context_m6A_sites.bed . --no-sign-request
```

We have resulting `m6A_detected.bam` from `m6A_DRACH.pod5`.

```bash
modkit validate \
  --bam-and-bed m6A_detected_sorted.bam drach_context_m6A_sites.bed \
  --log-filepath merged_validation.log \
  -o merged_validation_results.tsv
```

<img width="704" height="685" alt="image" src="https://github.com/user-attachments/assets/e8881086-5656-4e99-8f23-f5d3f88fb102" />

## EPI2ME Nextflow Workflow: epi2me-labs/wf-basecalling

In `--input` is the directory where your `pod5` file is located. The code works by filter all the files with the extension `.pod5 `. Seperate the `pod5` files for seperate basecalling.

```
nextflow run epi2me-labs/wf-basecalling \
  -c custom.config \
  -w /tmp/work \
  --input ./ \
  --ref ./reference/drach_context_strands.fa \
  --basecaller_cfg rna004_130bps_sup@v5.2.0 \
  --remora_cfg rna004_130bps_sup@v5.2.0_m6A_DRACH@v1 \
  --output_fmt bam
```
`--output_fmt bam` flag directs the command to make `bam` file rather than `cram` file by default. Which is alternative to `bam`.
<img width="918" height="779" alt="image" src="https://github.com/user-attachments/assets/e203b82f-ef5b-4a49-8501-1f5ce65c8ce2" />

After the nextflow is run successfully resulting `SAMPLE.pass.bam` will be in the `./output` directory by default. Validate the basecall by 
```
modkit validate --bam-and-bed SAMPLE.pass.bam drach_context_m6A_sites.bed
```

<img width="926" height="820" alt="image" src="https://github.com/user-attachments/assets/ef129f97-345e-4965-865d-facbd2604db2" />

EPI2ME report summary:

<img width="1591" height="895" alt="image" src="https://github.com/user-attachments/assets/692987f1-a81c-4719-82c6-4fea933348eb" />


Nextflow io provides execution timeline breakdown:

<img width="1888" height="788" alt="image" src="https://github.com/user-attachments/assets/ffdd7f79-9a08-4868-b240-36c81215d228" />


