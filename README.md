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
