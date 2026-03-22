# SNP Tables for the SNP Scan

The SNP scan (Step 5c) needs a table of per-founder allele frequencies at every
SNP. This is a **one-time** preparation step per population — you do not need to
repeat it for each experiment.

## Source file

`FREQ_SNPs.cM.txt.gz` — a gzipped tab-delimited table with one row per SNP.
Columns include chromosome, position, genetic map position (cM), and allele
frequencies for all 15 founders (A1–A7, AB8, B1–B7). This file is produced
from whole-genome sequencing of the founder lines and lives in `helpfiles/`.

## Splitting by population

The A-population and B-population use different founder sets, so each needs its
own table. `scripts/prep_snp_table.R` extracts the relevant columns:

```bash
Rscript scripts/prep_snp_table.R helpfiles/FREQ_SNPs.cM.txt.gz
```

This writes two files next to the input:

| File | Founders |
|------|----------|
| `helpfiles/FREQ_SNPs_Apop.cM.txt.gz` | A1 A2 A3 A4 A5 A6 A7 AB8 |
| `helpfiles/FREQ_SNPs_Bpop.cM.txt.gz` | AB8 B1 B2 B3 B4 B5 B6 B7 |

AB8 is shared between both populations and appears in both files.

## Using the SNP table in the pipeline

Pass the appropriate population file to `snp_scan.sh` with `--snp-table`:

```bash
sbatch --array=1-5 scripts/snp_scan.sh \
    --rfile     helpfiles/<project>/design.txt \
    --dir       process/<project> \
    --outdir    <scan_name> \
    --snp-table helpfiles/FREQ_SNPs_Apop.cM.txt.gz \
    --founders  A1,A2,A3,A4,A5,A6,A7,AB8
```

The `--founders` argument must list the same founders in the same order as the
columns in the SNP table (excluding CHROM, POS, and cM).
