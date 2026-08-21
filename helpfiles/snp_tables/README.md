# SNP tables for the SNP scan

The SNP scan (**Step 5b**) imputes each pool's ALT frequency at a SNP as

```
f_ALT(pool, SNP) = h(pool, window) · s(SNP)
```

`h` is the smoothed founder haplotype frequency vector. `s` is the vector of
per-founder ALT frequencies at that SNP, and is what this table supplies.

## Build it from your catalog

`h` is derived from `RefAlt`, which is counted against the SNP catalog. So `s`
must come from the same catalog — otherwise the two halves of that product are
different ascertainments of the same founders, and the scan tests SNPs the caller
rejected.

```bash
Rscript scripts/catalog2snptable.R \
    --catdir process/<project>/Catalog \
    --out    helpfiles/<project>_SNPs.cM.txt.gz
```

Inputs, all produced by Step 3:

| File | Supplies |
|---|---|
| `catalog.tsv.gz` | the positions — exactly the sites samples were counted at |
| `catalog.annot.tsv.gz` | per-founder read depths (`AD_<founder>`) |
| `catalog.founders.txt` | founder names, in catalog column order |

`cM` is interpolated from `helpfiles/flymap.r6.txt` using the same smoothing the
scans use. Output schema:

```
CHROM  POS  <founder1> ... <founderN>  cM
```

The script prints the `--snp-table` and `--founders` arguments to pass on. There
is no A/B split step: a catalog is built for one population, so its table already
carries the right founder set.

**Rebuild after any catalog change.** A re-cut with different `--min-dp`,
`--maxaf` or `--snpgap` changes which positions are in `catalog.tsv.gz`, and the
SNP table must change with it. It is cheap — no founder recall.

## Using it

```bash
bash scripts/run_snp_scan.sh \
    --design    helpfiles/<project>/design.txt \
    --dir       process/<project> \
    --scan      <scan_name> \
    --snp-table helpfiles/<project>_SNPs.cM.txt.gz \
    --founders  B1,B2,B3,B4,B5,B6,B7,AB8
```

`--founders` must list the same founders in the same order as the table's
columns. Use the same `--scan` name as Step 5a; `run_scan.sh` must have finished
first. `run_snp_scan.sh` submits the SLURM array itself — do not call
`snp_scan.sh` directly.

## Older projects: FREQ_SNPs_*.cM.txt.gz

Runs predating the founder-catalog caller used
`FREQ_SNPs_{A,B}pop.cM.txt.gz`, split out of a whole-genome founder table by
`scripts/prep_snp_table.R`. That table was not produced by this pipeline, its
filters are not recorded, and it is no longer tracked in the repo (commit
`3be8127`). Measured against the catalog's own founder filters, about 24% of its
SNPs are non-segregating among the founders and ~5% have a founder at
intermediate frequency. The scan drops non-segregating sites at run time
(`snp_scan.R:192`), so old results are not wrong on that account — but the SNP
set is not the catalog's.

Keep these files only to reproduce a pre-catalog analysis. For anything new, use
`catalog2snptable.R`.
