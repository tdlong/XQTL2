# XQTL2 Frequency-Smoothed Pipeline — Implementation Plan

**Author:** Anthony Long + Claude
**Date:** 2026-03-21
**Status:** Planning — do not modify production code

---

## 1. Motivation

The current `haps2scan` pipeline applies smoothing inconsistently:

| Data | Smoothed? | Notes |
|------|-----------|-------|
| Covariance matrices | Always (±250 kb) | Applied before Wald test |
| Haplotype frequencies (for Wald test) | Optional flag | Only sometimes |
| Heritability inputs | Never | Always uses raw frequencies |
| `meansBySample` output | Never | Always raw frequencies |

The proposed pipeline fixes this by smoothing haplotype frequencies **first** — before any statistics are computed — so that the Wald test, heritability estimates, and means output all operate on the same smoothed frequencies. This also naturally enables SNP-based QTL scans (Step 3) using the same smoothed frequency estimates.

---

## 2. What Will NOT Change

- The XQTL2 R package (visualization functions) is **untouched**
- Output file formats (`scan.txt`, `meansBySample.txt`) are **identical** — all downstream plotting works without modification
- The core `scan_functions.R` (Wald test, heritability) is **read-only referenced**, not modified
- No changes to the current production `haps2scan.stable.code.R`

---

## 3. File Structure

New scripts live in `scripts_freqsmooth/` inside the XQTL2 project, alongside the existing `scripts/` directory. The one-off test pipeline driver goes in `scripts_oneoffs/`, consistent with the project convention. Nothing in `scripts/` is touched.

```
XQTL2/
├── scripts/                           # UNCHANGED — current production scripts
├── scripts_freqsmooth/                # NEW
│   ├── smooth_haps.R                  # Reads haps out RDS → smoothed out RDS + smoothed means
│   ├── smooth_haps.sh                 # SLURM array wrapper
│   ├── freqsmooth_scan.R              # Reads smoothed out RDS → scan file (Wald + heritability)
│   ├── freqsmooth_scan.sh             # SLURM array wrapper
│   ├── snp_scan.R                     # Reads smoothed out RDS + SNP table → SNP scan file
│   └── snp_scan.sh                    # SLURM array wrapper
└── scripts_oneoffs/
    └── malathion_test_v2.sh           # Test driver (replaces steps 5+ of malathion_test_pipeline.sh)
```

**Input file**: `R.haps.<chr>.out.rds` — haplotype estimates produced by Step 4. Confirmed as the input to the existing `haps2scan.freqsmooth.R` (line 13). It is a tibble with one row per genomic window, and list-columns `CHROM`, `pos`, `sample` (pool names), `Haps` (founder frequency vector per pool), `Err` (reconstruction covariance matrix per pool), `Names` (founder names), `Groups`.

**Output files** — exact same names and formats as the current pipeline, so all downstream steps (concat, plot) work without modification:

| File | Description |
|------|-------------|
| `<scan_name>.meansBySample.<chr>.txt` | Per-founder smoothed frequency table — used by `XQTL_change_average()`, `XQTL_change_byRep()`, etc. |
| `<scan_name>.scan.<chr>.txt` | Genome scan statistics — used by `XQTL_Manhattan()`, `XQTL_zoom()`, `XQTL_region()`, etc. |
| `<scan_name>.snp_scan.<chr>.txt` | SNP-level scan — same column structure as scan file; same viz functions apply |

After `concat_Chromosome_Scans.sh`: `<scan_name>.scan.txt`, `<scan_name>.meansBySample.txt`, `<scan_name>.snp_scan.txt`.

Dependencies (read-only):
- `scripts/scan_functions.R` — Wald test and heritability functions
- `FREQ_SNPs.cM.txt.gz` — founder SNP state table (`snp_scan.R` only)

---

## 4. Shared Parameters

All scripts use the same named command-line arguments, matching the conventions already established in `malathion_test_pipeline.sh`:

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--chr` | string | required | Chromosome (chrX, chr2L, etc.) |
| `--dir` | path | required | Directory containing `R.haps.chrX.rds` files |
| `--outdir` | path | required | Output subdirectory name (created under `--dir`) |
| `--rfile` | path | required | Pool design file (same format as current pipeline) |
| `--smooth-kb` | integer | 125 | Half-window for smoothing in kb — **always expressed in kb, never in window counts** |

`snp_scan.sh` additionally requires:

| Parameter | Type | Description |
|-----------|------|-------------|
| `--snp-table` | path | Path to `FREQ_SNPs.cM.txt.gz` |
| `--founders` | string | Comma-separated founder column names (e.g. `A1,A2,A3,A4,A5,A6,A7,AB8`) |

The conversion from kb to window counts is an **internal implementation detail** handled transparently inside each R script: the step size in bp is inferred from the haplotype position data itself, so no step parameter is exposed or required.

---

## 5. Step 1 — Smooth Haplotype Frequencies

### Input
`R.haps.<chr>.out.rds` — haplotype frequency estimates per pool, one file per chromosome, containing per-pool haplotype frequency estimates and reconstruction covariance matrices at ~5 kb windows.

### What it does

For every pool × founder combination across all windows on a chromosome, apply an **edge-aware running mean** with half-window = `SMOOTH_HALF` positions:

```r
running_mean <- function(x, half_window) {
  n       <- length(x)
  is_val  <- !is.na(x)
  x_clean <- replace(x, !is_val, 0)
  cs      <- c(0, cumsum(x_clean))
  cn      <- c(0, cumsum(as.numeric(is_val)))
  lo      <- pmax(1L, seq_len(n) - half_window)
  hi      <- pmin(n,  seq_len(n) + half_window)
  totals  <- cs[hi + 1L] - cs[lo]
  counts  <- cn[hi + 1L] - cn[lo]
  ifelse(counts > 0, totals / counts, NA_real_)
}
```

Applied to:
1. **Haplotype frequencies** — each founder's frequency vector smoothed independently per pool. Because founder frequencies sum to 1 at each window and running mean is a convex combination, the smoothed frequencies still sum to 1 at every window.
2. **Reconstruction covariance matrices** — each (i,j) element of the nF × nF covariance matrix smoothed with the same half-window. This is justified by the same biological assumption: if true frequencies vary smoothly, so does estimation uncertainty.

### Outputs
- `<scan_name>.smoothed.out.<chr>.rds` — same structure as input, frequencies and covariances replaced with smoothed versions; input to Steps 2 and 3
- `<scan_name>.meansBySample.<chr>.txt` — long-format means table (chr, pos, TRT, REP, founder, freq) from smoothed frequencies; same format as current pipeline output, usable directly by all existing viz functions

---

## 6. Step 2 — Scan on Smoothed Frequencies

### Input
`<scan_name>.smoothed.out.<chr>.rds` from Step 1.

### What it does

Nearly identical logic to `haps2scan.stable.code.R`, with these changes:

| Current behavior | New behavior | Reason |
|-----------------|--------------|--------|
| `FREQ_SMOOTH_HALF` = optional | Not needed (set to 0 internally) | Smoothing already done in Step 1 |
| `SMOOTH_HALF = 50` for covariances | Not needed (set to 0 internally) | Covariance smoothing already done in Step 1 |
| Heritability from raw frequencies | Heritability from smoothed frequencies | Consistency — same data as Wald test |
| Means output = raw frequencies | Means output = smoothed frequencies | Consistency — Step 1 already wrote this |

The Wald test logic (eigendecomposition, ridge regularization, effective N, replicate pooling) is unchanged.

### Output
`<scan_name>.scan.<chr>.txt` — same columns as current pipeline:
`chr, pos, Wald_log10p, Falc_H2, Cutl_H2, cM`

---

## 7. Step 3 — SNP-Based QTL Scan

### Input
- `<scan_name>.smoothed.out.<chr>.rds` from Step 1 (smoothed haplotype frequencies and covariances)
- `FREQ_SNPs.cM.txt.gz` — founder SNP state table

### 7.1 Founder SNP State Table Format

`FREQ_SNPs.cM.txt.gz` columns (47 named columns + leading row index):

```
[row_idx]  CHROM  POS  A1 A2 A3 A4 A5 A6 A7 AB8  B1 B2 B3 B4 B5 B6 B7  C01 ... Z12  cM
```

- **A-population founders**: A1, A2, A3, A4, A5, A6, A7, AB8 (8 founders)
- **B-population founders**: B1, B2, B3, B4, B5, B6, B7 (7 founders)
- **Pool columns**: C01–C12B (control), Z01–Z12 (selected) — **not used** in Step 3
- **Founder values**: ALT allele frequency for that founder at that SNP, approximately 0 or 1 for inbred lines

For the **malathion test**, we extract columns `A1, A2, A3, A4, A5, A6, A7, AB8` (A population, 8 founders, as specified in `hap_params.R`). The B-population columns (B1–B7) are available for other experiments. The `--founders` argument selects which columns to use.

The file contains **1,675,709 SNPs** across all chromosomes — all pre-filtered to positions worth testing. No additional density filter is applied; we use all rows on the target chromosome.

### 7.2 Converting Smoothed Haplotype Frequencies to SNP Frequencies

For a given pool and SNP, the conversion is a dot product: **the pool's ALT frequency at a SNP equals the sum of each founder's contribution, weighted by how much of that founder is present in the pool.**

Concretely, for pool `p` and SNP at position `t`:

**Step A — Find the nearest haplotype window.**
Each SNP in `FREQ_SNPs.cM.txt.gz` is at a specific base-pair position. The smoothed haplotype estimates exist at ~5 kb window centers. For each SNP we find `w*(t)` = the haplotype window with the closest position to `t` (a nearest-join on position).

**Step B — Compute the imputed ALT frequency.**

Let `h(p, f)` = smoothed haplotype frequency of founder `f` in pool `p` at window `w*(t)` (8 values summing to 1, from the smoothed `Haps` column of the RDS).

Let `s(f, t)` = the ALT allele state of founder `f` at SNP `t` (approximately 0 or 1, from the A1–AB8 columns of `FREQ_SNPs.cM.txt.gz`).

Then:

```
f̂_ALT(p, t) = h(p, A1) × s(A1, t)  +  h(p, A2) × s(A2, t)  +  ···  +  h(p, AB8) × s(AB8, t)
             = h(p, w*(t))ᵀ · s(t)
```

`f̂_REF = 1 − f̂_ALT`. Clip to [0, 1] to guard against floating-point edge cases.

### 7.3 Propagating Reconstruction Covariance to SNP Space

The smoothed haplotype reconstruction covariance `Σ_h(p, w)` is `n_F × n_F`. This uncertainty propagates to the biallelic (REF, ALT) space through the `2 × n_F` founder state matrix **S**(t):

```
S(t) = | 1−s_1  1−s_2  ···  1−s_{n_F} |   ← REF row
       |   s_1    s_2  ···    s_{n_F}  |   ← ALT row
```

The propagated 2×2 reconstruction covariance is:

```
Σ_rec_SNP(p, t) = S(t) × Σ_h(p, w*(t)) × S(t)ᵀ
```

This is the reconstruction-only uncertainty (haplotype calling error mapped to SNP space). The multinomial sampling variance is added internally by `wald.test3`.

### 7.4 Modified Wald Test (df = 1)

With biallelic frequency vectors of length 2 (`[REF_freq, ALT_freq]`) and 2×2 covariance matrices, the existing `wald.test3` function is called directly:

```r
result <- wald.test3(
  p1    = cbind(f̂_REF_C, f̂_ALT_C),   # nrep × 2 matrix, control pools
  p2    = cbind(f̂_REF_Z, f̂_ALT_Z),   # nrep × 2 matrix, selected pools
  covar1 = Σ_rec_SNP_C,               # 2×2×nrep array, propagated covariances
  covar2 = Σ_rec_SNP_Z,
  nrepl  = nrep,
  N1     = N_control,
  N2     = N_selected
)
```

Because `length(p1[1,]) = 2`, `wald.test3` automatically uses `df = 2 − 1 = 1`. The eigendecomposition of the 2×2 total covariance (which has rank 1 due to the sum-to-1 constraint) yields one non-trivial eigenvalue. The test statistic is chi-squared with df = 1.

**No modifications to `wald.test3` are needed.**

### 7.5 Heritability at SNP Positions

Both Falconer H² and Cutler H² are computed per SNP using the biallelic adaptation of the same formulas applied in Step 2. The ALT allele plays the role of the "focal allele"; REF is redundant (1 − ALT) and omitted.

**Falconer H² (per SNP):**

For replicate `r`, let `f_C(r)` = ALT frequency in control pool, `f_Z(r)` = ALT frequency in selected pool:

```
H2_temp(r) = (f_Z(r) − f_C(r))² / max(f_C(r), af_cutoff)
Falc_H2 = 200 × mean_r(H2_temp) / i²
```

where `i = dnorm(qnorm(1 − Proportion)) / Proportion` (intensity of selection).

**Cutler H² (per SNP):**

```
Penetrance(r) = clamp( f_Z(r) × Proportion / f_C(r),  Proportion/2,  2×Proportion )
Affect(r)     = qnorm(1 − Proportion) − qnorm(1 − Penetrance(r))
Cutl_H2       = 200 × sum_r( Affect(r)² × f_C(r) )
```

These are identical in structure to the haplotype-scan formulas; the only difference is summing over 1 allele instead of 8 founders.

### 7.6 Output

`<scan_name>.snp_scan.<chr>.txt` — same core columns as scan file, plus one informational column:

```
chr  pos  Wald_log10p  Falc_H2  Cutl_H2  cM  n_informative_founders
```

`n_informative_founders` = number of founders that differ at this SNP. A SNP where all founders have the same state (no contrast possible) is skipped and omitted from output.

Downstream Manhattan plots work unchanged — they use `chr`, `pos`, `Wald_log10p`.

---

## 8. Test Pipeline — Malathion Data

`malathion_test_v2.sh` takes the haplotype RDS files produced by the existing pipeline (steps 1–4 unchanged) and runs the new steps 5a–5c in their place.

Design file format confirmed (`helpfiles/malathion_test/design.txt`):
```
bam      TRT  REP  REPrep  Num  Proportion
s_F_1    Z    1    1       500  0.01
c_F_1    C    1    1       500  NA
s_M_1    Z    2    1       500  0.01
c_M_1    C    2    1       500  NA
```
`Proportion` is present → Falconer and Cutler H² can be computed.

```bash
PROJECT=malathion_test
DESIGN=helpfiles/${PROJECT}/design.txt
SCAN=MALATHION_TEST_v2_smooth125
SMOOTH_KB=125
SNP_TABLE=FREQ_SNPs.cM.txt.gz
FOUNDERS=A1,A2,A3,A4,A5,A6,A7,AB8

# Step 5a: smooth haplotype frequencies → smoothed RDS + smoothed means file
jid_smooth=$(sbatch --parsable --dependency=afterok:${jid_haps} \
    --array=1-5 scripts_freqsmooth/smooth_haps.sh \
    --rfile ${DESIGN} --dir process/${PROJECT} --outdir ${SCAN} \
    --smooth-kb ${SMOOTH_KB})

# Step 5b: Wald scan on smoothed frequencies → scan file
jid_scan=$(sbatch --parsable --dependency=afterok:${jid_smooth} \
    --array=1-5 scripts_freqsmooth/freqsmooth_scan.sh \
    --rfile ${DESIGN} --dir process/${PROJECT}/${SCAN} --outdir ${SCAN})

# Step 5c: SNP scan → snp_scan file (can run in parallel with 5b)
jid_snp=$(sbatch --parsable --dependency=afterok:${jid_smooth} \
    --array=1-5 scripts_freqsmooth/snp_scan.sh \
    --rfile ${DESIGN} --dir process/${PROJECT}/${SCAN} --outdir ${SCAN} \
    --snp-table ${SNP_TABLE} --founders ${FOUNDERS})

# Step 6: concat all outputs across chromosomes
jid_concat=$(sbatch --parsable --dependency=afterok:${jid_scan}:${jid_snp} \
    -A tdlong_lab -p standard --mem=10G \
    --wrap="bash scripts/concat_Chromosome_Scans.sh process/${PROJECT}/${SCAN}")
```

Steps 5b and 5c both depend on 5a but are independent of each other, so they run in parallel.

---

## 9. Validation Checklist

After running the test pipeline on malathion data, compare against the current pipeline output:

- [ ] Manhattan plots show the same major peaks (malathion resistance locus should be prominent)
- [ ] -log10p values are comparable in magnitude (not dramatically inflated or deflated)
- [ ] Heritability estimates are in expected range
- [ ] `meansBySample` frequencies are smooth (visual check: no jagged rep-to-rep oscillation)
- [ ] SNP scan peaks overlap with haplotype scan peaks at the malathion locus
- [ ] Runtime is acceptable (Step 1 should be fast — O(n) per chromosome)

---

## 10. Integration Path

If the malathion test validates the new approach, the intent is:

1. **Replace** `scripts/haps2scan.freqsmooth.sh` + `scripts/haps2scan.freqsmooth.code.R` with the new three-script pipeline (`smooth_haps`, `freqsmooth_scan`, `snp_scan`)
2. **Move** the scripts from `scripts_freqsmooth/` into `scripts/` as the canonical implementation
3. **Update the README** to describe the new two-step workflow (smooth first, then scan + SNP scan) and retire the legacy `haps2scan.Apr2025` option
4. The SNP scan becomes a new standard Step 5c in the README alongside the haplotype scan

The new pipeline is a strict superset of the old smooth scan: it produces the same `.scan.txt` and `.meansBySample.txt` files (now consistently using smoothed frequencies throughout), plus a new `.snp_scan.txt`. All existing downstream steps and plotting functions continue to work unchanged.

## 11. Resolved Decisions

| Question | Decision |
|----------|----------|
| Founders for malathion | A population: A1, A2, A3, A4, A5, A6, A7, AB8 (from `hap_params.R`) |
| SNP density / filtering | Use all 1.67M SNPs in `FREQ_SNPs.cM.txt.gz`; skip only monomorphic-among-founders SNPs |
| Heritability for SNP scan | Both Falconer and Cutler H², biallelic adaptation (Section 7.5) |
| `design.txt` `Proportion` column | Confirmed present; heritability is computable for malathion test |
| Haplotype input file | `R.haps.<chr>.out.rds` (haplotype estimates, not SNP table) |
| Output file naming | `.scan.<chr>.txt`, `.meansBySample.<chr>.txt`, `.snp_scan.<chr>.txt` |
| New code location | `scripts_freqsmooth/` inside XQTL2 project; test driver in `scripts_oneoffs/` |
