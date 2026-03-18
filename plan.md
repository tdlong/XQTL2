# XQTL2 Pipeline Modernization Plan

---

## What Is Already Done

- [x] `.gitignore` expanded (excludes `*.tar`, `*.bam`, `.RData`, `*.out`, etc.)
- [x] `scripts/reorganize_cluster.sh` + fix scripts — one-offs moved to `scripts_oneoffs/`
- [x] `haps2scan.stable.*` and `haps2scan.freqsmooth.*` tracked in `scripts/`
- [x] Bug fix: `nrepl=1` crash in `wald.test3` / `mn.covmat`

---

## Pipeline Architecture (Two Parallel Pipelines)

Both pipelines are valid and should remain available. Users choose which to run.

```
fq2bam.sh  →  bam2bcf2REFALT.sh  →  REFALT2haps.R
                                          │
              ┌───────────────────────────┴────────────────────────────┐
              │                                                        │
   Pipeline 1: haps2scan.Apr2025.sh              Pipeline 2: haps2scan.freqsmooth.sh
   (raw scan, no smoothing)                      (covariance + optional freq smoothing)
              │                                                        │
              └───────────────────────────┬────────────────────────────┘
                                          │
                          concat_Chromosome_Scans.sh
                                          │
                               configs/<experiment>.R   ← experiment-specific plots
                                          │
                                    figures/slides/
```

**Experiment-specific submit scripts** (in `scripts_oneoffs/`) orchestrate the full
chain: submit scan arrays → depend concat → depend plots. `run_all_scans_125.sh` is
the model for this pattern.

---

## Pipeline 2 Details (freqsmooth / stable)

`haps2scan.freqsmooth.sh` calls `haps2scan.freqsmooth.R` which sources
`haps2scan.stable.code.R`. This code does:

1. **Covariance matrix smoothing** — running mean of reconstruction error covariance
   matrices across ±`SMOOTH_HALF` windows. Currently **hardcoded** to 50 (= ±250 kb
   at 5 kb/window).
2. **Eigenvalue regularization** — floors small eigenvalues at `RIDGE_FRACTION * mean`.
   Currently hardcoded to 0.01. Fine to leave as a documented constant.
3. **Frequency vector smoothing** (optional) — running mean of haplotype frequencies
   across ±`FREQ_SMOOTH_HALF` windows. Already an explicit argument. Preferred value:
   25 (= ±125 kb at 5 kb/window) for most experiments; 13 for 10 kb/window experiments.

### The only remaining problem: `SMOOTH_HALF` is hardcoded

`FREQ_SMOOTH_HALF` is already passed cleanly through `.sh` → `.R` → `.code.R`.
`SMOOTH_HALF` is not — it is hardcoded at 50 in `stable.code.R`. This needs to become
an explicit argument using the same pattern.

---

## Remaining Work

### Step 1 — Make `SMOOTH_HALF` an explicit argument

**Pattern** (mirror how `FREQ_SMOOTH_HALF` is already handled):

`haps2scan.freqsmooth.sh`:
```bash
COV_SMOOTH_HALF=$5
FREQ_SMOOTH_HALF=$6
Rscript scripts/haps2scan.freqsmooth.R $mychr $Rfile $mydir $myoutdir $COV_SMOOTH_HALF $FREQ_SMOOTH_HALF
```

`haps2scan.freqsmooth.R`:
```r
COV_SMOOTH_HALF  <- as.integer(args[5])
FREQ_SMOOTH_HALF <- as.integer(args[6])
```

`haps2scan.stable.code.R` — remove the hardcoded line:
```r
# REMOVE: SMOOTH_HALF <- 50
# It is now set by the driver .R file before this file is sourced.
if (!exists("SMOOTH_HALF"))      SMOOTH_HALF      <- 0L   # fallback: no cov smoothing
if (!exists("FREQ_SMOOTH_HALF")) FREQ_SMOOTH_HALF <- 0L   # fallback: no freq smoothing
```

**Note on units:** Window counts are used internally because the conversion depends on
the experiment's window step size (5 kb or 10 kb). The experiment-specific submit scripts
(in `scripts_oneoffs/`) are where kb → window-count conversion is documented. This is
appropriate — `run_all_scans_125.sh` already does this clearly in its header comment.

### Step 2 — Rename `pseudoscan` → `scan` in output filenames

Simpler. Apply in `haps2scan.stable.code.R` (and `haps2scan.Apr2025.code.R` for
consistency):

```r
fileout              <- paste0(mydirout, "/", myoutdir, ".scan.",          mychr, ".txt")
fileout_meansBySample<- paste0(mydirout, "/", myoutdir, ".meansBySample.", mychr, ".txt")
```

No need to encode smoothing parameters in the filename — the output directory name
already carries that information (e.g., `ZINC2_M_freqs125/`), following the pattern
established in `run_all_scans_125.sh`.

### Step 3 — `configs/` plotting scripts

The plotting scripts (`zinc2_125.R`, `aging_125.R`, etc.) live in `configs/` and are
called from experiment-specific submit scripts. These are experiment-specific and do
not need to be part of the core `scripts/` pipeline.

**Open question:** Should example/template plotting scripts in `configs/` be tracked?
They document the expected output format and serve as templates for new experiments.
Recommendation: track at least one as a template (e.g., a generic `configs/example_scan_plot.R`).

### Step 4 — Update `README.md`

Add section "Running the scan" with two subsections:

**Pipeline 1 (no smoothing):**
```bash
sbatch --array=1-5 scripts/haps2scan.Apr2025.sh $Rfile $mydir $myoutdir
```

**Pipeline 2 (stabilized, with smoothing):**
```bash
# COV_SMOOTH_HALF: ±windows for covariance smoothing (50 = ±250kb at 5kb/window)
# FREQ_SMOOTH_HALF: ±windows for frequency smoothing (25 = ±125kb at 5kb/window, 0 = off)
sbatch --array=1-5 scripts/haps2scan.freqsmooth.sh $Rfile $mydir $myoutdir $COV_SMOOTH_HALF $FREQ_SMOOTH_HALF
```

Then concatenate and plot following the pattern in `scripts_oneoffs/run_all_scans_125.sh`.

---

## Implementation Order

- [x] 1. Cleanup, gitignore, track stable/freqsmooth scripts
- [ ] 2. Add `SMOOTH_HALF` as explicit arg to `.sh` → `.R` → `.code.R`
- [ ] 3. Remove hardcoded `SMOOTH_HALF <- 50` from `stable.code.R`; add `if (!exists(...))` fallback
- [ ] 4. Rename `pseudoscan` → `scan` in output filenames (both pipelines)
- [ ] 5. Resolve `configs/` tracking question; add example plot template if appropriate
- [ ] 6. Update `README.md`
- [ ] 7. Push; pull on cluster and test

---

## Open Questions

1. Should any `configs/` plotting scripts be tracked as templates?
2. Should `haps2scan.stable.sh` (no freq smoothing) be kept, or is `haps2scan.freqsmooth.sh`
   with `FREQ_SMOOTH_HALF=0` sufficient?
3. `concat_Chromosome_Scans.Andreas.sh` — is this still the right concat script, or has
   it been updated? Should it be renamed to drop the `Andreas` label?
