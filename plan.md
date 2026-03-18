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

Experiment-specific submit scripts (in `scripts_oneoffs/`) orchestrate the full chain.
`run_all_scans_125.sh` is the model for this pattern.

---

## Remaining Work

### Step 1 — All smoothing parameters in kb

Users always specify smoothing in **kb**. Window counts are never exposed.
The code derives the window step size from the actual position spacing in the data,
then converts kb → windows internally.

**Submit call:**
```bash
sbatch --array=1-5 scripts/haps2scan.freqsmooth.sh \
    $Rfile $mydir $myoutdir $COV_SMOOTH_KB $FREQ_SMOOTH_KB
```

| Argument | Meaning | Default for new runs |
|----------|---------|----------------------|
| `COV_SMOOTH_KB`  | ±half-window for covariance smoothing | `125` |
| `FREQ_SMOOTH_KB` | ±half-window for frequency smoothing  | `125` (or `0` to disable) |

**`haps2scan.freqsmooth.sh`** — pass both kb values through:
```bash
COV_SMOOTH_KB=$5
FREQ_SMOOTH_KB=$6
Rscript scripts/haps2scan.freqsmooth.R $mychr $Rfile $mydir $myoutdir $COV_SMOOTH_KB $FREQ_SMOOTH_KB
```

**`haps2scan.freqsmooth.R`** — read both, pass to code:
```r
COV_SMOOTH_KB  <- as.integer(args[5])
FREQ_SMOOTH_KB <- as.integer(args[6])
```

**`haps2scan.stable.code.R`** — compute step size from data, derive window counts
internally, never expose them:
```r
# Remove hardcoded SMOOTH_HALF <- 50
# COV_SMOOTH_KB and FREQ_SMOOTH_KB are set by the driver before this file is sourced.
if (!exists("COV_SMOOTH_KB"))  COV_SMOOTH_KB  <- 0L
if (!exists("FREQ_SMOOTH_KB")) FREQ_SMOOTH_KB <- 0L

# Derive window step size from data (kb), convert smoothing distances to window counts
positions  <- sort(unique(windows$pos))
step_kb    <- median(diff(positions)) / 1000
COV_HALF   <- if (COV_SMOOTH_KB  > 0) round(COV_SMOOTH_KB  / step_kb) else 0L
FREQ_HALF  <- if (FREQ_SMOOTH_KB > 0) round(FREQ_SMOOTH_KB / step_kb) else 0L
```

`RIDGE_FRACTION` stays as a documented constant in the code (not a user-facing parameter).

---

### Step 2 — Rename `pseudoscan` → `scan` in output filenames

Apply in both `haps2scan.stable.code.R` and `haps2scan.Apr2025.code.R`:

```r
fileout               <- paste0(mydirout, "/", myoutdir, ".scan.",          mychr, ".txt")
fileout_meansBySample <- paste0(mydirout, "/", myoutdir, ".meansBySample.", mychr, ".txt")
```

The output directory name (e.g., `ZINC2_M_freqs125/`) already encodes the smoothing
parameters — no need to repeat them in the filename.

---

### Step 3 — `haps2scan.stable.sh`

Drop it. `haps2scan.freqsmooth.sh` with `FREQ_SMOOTH_KB=0` covers that case.
Delete the file (it just becomes one less thing to maintain).

---

### Step 4 — `configs/` plotting scripts

Plotting scripts (`zinc2_125.R`, `aging_125.R`, etc.) live in `configs/` and are called
from experiment-specific submit scripts. They are experiment-specific.

Track at least one as a documented template (e.g., `configs/example_scan_plot.R`) so new
users know the expected inputs and figure format.

**Open question:** Are the current `configs/` plotting scripts already pushed, or do they
need to be added? If they document the new figure style, they should be tracked.

---

### Step 5 — Update `README.md`

Add a "Running the scan" section with both pipelines:

**Pipeline 1 — raw scan (no smoothing):**
```bash
sbatch --array=1-5 scripts/haps2scan.Apr2025.sh $Rfile $mydir $myoutdir
```

**Pipeline 2 — stabilized scan with smoothing:**
```bash
# COV_SMOOTH_KB:  covariance smoothing half-window in kb (0 = off, 125 = recommended)
# FREQ_SMOOTH_KB: frequency smoothing half-window in kb  (0 = off, 125 = recommended)
sbatch --array=1-5 scripts/haps2scan.freqsmooth.sh \
    $Rfile $mydir $myoutdir 125 125
```

Then concatenate and plot following the pattern in `scripts_oneoffs/run_all_scans_125.sh`.

---

## Implementation Order

- [x] 1. Cleanup, gitignore, track stable/freqsmooth scripts
- [x] 2. `haps2scan.freqsmooth.sh` — add `COV_SMOOTH_KB` as arg4, `FREQ_SMOOTH_KB` as arg5
- [x] 3. `haps2scan.freqsmooth.R` — read both kb args, set globals before sourcing `.code.R`
- [x] 4. `haps2scan.stable.code.R` — remove hardcoded values; compute step size from data; use `COV_HALF`/`FREQ_HALF` internally
- [x] 5. Rename `pseudoscan` → `scan` in all output filenames and job names
- [x] 6. Delete `haps2scan.stable.sh` (redundant); rename concat scripts (drop `Andreas`)
- [x] 7. Update `README.md`
- [ ] 8. Push; pull on cluster and test

---

## Open Questions

1. Are the `configs/` plotting scripts already pushed to the repo, or do they need to be added?
2. Should `concat_Chromosome_Scans.Andreas.sh` be renamed to drop `Andreas`?
