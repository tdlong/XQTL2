# XQTL2 Pipeline Modernization Plan

## What Is Already Done

- [x] `.gitignore` expanded (excludes `*.tar`, `*.bam`, `.RData`, `.Rhistory`, `*.out`,
      `*_report.txt`, `data/`, `ref/`, `figures/`, `process/`, `scripts_oneoffs/`)
- [x] `scripts/reorganize_cluster.sh` — moves one-off scripts out of `scripts/`
- [x] `scripts/fix_move_to_oneoffs.sh` — corrects files accidentally moved to `analysis/`
- [x] `scripts/fix_move_freqsmooth_stable.sh` — returns freqsmooth/stable to `scripts/`
- [x] `haps2scan.stable.R`, `haps2scan.stable.code.R`, `haps2scan.stable.sh` tracked
- [x] `haps2scan.freqsmooth.R`, `haps2scan.freqsmooth.sh` tracked
- [x] Bug fix: `nrepl=1` crash in `wald.test3` / `mn.covmat` (`scan_functions.R`)

---

## What We Now Know About the Scripts

### `haps2scan.stable.code.R` supersedes `haps2scan.Apr2025.code.R`

`stable.code.R` refactors the per-window scan into three helpers and adds:

| Feature | stable.code.R | Apr2025.code.R |
|---------|--------------|----------------|
| Covariance matrix smoothing across windows | YES — ±`SMOOTH_HALF` running mean | No |
| Eigenvalue regularization (ridge) | YES — floor at `RIDGE_FRACTION * mean` | No |
| Optional frequency vector smoothing | YES — ±`FREQ_SMOOTH_HALF` | No |
| nrepl=1 handled correctly | YES — loop always uses `p1[i,]` | Fixed separately |

`Apr2025.code.R` is now superseded. `Andreas.code.R` is older still.

### `freqsmooth` and `stable` are the same code

`haps2scan.freqsmooth.R` = `haps2scan.stable.R` + reads `FREQ_SMOOTH_HALF` from `args[5]`.
Both source `haps2scan.stable.code.R`. They should be consolidated into one driver.

### The parameter disconnect (confirmed)

`SMOOTH_HALF = 50` and `RIDGE_FRACTION = 0.01` are hardcoded in `stable.code.R`.
- The `.sh` submits do **not** pass `SMOOTH_HALF` — it is never a runtime argument.
- `FREQ_SMOOTH_HALF` IS passed by `freqsmooth.R` but NOT by `stable.R`.
- Units are undefined in the argument name; "50" means ±250 kb (windows are ~5 kb apart).
- Output filenames never encode the smoothing window — two runs with different windows
  produce identically-named output files.

---

## Remaining Work

### Step 1 — Consolidate into one canonical `haps2scan` set

Replace `stable.R`, `stable.sh`, `freqsmooth.R`, `freqsmooth.sh` with a single set:

```
scripts/haps2scan.R        (driver: reads args, sets globals, sources .code.R)
scripts/haps2scan.code.R   (rename of stable.code.R, with param fixes)
scripts/haps2scan.sh       (SLURM submit, passes all args explicitly)
```

Keep `haps2scan.Andreas.*` and `haps2scan.Apr2025.*` in git history — no need to delete,
but they are no longer the canonical entry points.

---

### Step 2 — Parameters in kb, not window counts

Windows are spaced ~5 kb apart. Users should pass distances in **kb** — that is
interpretable. The R code converts internally.

**Submit interface:**
```bash
sbatch scripts/haps2scan.sh  $Rfile  $mydir  $myoutdir  $COV_SMOOTH_KB  $FREQ_SMOOTH_KB
```

| Argument | Meaning | Default | Example |
|----------|---------|---------|---------|
| `COV_SMOOTH_KB` | ±half-window for covariance smoothing | `0` (off) | `250` |
| `FREQ_SMOOTH_KB` | ±half-window for frequency smoothing | `0` (off) | `125` |

**Inside `haps2scan.R`:**
```r
COV_SMOOTH_KB  <- as.integer(args[5])   # 0 = no covariance smoothing
FREQ_SMOOTH_KB <- as.integer(args[6])   # 0 = no frequency smoothing

# Convert to window counts (windows ~5 kb apart)
SMOOTH_HALF      <- round(COV_SMOOTH_KB  / 5)
FREQ_SMOOTH_HALF <- round(FREQ_SMOOTH_KB / 5)
```

`RIDGE_FRACTION` stays as a documented constant in `haps2scan.code.R` — it is a
numerical tuning parameter, not a run-to-run variable. Document its value and meaning
with a comment.

**Preferred default for new runs:** `COV_SMOOTH_KB = 125` (±25 windows, ~±125 kb).
The old hardcoded value of 50 windows = 250 kb; the new preference is 125 kb.

---

### Step 3 — Rename output from `pseudoscan` → `scan`

Simpler, cleaner. Output filenames encode the smoothing so runs are distinguishable:

```r
cov_tag  <- if (COV_SMOOTH_KB  > 0) paste0(".cov",  COV_SMOOTH_KB)  else ""
freq_tag <- if (FREQ_SMOOTH_KB > 0) paste0(".freq", FREQ_SMOOTH_KB) else ""
tag      <- paste0(cov_tag, freq_tag)

fileout             <- paste0(mydirout, "/", myoutdir, tag, ".scan.", mychr, ".txt")
fileout_meansBySample <- paste0(mydirout, "/", myoutdir, tag, ".meansBySample.", mychr, ".txt")
```

Examples:
```
myoutdir.scan.chr2L.txt                    # no smoothing
myoutdir.cov125.scan.chr2L.txt             # covariance smoothed ±125 kb
myoutdir.cov250.freq125.scan.chr2L.txt     # both smoothed
```

---

### Step 4 — Standard plotting output

**Open question:** Which plots should be standard pipeline output (generated automatically
as part of `haps2scan`), vs. kept as one-off scripts in `scripts_oneoffs/`?

Options:
- A. Pipeline produces only the `.txt` scan tables; all plots stay in `scripts_oneoffs/`
- B. A single `plot_scan.R` script is added to `scripts/` that takes the scan output and
     produces standard manhattan plots — called optionally at the end of `haps2scan.sh`

**Action needed:** Share which plot formats from `scripts_oneoffs/plot_pseudoscan.R`
(or similar) should become standard. Once confirmed, we can decide A vs B.

---

### Step 5 — Update `README.md`

Add a section "Running the scan (current pipeline)" that shows the new `haps2scan.sh`
call with named arguments, explains the kb units, and gives examples for:
- No smoothing (fastest, baseline)
- Covariance smoothing only (recommended default)
- Both smoothing types

Mark `haps2scan.Andreas.*` and `haps2scan.Apr2025.*` as legacy in the README.

---

## Implementation Order

- [x] 1. Cleanup and gitignore
- [x] 2. Track freqsmooth and stable scripts
- [ ] 3. Consolidate freqsmooth + stable → `haps2scan.R / .code.R / .sh`
- [ ] 4. Implement kb-unit args with internal conversion
- [ ] 5. Rename `pseudoscan` → `scan` in output filenames
- [ ] 6. Resolve plotting question (Step 4 above)
- [ ] 7. Update `README.md`
- [ ] 8. Final push and cluster pull

---

## Open Questions

1. **Plotting:** Which formats from the one-off plot scripts should become standard
   pipeline output? (See Step 4)
2. **`configs/`:** Should example config/design files be tracked? Useful as templates.
3. **`helpfiles/`:** Anything in there that needs trimming or is it all still relevant?
4. **Legacy scripts:** Leave `Andreas.*` and `Apr2025.*` in `scripts/` (they just sit
   there, doing no harm), or move them to `scripts_oneoffs/` to reduce clutter?
