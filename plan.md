# XQTL2 Reorganization Plan

## Current Problems

1. **Root directory clutter** — large data archives, workspace dumps, run scripts, and report
   files sitting alongside tracked source files.
2. **`scripts/` is a mix of pipeline core and experiment-specific analysis** — makes it hard
   to know what is part of the reusable pipeline vs. a one-off analysis.
3. **Multiple overlapping `haps2scan` versions** (Andreas, Apr2025, freqsmooth, stable) with
   no clear canonical one.
4. **Smoothing parameter disconnect** — a window size is passed to the submit script but the
   actual window used downstream in R is not derived from that argument in a transparent way.
5. **`.gitignore` is too minimal** — large files, workspace dumps, and outputs are not excluded.

---

## Proposed Directory Structure

```
XQTL2/
├── scripts/          # Core, reusable pipeline (TRACKED)
│   ├── scan_functions.R
│   ├── XQTL_plotting_functions.R
│   ├── fq2bam.sh
│   ├── bam2bcf2REFALT.sh
│   ├── REFALT2haps.R          (consolidate Andreas versions)
│   ├── REFALT2haps.code.R
│   ├── REFALT2haps.sh
│   ├── haps2scan.R            (single canonical version, see §4)
│   ├── haps2scan.code.R
│   ├── haps2scan.sh
│   └── concat_Chromosome_Scans.R
│
├── scripts_oneoffs/         # Experiment-specific scripts (NOT TRACKED — in .gitignore)
│   ├── aging.R, malathion.R, zinc2.R, pupation.R
│   ├── run_*.sh, submit_*.sh
│   ├── plot_single_scan.R, plot_pseudoscan.R, plot_for_slides.R, etc.
│   ├── diagnose_scans.sh, cluster_check.sh, check_param_files.sh
│   └── *.out, *_report.txt
│
├── configs/          # Parameter/design files (TRACKED or selectively tracked)
├── data/             # Data directories (NOT TRACKED)
├── ref/              # Reference genome (NOT TRACKED)
├── helpfiles/        # Helper files (selectively tracked)
├── figures/          # Output figures (NOT TRACKED)
└── process/          # Intermediate pipeline outputs (NOT TRACKED)
```

---

## Step 1 — Update `.gitignore`

Add the following to `.gitignore`:

```
# large data/archives
*.tar
*.tar.gz
*.bam
*.bai

# R workspace
.RData
.Rhistory

# pipeline outputs / reports
*.out
*_report.txt
diagnose_*.out

# directories that are all output/data
data/
ref/
figures/
process/
scripts_oneoffs/

# temp
scripts/temp/
```

---

## Step 2 — Cluster Cleanup Script (`scripts/reorganize_cluster.sh`)

A script to be run **once** on the cluster to move files into the new layout.
It will:
- Create `scripts_oneoffs/` directory
- Move experiment-specific scripts from `scripts/` → `scripts_oneoffs/`
- Move root-level run scripts and reports → `scripts_oneoffs/`
- Leave core tracked scripts untouched

**Files to move from `scripts/` → `scripts_oneoffs/`:**
```
aging.R  malathion.R  zinc2.R  pupation.R  pupation.R
run_all_scans.sh  run_all_scans_125.sh  run_pupation_scan.sh
run_malathion_scan.sh  run_aging_scan.sh
submit_freqs100_scans.sh  submit_slide_plots.sh  submit_plots.sh
plot_pseudoscan.R  plot_single_scan.R  plot_zinc_overlay.R
plot_for_slides.R  make_scan_plot.R
diagnose_scans.sh  cluster_check.sh  check_param_files.sh
extract_downsampling.R  diag_freqsmooth.R
find_refalt_files.sh  fix_refalt_b3852.sh  merge_b3852.sh
fq2bam_illumina.sh  malathion_genes.txt
```

**Files to move from root → `scripts_oneoffs/` (or remove):**
```
diagnose_scans.out  refalt_report.txt  param_report.txt  cluster_report.txt
run_freqs250_all.sh  run_chr2R_freqsmooth.sh  run_stable_all.sh
```

**Files to keep at root (tracked):**
```
README.md  LICENSE  .gitignore  plan.md
scripts/   configs/  helpfiles/  data/  ref/  figures/  process/
```

**Large files to leave in place but NOT track:**
```
founders_bam_files.tar (115G)  AGE_SY.tar (151M)  .RData  .Rhistory
```

---

## Step 3 — Consolidate `haps2scan` Versions

| File | Status | Fate |
|------|--------|------|
| `haps2scan.Andreas.*` | Tracked, original | Archive → keep for git history, superseded |
| `haps2scan.Apr2025.*` | Tracked, current best | Becomes new canonical (rename or keep) |
| `haps2scan.freqsmooth.*` | Untracked (cluster only) | **Keep in scripts/, track, then integrate smooth option** |
| `haps2scan.stable.*` | Untracked (cluster only) | **Keep in scripts/, track — supersedes Apr2025?** |

**Action needed before this step:** Share `haps2scan.freqsmooth.R`,
`haps2scan.stable.R`, and `haps2scan.stable.code.R` so their logic can be
reviewed and integrated cleanly.

Goal: one canonical set of files:
```
scripts/haps2scan.R        (wrapper, reads args, sources .code.R)
scripts/haps2scan.code.R   (doscan2 and all logic)
scripts/haps2scan.sh       (SLURM submit, passes smooth_window arg)
```

---

## Step 4 — Smoothing as a Clear Pipeline Option

### Current disconnect
`haps2scan.freqsmooth.sh` passes a window size argument, but inside the R
script the window used for smoothing may be hardcoded or not clearly tied to
that argument.

### Fix
1. `haps2scan.sh` accepts an optional `smooth_window` argument (default = 0 = no smoothing):
   ```bash
   sbatch scripts/haps2scan.sh $Rfile $mydir $myoutdir $smooth_window
   ```
2. `haps2scan.R` passes `smooth_window` to `haps2scan.code.R`:
   ```r
   smooth_window <- as.integer(args[5])   # 0 = no smoothing
   ```
3. Inside `haps2scan.code.R`, smoothing is applied only if `smooth_window > 0`:
   ```r
   if (smooth_window > 0) {
     xx1 <- smooth_haps(xx1, window = smooth_window)
   }
   ```
4. Output filenames encode the window:
   ```r
   tag <- if (smooth_window > 0) paste0(".smooth", smooth_window) else ""
   fileout <- paste0(mydirout, "/", myoutdir, tag, ".pseudoscan.", mychr, ".txt")
   ```
   This way the window is always visible in the output file name.

### `smooth_haps()` function
Will be added to `scan_functions.R`. Logic to integrate from
`haps2scan.freqsmooth.R` / `haps2scan.stable.code.R` once those are shared.

---

## Implementation Order

- [ ] **1.** Update `.gitignore`
- [ ] **2.** Write and push `scripts/reorganize_cluster.sh`
- [ ] **3.** Pull on cluster, run `reorganize_cluster.sh`, verify layout
- [ ] **4.** Share untracked smoothing scripts (`freqsmooth`, `stable`) for review
- [ ] **5.** Consolidate `haps2scan` versions into single canonical set with smooth option
- [ ] **6.** Update `README.md` to reflect new structure and smoothing option
- [ ] **7.** Final push; archive or delete superseded `haps2scan.Andreas.*` and `haps2scan.Apr2025.*`

---

## Open Questions

1. Does `haps2scan.stable.code.R` supersede `haps2scan.Apr2025.code.R`, or
   is it a parallel approach?
2. What is the smoothing operation in `haps2scan.freqsmooth.R` — sliding
   window average on hap frequencies before the scan, or on LOD scores after?
3. Should `configs/` files be tracked? (design files are experiment-specific
   but useful as examples.)
4. Should `helpfiles/` content be trimmed or is it all needed?
