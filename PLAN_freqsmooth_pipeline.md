# XQTL2 Freqsmooth Pipeline — Living Plan

**Last updated:** 2026-03-22

---

## Vision

The existing `haps2scan.freqsmooth` pipeline applies smoothing inconsistently:
covariances are always smoothed, haplotype frequencies are optionally smoothed,
heritability inputs are never smoothed. The goal here is to build and validate a
cleaner "smooth first" replacement:

1. **`smooth_haps.R`** — smooth haplotype frequencies AND covariance matrices in one
   pass before any statistics are computed. All downstream steps then operate on
   the same smoothed data.
2. **`freqsmooth_scan.R`** — Wald test + Falconer/Cutler H² on smoothed frequencies.
   Produces the same `.scan.txt` and `.meansBySample.txt` as the existing pipeline.
3. **`snp_scan.R`** — NEW. Imputes per-SNP ALT frequencies from smoothed haplotype
   estimates via dot product, then runs a df=1 Wald test at every SNP. Produces
   `.snp_scan.txt`. This is an entirely new output that did not exist before.

**If the malathion validation passes**, these three scripts replace `haps2scan.freqsmooth`
as the canonical Step 5 in the pipeline. `scripts_freqsmooth/` moves to `scripts/`.
The SNP scan becomes a standard new Step 5c in the README. The legacy
`haps2scan.Apr2025` (raw scan, no smoothing) remains available as an alternative.

---

## Infrastructure Status — All Complete

All scripts are written, optimized, profiled, and pushed to git.

### New scripts in `scripts_freqsmooth/`
- [x] `smooth_haps.R` + `smooth_haps.sh` — edge-aware running mean, ±SMOOTH_KB
- [x] `freqsmooth_scan.R` + `freqsmooth_scan.sh` — vectorized (pre-built arrays, lapply); seff: 12:32 wall / 346 MB / 2 cores; mem set to 1G/core
- [x] `snp_scan.R` + `snp_scan.sh` — fully vectorized (no per-SNP loop); seff: 5:32 wall / 726 MB / 1 core; mem set to 3G
- [x] `concat_snp_scans.sh` — merges per-chromosome snp_scan files

### New plotting engines in `scripts/`
- [x] `plot_freqsmooth_H2.R` — 5-panel H² line plot (single estimator)
- [x] `plot_freqsmooth_snp.R` — 5-panel SNP scan dot plot
- [x] `plot_H2_overlay.R` — 5-panel overlay of Falconer H² + Cutler H² in one figure

### Malathion test drivers in `helpfiles/malathion_test/`
- [x] `malathion_test_v2.sh` — full pipeline from scratch (smooth → scan → snp_scan → concat → figures)
- [x] `malathion_test_v2_rescan.sh` — rescan from smooth onward (smooth already exists)
- [x] `malathion_rerun_snp.sh` — standalone SNP scan chain (snp_scan → concat → figure → seff)
- [x] `malathion_test_v2_smooth250.sh` — 250kb smooth pipeline
- [x] `malathion_rerun_snp_250.sh` — 250kb SNP scan chain

### Malathion figure scripts
- [x] `MALATHION_TEST_v2_smooth125.R` — Wald + Falconer H² + Cutler H² (125kb hap scan)
- [x] `MALATHION_TEST_v2_smooth125_snp.R` — Wald only (125kb SNP scan)
- [x] `MALATHION_TEST_v2_smooth250.R` — Wald + combined H² overlay + meansBySample QC (250kb hap scan)
- [x] `MALATHION_TEST_v2_smooth250_snp.R` — Wald only (250kb SNP scan)
- [x] `plot_meansBySample_QC.R` — local QC script: one founder × one chromosome

---

## Current State — Malathion Test

### 125kb results (downloaded, evaluated)
- [x] Hap scan figures: `MALATHION_TEST_v2_smooth125.hap/`
- [x] SNP scan figures: `MALATHION_TEST_v2_smooth125.snp/`
- Wald peaks confirmed, SNP scan looks good, H² estimates present but noisier than expected

### 250kb results (in progress)
- Hap scan jobs: concat/figure/seff in SLURM queue (50090743–50090745; time limits updated to 2h)
- SNP scan: **not yet submitted** — run after hap jobs complete:
  ```bash
  bash helpfiles/malathion_test/malathion_rerun_snp_250.sh
  ```
- SCP when ready:
  ```bash
  scp tdlong@hpc3.rcic.uci.edu:/dfs7/adl/tdlong/fly_pool/XQTL2/process/malathion_test/MALATHION_TEST_v2_smooth250/MALATHION_TEST_v2_smooth250.hap.tar.gz .
  scp tdlong@hpc3.rcic.uci.edu:/dfs7/adl/tdlong/fly_pool/XQTL2/process/malathion_test/MALATHION_TEST_v2_smooth250/MALATHION_TEST_v2_smooth250.snp.tar.gz .
  ```

---

## Validation Checklist

Goal: confirm new pipeline gives scientifically sensible results on malathion data
before promoting it to production.

- [x] Manhattan plots show expected major peaks (malathion resistance locus prominent) — confirmed on 125kb results
- [x] Wald -log10p values comparable in magnitude to existing pipeline — confirmed
- [x] SNP scan peaks overlap haplotype scan peaks at malathion locus — confirmed
- [x] Runtime acceptable — freqsmooth_scan 12:32 / 346 MB; snp_scan 5:32 / 726 MB
- [ ] **Smoothing window comparison**: do 250kb results look cleaner than 125kb? (meansBySample QC plot is key — check all 8 founders)
- [ ] **meansBySample smoothness**: no jagged rep-to-rep oscillation at chosen window size
- [ ] **Heritability behavior**: understand and document the Falconer vs Cutler divergence (see open questions below)

---

## Open Questions

### 1. Optimal smoothing window
125kb shows some oscillation in founder frequencies at shorter scales than expected.
250kb test underway. Compare `means_qc_chr3L.png` in both tarballs side by side.
Decision needed: which window to recommend as default?

### 2. Falconer vs Cutler H² divergence
Cutler is systematically lower than Falconer, especially at strong peaks. This is
by design: in `scan_functions.R` line 181, Cutler weights each founder's contribution
by `C` (control frequency), so rare founders with large responses contribute little.
Falconer uses `C` in the denominator, upweighting exactly those founders. The
Penetrance clamping at `[Proportion/2, 2*Proportion]` additionally limits Cutler
at strong peaks.

**This is not a bug — it is the intended behavior of the Cutler estimator.** But it
means the two estimators are measuring somewhat different things. Decision needed:
- Keep both in the output and overlay plot?
- Drop one?
- Document the difference for users?

### 3. H² in SNP scan
Decided to drop H² from the SNP scan figures (Wald only). The biallelic H² at
individual SNPs is not considered meaningful at this stage.

---

## Integration Plan (approved, ready to implement)

### Guiding principles
- README is the user-facing truth. Every script change must be reflected there.
- One start-to-finish pipeline template in the README — a user copies it, fills in variables, runs it.
- Simplify: the malathion test helper scripts are test artifacts; the production pattern is one `scripts_oneoffs/<project>_pipeline.sh` per experiment.
- Keep `haps2scan.Apr2025` as a documented alternative (genuinely different — no smoothing).

### Scripts — what changes

| Action | Scripts |
|--------|---------|
| **Move** `scripts_freqsmooth/` → `scripts/` | `smooth_haps.R/.sh`, `freqsmooth_scan.R/.sh`, `snp_scan.R/.sh`, `concat_snp_scans.sh`, `prep_snp_table.R` |
| **Retire** (delete) | `haps2scan.freqsmooth.sh`, `haps2scan.freqsmooth.R`, `haps2scan.freqsmooth.code.R` — fully replaced |
| **Keep unchanged** | `haps2scan.Apr2025.*`, `concat_Chromosome_Scans.*`, `scan_functions.R`, `XQTL_plotting_functions.R`, `fq2bam.sh`, `bam2bcf2REFALT.sh`, `REFALT2haps.*`, `show_project_layout.sh` |
| **Keep, not featured** | `plot_freqsmooth_H2.R` (single-estimator H² — available but not in main docs) |
| **Delete** | `scripts_freqsmooth/` directory (now empty) |

### README — new structure

**Step 5 — Run the scan** (replaces current Step 5)

- **5a** `smooth_haps.sh` — smooth haplotype frequencies and covariances (±smooth-kb)
- **5b** `freqsmooth_scan.sh` — Wald + H² on smoothed data; depends on 5a
- **5c** `snp_scan.sh` — SNP scan; depends on 5a, runs in parallel with 5b

Legacy option remains documented:
- `haps2scan.Apr2025.sh` — raw scan, no smoothing

**Step 6 — Concatenate chromosomes**

Two concat steps — one per output type:
- `concat_Chromosome_Scans.sh` — hap scan (unchanged)
- `concat_snp_scans.sh` — SNP scan (new)

**Step 7 — Generate publication figures**

Three plot engines now documented with examples:
- `plot_pseudoscan.R` — Wald -log10(p) from hap scan
- `plot_H2_overlay.R` — Falconer + Cutler H² overlaid on one figure
- `plot_freqsmooth_snp.R` — SNP Wald -log10(p)

**"Putting it all together" template** — updated to cover Steps 5a/5b/5c + both concats + figures.
The malathion test scripts serve as the worked example.

### SNP table documentation — new subfolder

`helpfiles/snp_tables/README.md` documents:
- What `FREQ_SNPs.cM.txt.gz` is and where it comes from (founder sequencing)
- How to run `prep_snp_table.R` to extract population-specific tables (A pop, B pop)
- That this is a **one-time step per population**, not per experiment
- Output files: `FREQ_SNPs_Apop.cM.txt.gz`, `FREQ_SNPs_Bpop.cM.txt.gz`

Main README Step 5c just says: "See `helpfiles/snp_tables/` for how SNP tables are prepared."

### Implementation order

1. Create `helpfiles/snp_tables/README.md`
2. Move scripts from `scripts_freqsmooth/` to `scripts/`
3. Delete `haps2scan.freqsmooth.*` and empty `scripts_freqsmooth/`
4. Rewrite `README.md` Steps 5–7 and the "Putting it all together" template
5. Commit and push; verify `git pull` on cluster works cleanly

---

## Cluster / Environment Notes

- Cluster: `hpc3.rcic.uci.edu`, user: `tdlong`
- Project root on cluster: `/dfs7/adl/tdlong/fly_pool/XQTL2`
- Standard partition: max 6G/core; use `--mem-per-cpu` not `--mem`
- SLURM dependency: always `afterok:JOBID` — never use array task ID suffixes
- Local project root: `/Users/anthonylong/Desktop/Scripts_Code/Cursor_projects/XQTL2` (git clone)
