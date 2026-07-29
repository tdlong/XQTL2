# scripts/legacy — superseded pipeline steps

Nothing in this directory is part of the current pipeline. These scripts are kept
so that **older analyses can be reproduced exactly**. New projects should not use
them. They are not maintained, and they are not covered by the resource tables or
the rerun-impact guidance in the main README.

If you are starting a new experiment, go to the main [README](../../README.md).

## Legacy joint QUAL caller (superseded by the founder-catalog caller)

| Script | What it did |
|---|---|
| `bam2bcf2REFALT.sh` | Called SNPs jointly across all samples + founders, one array task per chromosome, and wrote `RefAlt.<chr>.txt` |
| `run_refalt.sh` | Submitted the above as a 5-task array and printed the job ID |

**Why it was replaced.** It ascertained SNPs from the *pooled* call under a QUAL
gate, so the sample pools influenced which sites existed, and low-frequency
founder-segregating sites were lost. The founder-catalog caller decides the SNP
set from the founders alone and then counts each sample at fixed positions — same
peaks, more usable SNPs, and far cheaper to extend with new samples. See
**"Founder-catalog caller"** in the main README.

**Reproducing an old run:**

```bash
JID_REFALT=$(bash pipeline/scripts/legacy/run_refalt.sh \
    --bamlist helpfiles/<project>/bam_list.txt \
    --dir     process/<project>)
```

Note this writes `RefAlt.<chr>.txt` into the **top level** of `--dir`, which is
the pre-`Calls/Haps/Scans` layout. Run
`bash pipeline/scripts/reorganize_project.sh process/<project>` afterwards to move
it into `Calls/` before running the haplotype step.

## Experimental tiled caller (never adopted)

A scatter/gather variant of the legacy caller: tile the genome, call one tile per
array task, then concatenate. It was under validation when the founder-catalog
caller superseded the whole approach, so it was never promoted.

| Script | Role |
|---|---|
| `make_tiles.sh` | Tile a reference `.fai` into padded windows |
| `bam2bcf2REFALT.tiled.sh` | Scatter caller — one tile per array task |
| `reassemble_refalt.sh` | Gather per-tile tables into per-chromosome ones |
| `run_refalt.tiled.sh` | Wrapper that sizes the array to the tiling and chains the gather |
| `compare_refalt.sh` | Byte-identity check between two `RefAlt` directories |

`compare_refalt.sh` checks whether two `RefAlt` directories are *byte-identical*
— the right test for a refactor that should change nothing. To compare two
genuinely different callers, use `scripts/compare_refalt_calls.R`, which reports
SNP overlap and count agreement instead.

## Legacy scan (no haplotype smoothing)

| Script | Role |
|---|---|
| `haps2scan.Apr2025.sh` | Submits the unsmoothed scan as a 5-task array |
| `haps2scan.Apr2025.R` | Driver |
| `haps2scan.Apr2025.code.R` | Implementation |

Superseded by the smoothed scan (`run_scan.sh` → `smooth_haps.R` → `hap_scan.R`),
which masks unresolvable founder frequencies, fills gaps, and smooths before
testing. Kept for reference and for reproducing pre-smoothing results.

```bash
sbatch --array=1-5 pipeline/scripts/legacy/haps2scan.Apr2025.sh \
    --rfile  helpfiles/<project>/design.txt \
    --dir    process/<project> \
    --outdir <scan_name>
bash pipeline/scripts/concat_scans.sh process/<project>/Scans/<scan_name>
```

These read `Haps/R.haps.<chr>.out.rds` and write to `Scans/<scan_name>/`, the same
stage layout as the current scan.
