# Founder-Sites RefAlt Pipeline — Validation Plan

## Current Status — 2026-03-26

ZINC2 validation run submitted on HPC3. All 4 phases queued with dependencies.

| Phase | Job ID | Status |
|-------|--------|--------|
| 1 — Founder site catalog | 50234868 (array 1-5) | submitted |
| 2 — Per-sample RefAlt (60 samples) | see logs/ZINC2_fromsites/jobs.log | submitted |
| 3 — Consolidation | 50234934 (array 1-5) | waiting on Phase 2 |
| 4 — Comparison | 50234935 | waiting on Phase 3 |

**When returning:** check results with:
```bash
cat logs/ZINC2_fromsites/compare.50234935.out
bash scripts_oneoffs/profile_jobs.sh logs/ZINC2_fromsites/jobs.log
```

Key question: are counts exactly identical at shared sites, and how many sites drop out of the founder-only catalog vs the original joint-called RefAlt?

---

## Motivation

The current BAM → RefAlt step (bam2bcf2REFALT.sh) runs joint variant calling across
all samples simultaneously. This requires a full rerun whenever new samples are added
to a project, taking over a day. Because all samples are derived from exactly 8 inbred
B-population founders, the complete variant catalog is fixed by the founders. This
experiment tests whether we can:

1. Call variants once from founders only → permanent site catalog
2. For any new sample: genotype at catalog sites only (hours, not days)
3. Add new samples by running per-sample jobs and consolidating — no reprocessing of existing samples

## Validation Strategy

Run the new pipeline on ZINC2 (which already has original RefAlt files in process/ZINC2/)
and compare results. At shared sites, counts must be **exactly identical** — not just
correlated. The only expected differences are in which sites are included.

---

## Scripts

| Script | Purpose | Array | Partition | CPUs | Mem/CPU | Time | Module |
|--------|---------|-------|-----------|------|---------|------|--------|
| bam2founder_sites.sh | Build B-founder site catalog (one-time) | 1–5 (chr) | standard | 2 | 6G | 24:00:00 | bcftools/1.21 |
| bam2bcf2REFALT_fromsites.sh | Per-sample RefAlt at founder sites | 1–5 (chr) | standard | 1 | 4G | 12:00:00 | bcftools/1.21 |
| consolidate_refalt_fromsites.sh | Merge per-sample files → wide RefAlt | 1–5 (chr) | standard | 1 | 5G | 1:00:00 | R/4.2.2 |
| compare_refalt_zinc2.sh | Validate new vs original RefAlt | none | standard | 1 | 5G | 1:00:00 | R/4.2.2 |

Consolidation and comparison read 89–155MB files; standard/5G is well above what is
needed. See `profile_jobs.sh` output after the first run to calibrate.

---

## Run Order

```bash
# ── Phase 1: Build founder site catalog (one-time for B-population) ──────────
mkdir -p process/B_founder_sites
jid1=$(sbatch scripts_oneoffs/bam2founder_sites.sh \
         helpfiles/B_founders.bams.txt \
         process/B_founder_sites \
       | awk '{print $4}')
echo "Phase 1 job: $jid1"

# ── Phase 2: Per-sample RefAlt for ZINC2 experimental pools ──────────────────
mkdir -p process/ZINC2_fromsites
mkdir -p logs/ZINC2_fromsites
jid2_list=()
while read bam; do
  [[ $bam == /dfs7* ]] && continue   # skip founder lines
  jid=$(sbatch --dependency=afterok:$jid1 \
               --output=logs/ZINC2_fromsites/$(basename $bam .bam).%A_%a.out \
               scripts_oneoffs/bam2bcf2REFALT_fromsites.sh \
               $bam \
               process/B_founder_sites \
               process/ZINC2_fromsites \
         | awk '{print $4}')
  jid2_list+=($jid)
  echo -e "$(basename $bam .bam)_phase2\t$jid" >> logs/ZINC2_fromsites/jobs.log
done < helpfiles/ZINC2/ZINC2.bams

# Colon-separated list for dependency
jid2_dep=$(IFS=:; echo "${jid2_list[*]}")

# ── Phase 3: Consolidate per-sample files into wide RefAlt ───────────────────
jid3=$(sbatch --dependency=afterok:$jid2_dep \
              --output=logs/ZINC2_fromsites/consolidate.%A_%a.out \
              scripts_oneoffs/consolidate_refalt_fromsites.sh \
              process/ZINC2_fromsites \
              process/B_founder_sites \
        | awk '{print $4}')
echo -e "consolidate\t$jid3" >> logs/ZINC2_fromsites/jobs.log
echo "Phase 3 job: $jid3"

# ── Phase 4: Compare against original ZINC2 RefAlt ───────────────────────────
jid4=$(sbatch --dependency=afterok:$jid3 \
              --output=logs/ZINC2_fromsites/compare.%A.out \
              scripts_oneoffs/compare_refalt_zinc2.sh \
              process/ZINC2 \
              process/ZINC2_fromsites \
        | awk '{print $4}')
echo -e "compare\t$jid4" >> logs/ZINC2_fromsites/jobs.log
echo "Phase 4 job: $jid4"
```

---

## Resource Profiling

After all jobs complete, profile with:

```bash
bash scripts_oneoffs/profile_jobs.sh logs/ZINC2_fromsites/jobs.log
```

The job log is built automatically by the run commands above. For array jobs, seff
reports aggregate stats across all tasks. To profile individual array tasks, append
`_<task_id>` to the job ID (e.g. `12345_3` for task 3).

---

## What to Look For in the Comparison

| Metric | Expected | Red flag |
|--------|----------|----------|
| Counts at shared sites | Exactly identical | Any difference |
| Sites retained (%) | >95% | <90% suggests QUAL threshold needs review |
| Sites gained | 0 | Any gained sites are unexpected |
| Dropped-site founder alt freq | Low (marginal in founders too) | High founder freq + dropped = threshold issue |
| Dropped-site pool alt freq | Low (pools were pushing QUAL over 59) | High pool freq = real signal being lost |

---

## Future Use (If Validation Passes)

For any new B-population project:

```bash
# Run Phase 2 only for new samples (Phase 1 output is permanent)
sbatch scripts_oneoffs/bam2bcf2REFALT_fromsites.sh <new.bam> process/B_founder_sites <project_fromsites_dir>

# Re-consolidate (cheap — just re-joins per-sample files including new ones)
sbatch scripts_oneoffs/consolidate_refalt_fromsites.sh <project_fromsites_dir> process/B_founder_sites
```

Existing sample files in the project directory are untouched.
