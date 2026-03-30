# CLAUDE.md — XQTL2

## What this repo is

This is the **public XQTL2 pipeline**. Outside users clone this to run their own
experiments. It contains only generic pipeline code and shared reference data.

## Rules — read before doing anything

1. **Never add a second git remote.** This repo has one remote: `origin` pointing
   to `https://github.com/tdlong/XQTL2.git`. Never add `dev` or any other remote.

2. **Never commit project-specific files.** The following are gitignored on purpose
   and must stay that way:
   - `scripts_oneoffs/` — project submission scripts belong in the user's project repo
   - `helpfiles/<PROJECT>/` — project configs belong in the user's project repo
   - `output/`, `logs/`, `configs/`, `data/`, `process/`, `figures/`

3. **Never commit debug output, diagnostic files, or planning documents.**

4. **If asked to do project-specific work** (run a scan, debug a project, write
   a submission script for ZINC2/AGE_SY/etc.) — that work belongs in XQTL2-dev,
   not here. Say so and stop.

5. **scripts/ contains the pipeline.** Do not duplicate pipeline functionality
   elsewhere. If a script already exists in `scripts/`, improve it there.

6. **helpfiles/ contains only shared reference data** (flymap, founder bams,
   SNP frequency tables, het_bounds) and generic templates. Nothing project-specific.

## What belongs here

- `scripts/` — pipeline scripts
- `helpfiles/flymap.r6.txt`, `helpfiles/founder.bams.txt`, `helpfiles/het_bounds.txt`
- `helpfiles/FREQ_SNPs_*.cM.txt.gz`
- `helpfiles/generic_haplotype_parameters.R`
- `README.md`, `Slurm.md`, `LICENSE`, `.gitignore`, `CLAUDE.md`

## User setup

Users who want to run their own projects should follow the "Setting up your
project repo" instructions in README.md. Their project files live in a separate
repo with a `pipeline` symlink pointing at this XQTL2 clone.
