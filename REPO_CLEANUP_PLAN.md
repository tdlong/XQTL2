# XQTL2 Repository Cleanup Plan

## Goal

Two clean, separate repos everywhere (this machine, server, laptop):

- **XQTL2** — public pipeline. Contains only generic pipeline code and reference
  data. Any user clones this and never modifies it. Updated via `git pull`.
- **XQTL2-dev** (or MyLab-XQTL for other users) — private project repo. Contains
  all project-specific work: helpfiles, scripts_oneoffs, configs, data, process.
  References XQTL2 via a symlink called `pipeline`. Never pushes to XQTL2.

The symlink is the key architectural decision: the project repo contains a symlink
called `pipeline` pointing at the local XQTL2 clone. All project scripts use
`pipeline/scripts/run_scan.sh`, `pipeline/helpfiles/het_bounds.txt` etc.
Paths work identically on every machine — you just point the symlink at the right place.

A third repo, **XQTL2-project**, serves as a public template that users instantiate
to create their own project repo. It contains the malathion worked example so users
can verify their full setup works before analyzing their own data.

---

## Constraints

- XQTL2 has outside users — we cannot rewrite its git history (no force push)
- TBs of data on the server must not be touched
- Every step must be reversible or at least non-destructive
- The solution must work for outside users setting up their own projects too

---

## Three Repos

### XQTL2 (public, pipeline only)
What a user sees when they clone this:
- `scripts/` — all pipeline scripts
- `helpfiles/flymap.r6.txt`
- `helpfiles/founder.bams.txt`
- `helpfiles/B_founders.bams.txt`
- `helpfiles/FREQ_SNPs_Apop.cM.txt.gz`
- `helpfiles/FREQ_SNPs_Bpop.cM.txt.gz`
- `helpfiles/het_bounds.txt`
- `helpfiles/snp_tables/README.md`
- `helpfiles/generic_haplotype_parameters.R`
- `helpfiles/A_generic_haplotype_parameters.R`
- `README.md` — pipeline docs + instructions for setting up a project repo
- `Slurm.md`, `LICENSE`, `.gitignore`, `CLAUDE.md`

### XQTL2-project (public, template project repo)
Users instantiate this as their own private repo on GitHub (one click).
Contains a fully working malathion example so users can test their setup
end-to-end before running their own data. Structure:
```
XQTL2-project/
├── .gitignore              ← ignores data/, process/, ref/, figures/, output/, logs/, pipeline
├── README.md               ← setup instructions + how to run malathion test
├── helpfiles/
│   └── malathion_test/     ← barcodes, design file, hap params for malathion
└── scripts_oneoffs/
    └── malathion_test/     ← submission script calling pipeline/scripts/
```
The `pipeline` symlink is NOT committed — it is created by the user during setup.

### XQTL2-dev (private, Long lab projects)
Same structure as XQTL2-project but contains all lab projects:
```
XQTL2-dev/
├── pipeline -> ../XQTL2    ← symlink (not tracked in git)
├── helpfiles/
│   ├── ZINC2/
│   ├── AGE_SY/
│   └── ...
├── scripts_oneoffs/
│   ├── ZINC2/
│   ├── AGE_SY/
│   └── ...
├── data/                   ← gitignored
├── process/                ← gitignored
├── figures/                ← gitignored
└── CLAUDE.md
```

---

## Target Directory Layout (every machine)

```
fly_pool/                        (or ~/Desktop/Cursor/ on laptop/desktop)
├── XQTL2/                       ← clean public pipeline clone, read-only
└── XQTL2-dev/                   ← private lab project repo
    ├── pipeline -> ../XQTL2     ← symlink
    └── ...
```

---

## Note on Symlinks and SLURM

Symlinks work correctly on shared filesystems (NFS, GPFS, DFS) with SLURM.
The OS resolves symlinks transparently. No special handling needed.

---

## Step-by-Step Plan

### Phase 0 — Update scripts_oneoffs to use pipeline/ paths

Before cleaning XQTL2, audit and update all scripts in scripts_oneoffs/ that
call pipeline scripts. Change `scripts/` → `pipeline/scripts/` and
`helpfiles/het_bounds.txt` → `pipeline/helpfiles/het_bounds.txt` etc.
This way they work correctly once the symlink is in place.

---

### Phase 1 — Clean up XQTL2 on this machine and push to GitHub

**Step 1.1** — Remove project-specific files from XQTL2 in a single commit:
- All project helpfiles: `helpfiles/AGE_SY/`, `helpfiles/ZINC2/`, `helpfiles/MTX/`,
  `helpfiles/JUICE/`, `helpfiles/STARVE/`, `helpfiles/YW/`, `helpfiles/B6885/`,
  `helpfiles/MALATHION/`, `helpfiles/pupalHeight2/`, `helpfiles/AGE_Aug13_24/`,
  `helpfiles/ZINC_Hanson/`, `helpfiles/malathion_test/`
- `scripts_oneoffs/` — entire directory
- `output/chrX_diag.txt`
- `plan.md`, `PLAN_freqsmooth_pipeline.md`
- `configs/`

**Step 1.2** — Update `.gitignore` to permanently block these paths from
ever being tracked in XQTL2 again.

**Step 1.3** — Write `CLAUDE.md` for XQTL2.

**Step 1.4** — Update `README.md`:
- Remove project-specific content
- Add "Setting up your project repo" section pointing users to XQTL2-project
- Document the symlink setup
- Document how to run the malathion test from XQTL2-project

**Step 1.5** — Push to origin.

**Step 1.6** — Remove `dev` remote from XQTL2 on this machine:
```bash
git remote remove dev
```

---

### Phase 2 — Create XQTL2-project repo on GitHub

**Step 2.1** — Create new public repo `XQTL2-project` on GitHub and mark it
as a template repository (Settings → check "Template repository").

**Step 2.2** — Initialize locally with the correct structure:
```
helpfiles/malathion_test/   ← moved from XQTL2
scripts_oneoffs/malathion_test/  ← submission script using pipeline/ paths
.gitignore
README.md
```

**Step 2.3** — The malathion_test submission script calls:
```bash
bash pipeline/scripts/run_scan.sh ...
```

**Step 2.4** — Push to GitHub.

---

### Phase 3 — Set up XQTL2-dev correctly on this machine (Desktop)

**Step 3.1** — Remove `origin` remote from current XQTL2 folder and rename:
```bash
cd /Users/tdlong/Desktop/Cursor/XQTL2
git remote remove origin
cd ..
mv XQTL2 XQTL2-dev
```

**Step 3.2** — Clone fresh XQTL2:
```bash
git clone https://github.com/tdlong/XQTL2.git /Users/tdlong/Desktop/Cursor/XQTL2
```

**Step 3.3** — Create the symlink:
```bash
cd /Users/tdlong/Desktop/Cursor/XQTL2-dev
ln -s ../XQTL2 pipeline
```

**Step 3.4** — Write `CLAUDE.md` for XQTL2-dev.

**Step 3.5** — Verify:
```bash
git -C XQTL2 remote -v      # only origin -> XQTL2.git
git -C XQTL2-dev remote -v  # only dev -> XQTL2-dev.git
ls -la XQTL2-dev/pipeline   # symlink -> ../XQTL2
```

---

### Phase 4 — Fix the server

**Step 4.1** — Rename current XQTL2 to XQTL2-dev:
```bash
mv /dfs7/adl/tdlong/fly_pool/XQTL2 /dfs7/adl/tdlong/fly_pool/XQTL2-dev
```

**Step 4.2** — Remove origin remote, switch to main:
```bash
cd /dfs7/adl/tdlong/fly_pool/XQTL2-dev
git remote remove origin
git checkout main
```

**Step 4.3** — Clone fresh XQTL2 (Phase 1 must be complete first):
```bash
git clone https://github.com/tdlong/XQTL2.git /dfs7/adl/tdlong/fly_pool/XQTL2
```

**Step 4.4** — Create the symlink:
```bash
cd /dfs7/adl/tdlong/fly_pool/XQTL2-dev
ln -s ../XQTL2 pipeline
```

**Step 4.5** — Verify:
```bash
git -C /dfs7/adl/tdlong/fly_pool/XQTL2 remote -v      # only origin
git -C /dfs7/adl/tdlong/fly_pool/XQTL2-dev remote -v  # only dev
ls -la /dfs7/adl/tdlong/fly_pool/XQTL2-dev/pipeline   # symlink -> ../XQTL2
```

---

### Phase 5 — Fix the laptop

**Step 5.1** — Check current state:
```bash
git remote -v   # in whatever XQTL2 folder exists on laptop
```

**Step 5.2** — Apply same steps as Phase 3 (Desktop).

---

### Phase 6 — CLAUDE.md files

**XQTL2/CLAUDE.md:**
- This is the public pipeline repo. Outside users clone this.
- Never add a second remote.
- Never commit project-specific files.
- scripts_oneoffs/, helpfiles/<PROJECT>/, output/, configs/ are gitignored on purpose.
- If asked to do project-specific work, it belongs in XQTL2-dev not here.

**XQTL2-dev/CLAUDE.md:**
- This is the private lab project repo.
- The pipeline is at `pipeline/` (symlink to ../XQTL2).
- Read `pipeline/README.md` before writing any script.
- Never rewrite pipeline scripts — call them via `pipeline/scripts/`.
- scripts_oneoffs/<project>/ contains only submission scripts.
- Never push to XQTL2.git. Only remote is dev (XQTL2-dev.git).

---

## Order of Operations

1. **Phase 0** — update scripts_oneoffs/ paths to use `pipeline/` (on this machine)
2. **Phase 1** — clean XQTL2 and push to GitHub
3. **Phase 2** — create XQTL2-project template repo
4. **Phase 3** — split Desktop into XQTL2 + XQTL2-dev with symlink
5. **Phase 6** — write CLAUDE.md files
6. **Phase 4** — fix server (requires Phase 1 complete)
7. **Phase 5** — fix laptop

**Do not start Phase 4 until Phase 1 is complete and pushed to GitHub.**

---

## What This Does NOT Fix

- The dirty git history in XQTL2. Stays forever, files will be gone after Phase 1.
- The old newpipeline_* directories on the server. Legacy, untouched.
