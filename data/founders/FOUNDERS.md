# Founder BAM provenance

The founder BAMs in this directory are what the XQTL2 pipeline uses to define and
call SNPs (via `helpfiles/A_founders.bams.txt` and `helpfiles/B_founders.bams.txt`).
This file records where they come from and how they were built, so the callset's
foundation is documented rather than an opaque tarball.

To verify any statement here yourself:
```bash
samtools view -H <founder>.bam | grep -E '^@RG|^@PG'
```
`@RG` gives the sample / flowcell / library; the `@PG` chain gives the exact
bwa / samtools / picard commands and versions used to build the file.

## Summary

- **What:** pre-aligned BAMs for the DSPR founders (A-population A1–A7,AB8;
  B-population AB8,B1–B7).
- **Source:** **resequenced PE150** libraries of **archived founder DNA** — the
  `arc.<founder>` in the BAM `@PG` paths stands for *archived*. These are the
  RESEQUENCED founders, **not** the original early-days Illumina BAMs.
- **Reference:** dm6 (FlyBase release 6 / UCSC `dm6`).
- **Packaging:** distributed as `founders_bam_files.tar` (README §2); current
  build labelled July 2025 R6 (`founders.July2025.R6.tar`, `July_2025_list.txt`).
- **Dates:** main resequenced set — *TBD* (no `DT` tag in the headers); the B5
  fix is December 2024 (see below).

## Build recipe — all founders except B5

Reconstructed from the `@PG` chain (example: `A1.dedup.bam`):

1. `bwa mem 0.7.17-r1188 -M` vs `ref/dm6.fa`, per library (samtools 1.10)
2. `samtools merge` the libraries for each founder
3. re-sort (samtools 1.15.1)
4. add read groups (`ID` = `SM` = founder name) → `<name>.RG.bam` intermediate
5. Picard `MarkDuplicates REMOVE_DUPLICATES=true` → **`<name>.dedup.bam`** (shipped)

**Read groups:** `ID` and `SM` = the founder name (`A1`, `B3`, `AB8`, …). This
`SM` tag is how the pipeline identifies founders — `catalog_build.sh` resolves
the founder set from `hap_params.R` names by matching each BAM's `SM` tag.

Example — `A1.dedup.bam`: Illumina flowcell `D109LACXX`, run `nR331` lane L3, two
barcoded libraries (`P14`, `P46`) merged.

## name → file (what the pipeline points at)

| population | founders | file |
|-----------|----------|------|
| A (`A_founders.bams.txt`) | A1–A7, AB8 | `<name>.dedup.bam` |
| B (`B_founders.bams.txt`) | AB8, B1–B4, B6, B7 | `<name>.dedup.bam` |
| B | **B5** | **`B5.RG.bam`** (composite inversion-fix file, below) |

## B5 — the exception (December 2024)

B5's archived DNA is **heterozygous for a chr2L inversion**, so a naïve
resequence mixes two haplotypes on 2L. `B5.RG.bam` is a composite built to keep
only the **standard (non-inverted) 2L haplotype**. Procedure (from
`/dfs7/adl/tdlong/DSPR/DNAseq/scripts/create_new_founders.sh` and the working
notes in `../NOT_founders/B5_stuff.txt`):

1. Standard reseq of `arc.B5` vs dm6 → `tempB5.dedup.bam`.
2. Build a B5-standard custom genome `B5.fa`:
   `bcftools consensus -s B5Xiso1 -H A` substitutes B5's **ALT** alleles (from the
   B5×iso-1 cross, `B5Xiso1.dedup.bam`) into dm6.
3. Re-align `arc.B5` to `B5.fa` → `B5_custom.bam`. Reads from the standard
   haplotype match the alt alleles **perfectly** (`NM:i:0`); inverted-haplotype
   reads do not.
4. Composite `B5.bam`:
   - **chr2L:** only `-q30` `NM:i:0` reads from `B5_custom.bam` (standard
     haplotype; inverted-haplotype reads discarded)
   - **all other chromosomes:** from `tempB5.dedup.bam` (normal reseq)
5. read-group stamped → **`B5.RG.bam`** (shipped).

Tool versions for the B5 fix differ from the main pipeline: bwa 0.7.17,
samtools 1.10, bcftools 1.10.2, picard-tools 1.87, java 17.

**Consequence — B5's chr2L is shallow *and* circular.** Coverage from the notes:
chr2L ~9.5× vs chr2R ~22×, because the inverted-haplotype reads were dropped by
design. Worse, the surviving chr2L reads were *selected* to match B5's own alleles
(`NM:i:0` vs `B5.fa`), so B5's chr2L "fixation" is guaranteed by construction, not
independent evidence. B5 is a normal founder on every other chromosome. This is why
`catalog_build.sh` drops B5 from the SNP rules **on chr2L only** by default
(`--exempt-founders B5:chr2L`) — its counts are still written to `RefAlt`
everywhere. See XQTL2 issue #15 and the catalog rules in the README.

## `../NOT_founders/` — B5 fix intermediates (not used by the pipeline)

| file | what it is |
|------|-----------|
| `B5.fa` (+ bwa index) | B5-standard custom genome (step 2) |
| `B5_custom.bam` | `arc.B5` aligned to `B5.fa` (step 3) |
| `B5Xiso1.dedup.bam` | the B5×iso-1 cross used to build `B5.fa` |
| `tempB5.dedup.bam` | standard reseq of `arc.B5` (step 1) |
| `B5_stuff.txt` | the working notes / commands |

## Unverified / TODO

- Exact sequencing date(s) of the main resequenced set (no `DT` tag in headers).
- `b3852.merged.bam` and `B5.bam` are present in this folder but are **not**
  referenced by `A_/B_founders.bams.txt`; status unconfirmed.
