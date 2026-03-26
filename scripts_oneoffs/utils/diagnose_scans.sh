#!/bin/bash
###############################################################################
# diagnose_scans.sh
# Run from: /dfs7/adl/tdlong/fly_pool/XQTL2
# Output: diagnose_scans.out
###############################################################################

OUT="diagnose_scans.out"
SARAH="/dfs7/adl/sruckman/XQTL/XQTL2"
BASE="process/XQTL_talk"

module load R/4.2.2

exec > "${OUT}" 2>&1

sep() { echo; echo "══════════════════════════════════════════════"; echo "  $1"; echo "══════════════════════════════════════════════"; }

# ── Test / parameter files ────────────────────────────────────────────────────
sep "PUPATION: TBparameter.txt"
cat "${SARAH}/helpfiles/TBparameter.txt"

sep "ZINC2: Zinc2.test.M.txt"
cat helpfiles/ZINC2/Zinc2.test.M.txt

sep "ZINC2: Zinc2.test.F.txt"
cat helpfiles/ZINC2/Zinc2.test.F.txt

sep "AGING: AGE_SY20_M.test.txt"
cat helpfiles/AGE_SY/AGE_SY20_M.test.txt

sep "MALATHION: malathion.test.txt"
cat helpfiles/MALATHION/malathion.test.txt

# ── Hap files: inventory ──────────────────────────────────────────────────────
for EXP in ZINC2 AGING PUPATION MALATHION; do
    sep "HAP FILES: ${EXP}"
    ls -lh ${BASE}/${EXP}/R.haps.*.rds 2>/dev/null || echo "NONE FOUND"
done

# ── Hap files: pool names inside each experiment ──────────────────────────────
sep "POOL NAMES IN HAPS (chrX)"
for EXP in ZINC2 AGING PUPATION MALATHION; do
    echo "--- ${EXP} ---"
    RDS="${BASE}/${EXP}/R.haps.chrX.rds"
    [ -f "${RDS}" ] || { echo "missing"; continue; }
    Rscript -e "
h <- readRDS('${RDS}')
cat('class:', class(h), '\n')
cat('length:', length(h), '\n')
if (is.list(h)) cat('names:', paste(names(h), collapse=' '), '\n')
if (is.data.frame(h)) { cat('cols:', paste(names(h), collapse=' '), '\n'); cat('rows:', nrow(h), '\n') }
"
done

# ── Pseudoscan files: existence and Wald range ────────────────────────────────
sep "PSEUDOSCAN SUMMARIES"
declare -A SCANS
SCANS[ZINC2_M]="${BASE}/ZINC2/ZINC2_M_freqs250/ZINC2_M_freqs250.pseudoscan.txt"
SCANS[ZINC2_F]="${BASE}/ZINC2/ZINC2_F_freqs250/ZINC2_F_freqs250.pseudoscan.txt"
SCANS[AGING]="${BASE}/AGING/AGE_SY20_M_freqs250/AGE_SY20_M_freqs250.pseudoscan.txt"
SCANS[PUPATION]="${BASE}/PUPATION/PUPATION_TB_freqs250/PUPATION_TB_freqs250.pseudoscan.txt"
SCANS[MALATHION]="${BASE}/MALATHION/MALATHION_freqs250/MALATHION_freqs250.pseudoscan.txt"

for KEY in "${!SCANS[@]}"; do
    F="${SCANS[$KEY]}"
    echo "--- ${KEY} ---"
    if [ ! -f "${F}" ]; then echo "MISSING: ${F}"; continue; fi
    Rscript -e "
d <- read.table('${F}', header=TRUE)
cat('rows:', nrow(d), '\n')
cat('cols:', paste(names(d), collapse=' '), '\n')
cat('max Wald_log10p:', round(max(d\$Wald_log10p, na.rm=TRUE), 2), '\n')
cat('top 5 peaks:\n')
print(head(d[order(-d\$Wald_log10p), c('chr','pos','Wald_log10p')], 5))
"
done

# ── Pool names in hap files ───────────────────────────────────────────────────
sep "UNIQUE POOL NAMES IN HAPS (name column, chrX.rds)"
for EXP in ZINC2 AGING PUPATION MALATHION; do
    echo "--- ${EXP} ---"
    RDS="${BASE}/${EXP}/R.haps.chrX.rds"
    [ -f "${RDS}" ] || { echo "missing"; continue; }
    Rscript -e "h=readRDS('${RDS}'); cat(unique(h\$name),'\n')"
done

sep "UNIQUE POOL NAMES AND FOUNDERS IN HAPS (chrX.out.rds)"
for EXP in ZINC2 AGING PUPATION MALATHION; do
    echo "--- ${EXP} ---"
    RDS="${BASE}/${EXP}/R.haps.chrX.out.rds"
    [ -f "${RDS}" ] || { echo "missing"; continue; }
    Rscript -e "
h <- readRDS('${RDS}')
cat('cols:', paste(names(h), collapse=' '), '\n')
cat('samples:', paste(sort(unique(h\$sample)), collapse=' '), '\n')
row1 <- h[1,]
if ('Names' %in% names(h)) cat('founders (row1):', paste(unlist(row1\$Names), collapse=' '), '\n')
"
done

# ── Known-good pseudoscan summaries ──────────────────────────────────────────
sep "KNOWN-GOOD: ZINC2_M (process/ZINC2_v2/ZINC2_M)"
ls process/ZINC2_v2/ZINC2_M/*.pseudoscan.txt 2>/dev/null || echo "NOT FOUND"
for f in process/ZINC2_v2/ZINC2_M/*.pseudoscan.txt; do
    [ -f "$f" ] || continue
    echo "--- $f ---"
    Rscript -e "
d <- read.table('${f}', header=TRUE)
cat('rows:', nrow(d), '\n')
cat('max Wald_log10p:', round(max(d\$Wald_log10p, na.rm=TRUE), 2), '\n')
print(head(d[order(-d\$Wald_log10p), c('chr','pos','Wald_log10p')], 3))
"
done

sep "KNOWN-GOOD: ZINC2_F (process/ZINC2_v2/ZINC2_F)"
ls process/ZINC2_v2/ZINC2_F/*.pseudoscan.txt 2>/dev/null || echo "NOT FOUND"
for f in process/ZINC2_v2/ZINC2_F/*.pseudoscan.txt; do
    [ -f "$f" ] || continue
    echo "--- $f ---"
    Rscript -e "
d <- read.table('${f}', header=TRUE)
cat('rows:', nrow(d), '\n')
cat('max Wald_log10p:', round(max(d\$Wald_log10p, na.rm=TRUE), 2), '\n')
print(head(d[order(-d\$Wald_log10p), c('chr','pos','Wald_log10p')], 3))
"
done

sep "KNOWN-GOOD: PUPATION TopBot pseudoscan"
SARAH="/dfs7/adl/sruckman/XQTL/XQTL2"
ls ${SARAH}/process/pupal_dynamic/TopBot/*.pseudoscan.txt 2>/dev/null || echo "NOT FOUND"
for f in ${SARAH}/process/pupal_dynamic/TopBot/*.pseudoscan.txt; do
    [ -f "$f" ] || continue
    echo "--- $f ---"
    Rscript -e "
d <- read.table('${f}', header=TRUE)
cat('rows:', nrow(d), '\n')
cat('max Wald_log10p:', round(max(d\$Wald_log10p, na.rm=TRUE), 2), '\n')
print(head(d[order(-d\$Wald_log10p), c('chr','pos','Wald_log10p')], 3))
"
done

sep "KNOWN-GOOD: PUPATION hap pool names (chrX.out.rds)"
RDS="${SARAH}/process/pupal_dynamic/R.haps.chrX.out.rds"
if [ -f "${RDS}" ]; then
    Rscript -e "
h <- readRDS('${RDS}')
cat('class:', class(h), '\n')
cat('names:', paste(names(h), collapse=' '), '\n')
if ('sample' %in% names(h)) cat('unique sample:', paste(unique(h\$sample), collapse=' '), '\n')
"
else
    echo "NOT FOUND: ${RDS}"
fi

sep "OUR PUPATION hap pool names (chrX.out.rds)"
RDS="process/XQTL_talk/PUPATION/R.haps.chrX.out.rds"
if [ -f "${RDS}" ]; then
    Rscript -e "
h <- readRDS('${RDS}')
cat('class:', class(h), '\n')
cat('names:', paste(names(h), collapse=' '), '\n')
if ('sample' %in% names(h)) cat('unique sample:', paste(unique(h\$sample), collapse=' '), '\n')
"
else
    echo "NOT FOUND: ${RDS}"
fi

sep "PUPATION parfile"
cat "/dfs7/adl/sruckman/XQTL/XQTL2/helpfiles/hap_pupal_parameters.R"

sep "haps2scan.freqsmooth.sh"
cat scripts/haps2scan.freqsmooth.sh

sep "haps2scan.freqsmooth.R"
cat scripts/haps2scan.freqsmooth.R

sep "haps2scan.stable.code.R"
cat scripts/haps2scan.stable.code.R

sep "NFOUNDERS CHECK (from .out.rds Groups column)"
for EXP in ZINC2 AGING PUPATION MALATHION; do
    echo "--- ${EXP} ---"
    RDS="${BASE}/${EXP}/R.haps.chrX.out.rds"
    [ -f "${RDS}" ] || { echo "missing"; continue; }
    Rscript -e "
xx1 <- readRDS('${RDS}')
cat('Nfounders:', length(xx1\$Groups[[1]][[1]]), '\n')
cat('founder names:', paste(unlist(xx1\$Names[[1]]), collapse=' '), '\n')
cat('samples:', paste(sort(unique(xx1\$sample)), collapse=' '), '\n')
"
done

sep "ZINC2_v2 scan directory listing"
ls process/ZINC2_v2/

# ── ZINC2_v2 freqs250 comparison (same stable code as ours) ──────────────────
sep "ZINC2_v2 ZINC2_M_freqs250 (stable code, same as ours)"
F="process/ZINC2_v2/ZINC2_M_freqs250/ZINC2_M_freqs250.pseudoscan.txt"
if [ -f "${F}" ]; then
    Rscript -e "
d <- read.table('${F}', header=TRUE)
cat('rows:', nrow(d), '\n')
cat('max Wald_log10p:', round(max(d\$Wald_log10p, na.rm=TRUE), 2), '\n')
print(head(d[order(-d\$Wald_log10p), c('chr','pos','Wald_log10p')], 3))
"
else
    echo "NOT FOUND: ${F}"
fi

sep "ZINC2_v2 ZINC2_F_freqs250 (stable code, same as ours)"
F="process/ZINC2_v2/ZINC2_F_freqs250/ZINC2_F_freqs250.pseudoscan.txt"
if [ -f "${F}" ]; then
    Rscript -e "
d <- read.table('${F}', header=TRUE)
cat('rows:', nrow(d), '\n')
cat('max Wald_log10p:', round(max(d\$Wald_log10p, na.rm=TRUE), 2), '\n')
print(head(d[order(-d\$Wald_log10p), c('chr','pos','Wald_log10p')], 3))
"
else
    echo "NOT FOUND: ${F}"
fi

sep "ZINC2_v2 ZINC2_M_stable (stable code, no freq smoothing)"
F="process/ZINC2_v2/ZINC2_M_stable/ZINC2_M_stable.pseudoscan.txt"
if [ -f "${F}" ]; then
    Rscript -e "
d <- read.table('${F}', header=TRUE)
cat('rows:', nrow(d), '\n')
cat('max Wald_log10p:', round(max(d\$Wald_log10p, na.rm=TRUE), 2), '\n')
print(head(d[order(-d\$Wald_log10p), c('chr','pos','Wald_log10p')], 3))
"
else
    echo "NOT FOUND: ${F}"
fi

sep "DONE"
echo "Output written to ${OUT}"
