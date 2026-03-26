#!/bin/bash
###############################################################################
# Fix B6885 REFALT files: merge b3852 + b3852_r columns into SnLz
#
# The samtools merge preserved separate read group sample names, so the
# REFALT tables have REF_b3852/ALT_b3852 and REF_b3852_r/ALT_b3852_r as
# separate columns. This script sums them into REF_SnLz/ALT_SnLz and drops
# the originals.
#
# Usage: bash scripts/fix_refalt_b3852.sh process/B6885
###############################################################################

DIR=${1:?Usage: bash fix_refalt_b3852.sh <process_dir>}

for CHR in chr2L chr2R chr3L chr3R chrX; do
    FILE="${DIR}/RefAlt.${CHR}.txt"
    if [ ! -f "$FILE" ]; then
        echo "SKIP: $FILE not found"
        continue
    fi

    echo -n "Fixing $FILE ... "
    cp "$FILE" "${FILE}.bak"

    # Last 4 columns are REF_b3852, ALT_b3852, REF_b3852_r, ALT_b3852_r
    # Replace with 2 columns: REF_SnLz, ALT_SnLz (summed)
    awk 'BEGIN{OFS="\t"}
    NR==1 {
        for (i=1; i<=NF-4; i++) printf "%s\t", $i
        print "REF_SnLz\tALT_SnLz"
        next
    }
    {
        for (i=1; i<=NF-4; i++) printf "%s\t", $i
        printf "%d\t%d\n", $(NF-3)+$(NF-1), $(NF-2)+$NF
    }' "$FILE" > "${FILE}.fixed"

    mv "${FILE}.fixed" "$FILE"
    echo "done (backup: ${FILE}.bak)"
done

echo ""
echo "All REFALT files fixed. Now re-run from step 3:"
echo "  bash run_B6885.sh 3"
