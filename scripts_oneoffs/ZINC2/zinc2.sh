#!/bin/bash
# Run from repo root: bash scripts_oneoffs/ZINC2/zinc2.sh
Rscript scripts/plot_pseudoscan.R \
  --scan    process/XQTL_talk/ZINC2/ZINC2_M_freqs250/ZINC2_M_freqs250.pseudoscan.txt \
  --scan    process/XQTL_talk/ZINC2/ZINC2_F_freqs250/ZINC2_F_freqs250.pseudoscan.txt \
  --label   Males \
  --label   Females \
  --colour  "#1F78B4" \
  --colour  "#E31A1C" \
  --threshold 10 \
  --out     figures/slides/zinc2.png \
  --format  powerpoint \
  --height  5.5
