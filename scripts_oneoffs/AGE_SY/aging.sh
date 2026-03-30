#!/bin/bash
# Run from repo root: bash scripts_oneoffs/AGE_SY/aging.sh
Rscript scripts/plot_pseudoscan.R \
  --scan    process/XQTL_talk/AGING/AGE_SY20_M_freqs250/AGE_SY20_M_freqs250.pseudoscan.txt \
  --label   "Males day 20 selected vs control" \
  --colour  black \
  --threshold 10 \
  --out     figures/slides/aging.png \
  --format  powerpoint \
  --height  5.8
