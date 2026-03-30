#!/bin/bash
# Run from repo root: bash scripts_oneoffs/MALATHION/malathion.sh
Rscript scripts/plot_pseudoscan.R \
  --scan    process/XQTL_talk/MALATHION/MALATHION_freqs250/MALATHION_freqs250.pseudoscan.txt \
  --label   "Malathion selected vs control" \
  --colour  black \
  --threshold 10 \
  --genes   scripts_oneoffs/MALATHION/malathion_genes.txt \
  --out     figures/slides/malathion.png \
  --format  powerpoint \
  --height  6.0
