# Generic haplotype parameters template. Copy to helpfiles/<project>/hap_params.R
# and add a names_in_bam vector (see README, "Create the haplotype parameters file").

# Founders. The B panel is B1-B7 + AB8. Every founder in the reconstruction must
# be listed here -- in a test cross that includes the tester line, not just the
# panel, or its share of each genome is redistributed across the panel founders.
founders=c("B1","B2","B3","B4","B5","B6","B7","AB8")

# Window step size in bp (5000 typical; 10000 for very large experiments)
step = 5000

# Base half-window in bp. In low-recombination regions the window expands
# proportional to max_RR / local_RR so each window spans similar recombination
# distances regardless of genomic position. This is the haplotype estimation
# window and is separate from the scan's --smooth setting.
size = 50000

# Tree height cutoff for founder distinguishability (2.5 is default)
h_cutoff=2.5
