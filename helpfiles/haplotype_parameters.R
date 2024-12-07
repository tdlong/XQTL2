##########	
#  R project specific parameters
#  haplotype.parameters.R
#  note that the number are here for when I move to 6 replicates
##########

# list of founders for the specific population of the experiment
# take from the RGs of the bam, and generally abbreviate the standard DSPR way
founders=c("B1","B2","B3","B4","B5","B6","B7","AB8")

# list of samples
# used in REFALT2hap 
# the list of sample names to consider, they should match up with the 
# names given to the earlier fq2bam, as these are taken from the RGs of the bam 
names_in_bam=c("R1con","R1age","R2con","R2age","R3con","R3age","R4con","R4age","R5con","R5age","R6con","R6age")

# step_size (bp)
# haplotypes will be imputed every step/1000 kb at the kb (each @ 10,20,30,etc)
step = 10000

# +/- window_size (bp)
# the windows are centered on the steps above, but of width +/- size/1000 kb
# if Numflies is generally large (>>200) and seq coverage high I think 50kb is good here
# but this is for sure a tuning parameter in that bigger windows lead to better
# haplotype inference, but poorer localization.  Smaller windows lead to better
# localization (in theory), but that is offset by poor haplotype inference.  One
# indication the window is too small is noisy neighbouring -log10p values from
# the scan.
size = 50000

# tree cutoff height to claim founders cannot be distinguished from one another
# this parameter is somewhat mysterious if I am being honest
# its units are Euclidean distance = sqrt(sum((Fi-Fj)^2) between two founders
# so over a window of 500 SNPs a distance of <2.5 implies only
# 6 SNPs being different, or 1.2% of SNPs, pretty similar haplotypes
# a cutoff of 5 implies 25 SNPs different, or 5% of SNPs
# a problem is the distance is not corrected for the number of SNPs, so if a window
# has 1000 as opposed to 500 SNPs ... then a given distance implies fewer fixed
# difference SNPs. I would tend to leave this at 2.5 or perhaps 5.
h_cutoff=2.5


