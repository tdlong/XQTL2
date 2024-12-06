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

