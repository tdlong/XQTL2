##########	
#  R project specific parameters
#  testing.parameters.R
#  note that the number are here for when I move to 6 replicates
##########

# list of samples
# used in REFALT2hap 
# the list of sample names to consider, they should match up with the 
# names given to the earlier fq2bam, as these are taken from the RGs of the bam 
names_in_bam=c("R1con","R1age","R2con","R2age","R3con","R3age","R4con","R4age","R5con","R5age","R6con","R6age")

# note the naming convention has three fields
#  Con vs Treatment -- must have two levels
#  Replicate
#  Possibly replicate within replicate, often "1"
samples=c("Con_1_1","Age_1_1","Con_2_1","Age_2_1","Con_3_1","Age_3_1","Con_4_1","Age_4_1","Con_5_1","Age_5_1","Con_6_1","Age_6_1")

# Numflies
# The number of flies in each pool
Numflies = data.frame(pool=samples,Num=c(570,1177,520,814,610,482,580,997,640,542,610,647))

# Proportion of Flies selected per replicate
ProportionSelect = data.frame(REP=c(1,2,3,4,5,6),Proportion=c(0.113,0.087,0.040,0.080,0.045,0.053))

# Mapping of Treatments to Control versus Selected
# Prefixes must be mapped to C for controls or Z for selected
TreatmentMapping = data.frame(longTRT=c("Con","Age"),TRT=c("C","Z"))


