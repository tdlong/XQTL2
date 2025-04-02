# Current XQTL pipeline

## get sequences from core...

```bash
# save the email from Yuzo as blah.txt and then extrat the files
# you will have to manually edit the resulting file a little
cat blah.txt | grep http | cut -f1 -d' ' | awk '{printf("wget %s\n",$0)}' >get_data.sh

# get_data.sh lines should look like this
wget https://hts.igb.uci.edu/tdlong24102845/xR019-L8-G2-P045-TGGCTATG-TTGTCAGC-READ1-Sequences.txt.gz
# make a directory to store your raw data
mkdir data/raw/Oct28_24
# cd to folder, put get_data.sh script in folder, then run it as a slurm job
sbatch get_data.sh
# you will have to add the 4 lines below to the top of the get_data.sh file
#  so it is a slurm script, for example something like the below would allow 
#  it to run, all my scripts are designed to run on a slurm cluster and have
#  lines like the below in them

#!/bin/bash
#SBATCH --job-name=getdata
#SBATCH -A ???        ## account to charge 
#SBATCH -p standard   ## partition/queue name
```
## create a file with read name mapping to fastq files.  And map the reads to the reference genome

```bash
# first you need indexed genomes to align to
# my script expect them to be in ref/ and have the root filename dm6.fa
# you have to be careful to map your reads to the same genome as the founders!
# you can get an indexed reference genome from here
cp /dfs7/adl/tdlong/fly_pool/newpipeline_aging/XQTL_pipeline/ref/* ref/.
# here is what is now in your directory (standard indexing)
# the files are big enough to not be on git
ls ref
dm6.dict
dm6.fa
dm6.fa.amb
dm6.fa.ann
dm6.fa.bwt
dm6.fa.fai
dm6.fa.pac
dm6.fa.sa
dm6.fa.sizes
README_ref.txt
README.txt

# now you need a file that lets your name the bam files from the raw read barcodes
# There are many ways to do this, I map barcodes to read names
# In the file below the first column is F barcode, 2nd column the R barcode,
#  and 3rd is sample name (tab delimited)
# put this file in the helperfiles directory
# note the structure of the raw file names above, depending on how your
# raw files are named my scripts may need to be modified
cat helperfiles/readname.mapping.Oct28.txt
TGGCTATG	TTGTCAGC	R3con
GTCCTAGA	TTGTCAGC	R3age
ACTTGCCA	TTGTCAGC	R5con
TCTTCGTG	TTGTCAGC	R5age
TCCACTCA	TTGTCAGC	R6con
CTGTGCTT	TTGTCAGC	R6age

# now run the alignments (fq -> bam)
# you may need to define the dir where the raw data is and where you want
#  to write bams to in the script (data/raw/Oct28_24 below)
# you will need to make the array job as big as the number of samples being aligned
# (the number of line in readname mapping file above)
mkdir data/bam/Oct28_24
NN=`wc -l helperfiles/readname.mapping.Oct28.txt | cut -f1 -d' '`
sbatch --array=1-$NN scripts/fq2bam.sh helperfiles/readname.mapping.Oct28.txt data/raw/Oct28_24 data/bam/Oct28_24 

# after it finishes (it could take overnight)
ls data/bam/Oct28_24

14G  Aug 13 23:58 data/.../R1age.bam
18G  Aug 14 02:39 data/.../R1con.bam
16G  Aug 14 00:26 data/.../R2age.bam
14G  Aug 13 23:40 data/.../R2con.bam
9.0G Oct 28 18:26 data/.../R3age.bam
17G  Oct 28 23:37 data/.../R3con.bam
14G  Aug 13 23:11 data/.../R4age.bam
18G  Aug 14 02:29 data/.../R4con.bam
5.1G Oct 28 15:58 data/.../R5age.bam
13G  Oct 28 21:31 data/.../R5con.bam
5.2G Oct 28 16:03 data/.../R6age.bam
13G  Oct 28 20:22 data/.../R6con.bam
```
It may be helpful to look at the alignment script above. We align (bwa mem), sort, add readgroups, index.  The readgroups are important. 

## Go from bams to bcf to REFALT 
The scripts in this section produce files that tabulates counts of REF and ALT alleles (for well behaved SNPs) for each SNP and sample

At this step the "helpfiles/Oct28_24.bams" provides paths to all the bam files.  This file is important as it contains paths to all your poolseq samples, plus the pre-aligned founder bams.  Since the pre-aligned bams are big, I just give paths to where I store them on hpc3, clearly this requires you are part of my "group" and can see them.  If you not doing this at UCI, then you would need to download these bam files from somewhere and they would take several TBs of space.
```bash
# This could take 20+ hours
# I will put processed stuff here
mkdir process
mkdir process/Oct28_24
# add bam files, then add founders
find data/bam/Oct28_24 -name "*.bam" -size +1G > helpfiles/Oct28_24.bams
cat helpfiles/founder.bams.txt | grep "B" >>helpfiles/Oct28_24.bams
# now generate the REFALT files
sbatch scripts/bam2bcf2REFALT.sh helpfiles/Oct28_24.bams process/Oct28_24
```

## edit haplotype.parameters.R to reflect your data
This file contains a bunch of information that lets scripts that call haplotypes and run the "GWAS" scan without intervention.  So you will probably spend some time tweaking this file.  If you do sub-analyses (on subsets of your data) your would have multiple version of this file that are passed to scripts.
```bash
##########	
#  R project specific parameters
#  haplotype.parameters.R
#  note that the number are here for when I move to 6 replicates
##########
# RGs = the readgroups add to the bams earlier 
# list of founders for the specific population of the experiment
# take from the RGs of the founder bams, and generally abbreviate the standard DSPR way
# different or custom founder populations would have different founders
# the software estimates the frequency of each founder haplotype in each pooled sample
founders=c("B1","B2","B3","B4","B5","B6","B7","AB8")

# list of samples to process
# used in REFALT2hap 
# the list of sample names to consider, they should match up with the 
# names given to the earlier fq2bam script, as these are taken from the RGs of the resulting bams 
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

```

## call the haplotypes (30-60min)
You should generally at this step just call the haplotypes for all the samples you have and deal with the GWAS separately.
```bash
# define output directory
# note the libraries needed in REFALT2haps.Andreas.R!!
sbatch scripts/REFALT2haps.Andreas.sh helpfiles/haplotype_parameters.R "process.Oct28"
```

## run the scan (15min)
```bash
# same folder as above
# note libraries needed
# note the prefix that defines a folder for output (it is created by the script)
# note the editted testing parameters file for different analyses
sbatch scripts/haps2scan.Andreas.sh helperfiles/testing.parameters.A.R "process.Oct28" "TEST_A"
sbatch scripts/haps2scan.Andreas.sh helperfiles/testing.parameters.B.R "process.Oct28" "TEST_B"
```

## concatenate and summarize (10 min)
Up until this point all analyses are done chromosome-by-chromosome for speed.  Now we concatenate and generate some summary figures
```bash
# note the path to the results for each scan above
bash scripts/concat_Chromosome_Scans.Andreas.sh "process.Oct28/TEST_A"
bash scripts/concat_Chromosome_Scans.Andreas.sh "process.Oct28/TEST_B"
```

## quick and dirty summary plots
```bash
# note the file names are hardwired, so be careful where to write them
scp 'tdlong@hpc3.rcic.uci.edu:/dfs7/adl/tdlong/.../process.Oct28/TEST_A/Age*.png' TEST_A/. 
scp 'tdlong@hpc3.rcic.uci.edu:/dfs7/adl/tdlong/.../process.Oct28/TEST_B/Age*.png' TEST_B/. 
```

