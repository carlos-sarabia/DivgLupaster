#!/bin/sh

#####################################
## Genotype likelihoods with ANGSD ##
#####################################

# This script will be used to generate genotype likelihoods of our 23 genomes of wild and domestic canids following 
# Matteo Fumagalli's tutorial: https://github.com/mfumagalli/ngsTools/blob/master/TUTORIAL.md#basic-filtering-using-angsd
# These genotype likelihoods will conserve the uncertainty of called genotypes in subsequent analyses. 

###################################
## Variable and PATHs definition ##
###################################
# We have four PATHs:

# PATH to a list of the 23 autosomal genomes in .bam resulting preprocessing pipeline
BAM='/path/to/bamfolder/list'
# PATH where our program executables are stored (can also be defined as variables)
PROG='/path/to/programs/folder'
# PATH with the fasta reference genome CanFam3.1 - Lindblad-Toh et al., (2005) or Lycaon pictus - Campana et al., (2016)
REF='/path/to/ref/CanFam3.1'
# PATH to output folder where genotype likelihoods will be stored:
OUT='/path/to/ANGSD/output/1.angsd/'
# (Only to be used in analyses with neutral sites) PATH with a set of genomic regions at 10-kb distances upwards and downwards 
# from genes. 
SITES='/path/to/neutral/regions/'

#####################################
## Estimating genotype likelihoods ##
#####################################

# We estimated a lower limit of coverage per genome of 5X (minIndDepth 5). 
# Max depth of all genomes combined was an estimation of the cumulative mean coverage of all 23 genomes multiplied by two.
# For some analyses, we are interested in having the neutral regions trimmed; for some analyses not (-sites option).
# Calling genotype likelihoods at the whole autosomal genome would require 5 days. We have subdivided the task and parallelized it
# across the 38 chromosomes. The longest time to run a script was around 12 hours for chromosome 1. 

$PROG/angsd/angsd -P 16 -bam $BAMFOLD/23genomes.chr01.filelist -ref $REF/chr01.fasta -out $OUT/chr01.canfam -r chr01 \
	-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
	-minMapQ 20 -minQ 20 -minInd 23 -setMinDepth 46 -setMaxDepth 742 -doCounts 1 \
	-GL 1 -doGlf 1 -minIndDepth 5

# For those analyses that require neutral sites, we added -sites $SITES at the end of this command.

##################################
## SNP discovery (PLINK format) ##
##################################

# ANGSD can also call SNPs in PLINK format for programs like Admixture or flashPCA to work. 
# We called SNPs separately for each chromosome. Each PLINK format file was simply concatenated later with the cat command. 

$PROG/angsd/angsd -P 16 -bam $BAMFOLD/23genomes.chr01.filelist -ref $REF/chr01.fasta -out $OUT/chr01.canfam.plink -r chr01 \
	-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
	-minMapQ 20 -minQ 20 -minInd 23 -setMinDepth 46 -setMaxDepth 742 -doCounts 1 \
	-GL 1 -minIndDepth 5 -doMajorMinor 1 -doMaf 1 -skipTriallelic 1 -SNP_pval 0.00001 -doPlink 2 -doGeno -4 -doPost 1

################
## REFERENCES ##
################

# Lindblad-Toh, K., et al. (2005). Genome sequence, comparative analysis and haplotype structure of the domestic dog. Nature, 438(7069), 803–819. https://doi.org/10.1038/nature04338
# Campana, M. G., et al. (2016). Genome sequence, population history, and pelage genetics of the endangered African wild dog (Lycaon pictus). BMC Genomics, 17(1), 1–10. https://doi.org/10.1186/s12864-016-3368-9



