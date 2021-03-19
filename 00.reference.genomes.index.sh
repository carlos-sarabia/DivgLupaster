#!/bin/bash

##################################################################
## Indexing and creating a dictionary for our reference genomes ##
##################################################################

# Our reference genomes are in .fasta format and require being indexed for some programs of our pipeline to run. Also, we will create a dictionary using Picard. This step must not take longer than a couple of hours for a 3.4-Gbases genome.

###################################
## Variable and PATHs definition ##
###################################
# We will have two PATHs:

# PATH where our program executables are stored (can also be defined as variables)
PROG='/path/to/programs/folder'
# PATH with the fasta reference genome (in this case domestic dog CanFam3.1 - Lindblad-Toh et al., (2005), 
# but can be also Lycaon pictus - Campana et al., (2016)
REF='/path/to/ref/CanFam3.1'

##############################################################
## Indexing and creating dictionary of the reference genome ##
##############################################################

$PROG/bwa/bwa index $REF

$PROG/samtools/samtools faidx $REF

java -jar $PROG/picard-tools/picard.jar CreateSequenceDictionary R=/path/to/ref/CanFam3.1.fa O=/path/to/ref/CanFam3.1.dict

################
## REFERENCES ##
################

# Lindblad-Toh, K., et al. (2005). Genome sequence, comparative analysis and haplotype structure of the domestic dog. Nature, 438(7069), 803–819. https://doi.org/10.1038/nature04338
# Campana, M. G., et al. (2016). Genome sequence, population history, and pelage genetics of the endangered African wild dog (Lycaon pictus). BMC Genomics, 17(1), 1–10. https://doi.org/10.1186/s12864-016-3368-9
