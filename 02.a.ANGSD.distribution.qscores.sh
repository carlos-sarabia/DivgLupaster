#!/bin/sh

####################################
## Distribution of quality scores ##
####################################

# This script will be used to print a distribution of quality scores along our 23 genomes of wild and domestic canids following 
# Matteo Fumagalli's tutorial: https://github.com/mfumagalli/ngsTools/blob/master/TUTORIAL.md#basic-filtering-using-angsd
# Printing such a distribution is important to define whether the quality of our reads is uniformly distributed and if
# we will lose much information at setting different upper and lower filters for quality. 

###################################
## Variable and PATHs definition ##
###################################
# We have four PATHs:

# PATH to a list of the 23 autosomal genomes in .bam resulting preprocessing pipeline
BAM='/path/to/bamfolder/list'
# PATH where our program executables are stored (can also be defined as variables)
PROG='/path/to/programs/folder'
# PATH with the fasta reference genome (in this case domestic dog CanFam3.1 - Lindblad-Toh et al., (2005), 
# but can be also Lycaon pictus - Campana et al., (2016)
REF='/path/to/ref/CanFam3.1'
# PATH to output folder where genotype likelihoods will be stored:
OUT='/path/to/ANGSD/output/test.qc.distribution'

#################################################
## Generating a distribution of quality scores ##
#################################################

# In a previous step of the preprocessing pipeline, we have calculated both the mean average and standard deviation of read depths
# We can calculate an upper limit of the cumulative standard deviation. We have set up an upper limit of maxDepth 800 to see where is the main portion of the area under the curve. 

$PROG/angsd/angsd -P 16 -bam $BAMFOLD/23genomes.filelist -ref $REF -out $OUT/genoplot.qc \
	-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
	-minMapQ 20 -doQsDist 1 -doDepth 1 -doCounts 1 -maxDepth 800

# To look at the files generated, move to the $OUT folder (cd $OUT).
# To see counts of quality scores:
# less -S $OUT/genoplot.qc.qs
# To see counts of per-sample depth:
# less -S $OUT/genoplot.qc.depthSample
# wc -l $OUT/genoplot.qc.depthSample
# To see counts of global depth:
# less -S $OUT/genoplot.qc.depthGlobal

################
## REFERENCES ##
################

# Lindblad-Toh, K., et al. (2005). Genome sequence, comparative analysis and haplotype structure of the domestic dog. Nature, 438(7069), 803–819. https://doi.org/10.1038/nature04338
# Campana, M. G., et al. (2016). Genome sequence, population history, and pelage genetics of the endangered African wild dog (Lycaon pictus). BMC Genomics, 17(1), 1–10. https://doi.org/10.1186/s12864-016-3368-9







