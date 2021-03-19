#!/bin/sh

############################################################
## PCA using genotype likelihoods with ANGSD and ngsTools ##
############################################################

# This script will be used to generate genotype likelihoods of our 23 genomes of wild and domestic canids following 
# Matteo Fumagalli's tutorial: https://github.com/mfumagalli/ngsTools/blob/master/TUTORIAL.md#basic-filtering-using-angsd
# These genotype likelihoods will conserve the uncertainty of called genotypes in subsequent analyses. They will be used to
# generate a PCA with ngsTools.

###################################
## Variable and PATHs definition ##
###################################
# We have several PATHs:

# PATH to a list of the 23 autosomal genomes in .bam resulting preprocessing pipeline
BAM='/path/to/bamfolder/list'
# PATH where our program executables are stored (can also be defined as variables)
PROG='/path/to/programs/folder'
# PATH with the fasta reference genome (in this case domestic dog CanFam3.1 - Lindblad-Toh et al., (2005), 
# but can be also Lycaon pictus - Campana et al., (2016)
REF='/path/to/ref/CanFam3.1'
# PATH to output folder where our output will be stored:
OUT='/path/to/ANGSD/output/2.PCA/'
# PATH with a set of genomic regions at 10-kb distances upwards and downwards from genes. 
SITES='/path/to/neutral/regions/'

#############################################
## Estimating genotype likelihoods for PCA ##
#############################################

# We estimated a lower limit of coverage per genome of 5X (minIndDepth 5). 
# Max depth of all genomes combined was an estimation of the cumulative mean coverage of all 23 genomes multiplied by two.
# For some analyses, we are interested in having the neutral regions trimmed; for some analyses not. 
# We observed that running the whole autosomal genome would require more than 5 days (real time), so instead we subdivided
# the program in 38 instances -for 38 chromosomes- and ran them in parallel. The program took max. 9 hours (for chromosome 1).

##########################
# 1st step. We estimate genotype likelihoods (repeat for every autosome).

$PROG/angsd/angsd -P 16 -bam $BAM/23genomes.chr01.filelist -ref $REF -out $OUT/chr01.canfam -r chr01 \
	-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
	-minMapQ 20 -minQ 20 -minInd 23 -setMinDepth 46 -setMaxDepth 742 -doCounts 1 \
	-GL 1 -doMajorMinor 1 -doMaf 1 -skipTriallelic 1 -minIndDepth 5 -sites $SITES -SNP_pval 1e-5 -doGeno 32 -doPost 1

#########################
# 2nd step. Run ngsCovar to perform a PCA

# First we need to unzip all .geno.gz and .mafs.gz files and concatenate them (simply with command cat)
# Then we need to estimate the number of sites: 
NSITES = `cat $OUT/autosomes.canfam.mafs | tail -n+2 | wc -l`

# Then we can run ngsCovar. As in Fumagalli's tutorial, we do not want to perform genotype calling (-call 0). 
# We also do not want to normalise by allele frequency (-norm 0)
$PROG/ngsTools/ngsPopGen/ngsCovar -probfile $OUT/autosomes.canfam.geno -outfile $OUT/autosomes.canfam.covar \
	-nind 23 -nsites $NSITES -call 0 -norm 0 

########################
# 3rd. Plot the PCR. 

# We first create a plink cluster-like file defining the labelling for each sample.
Rscript -e 'write.table(cbind(seq(1,23), rep(1,23), c(rep("afr_wolf",7), rep("gr_wolf",6), rep("dog",7), rep("eth_wolf",1), rep("eu_gjackal",2))), row.names=F, sep="\t", col.names=c("FID","IID","CLUSTER"), file="autosomes.canfam.clst", quote=F)'

# Then, we run and plot the PCA. 
Rscript /home/carlos/Desktop/programs/ngsTools/Scripts/plotPCA.R -i $OUT/autosomes.canfam.covar -c 1-2 -a $OUT/autosomes.canfam.clst -o $OUT/autosomes.canfam.pca.pdf

################
## REFERENCES ##
################

# Lindblad-Toh, K., et al. (2005). Genome sequence, comparative analysis and haplotype structure of the domestic dog. Nature, 438(7069), 803–819. https://doi.org/10.1038/nature04338
# Campana, M. G., et al. (2016). Genome sequence, population history, and pelage genetics of the endangered African wild dog (Lycaon pictus). BMC Genomics, 17(1), 1–10. https://doi.org/10.1186/s12864-016-3368-9
# Liu, Y. H., et al. (2018). Whole-Genome sequencing of African dogs provides insights into adaptations against tropical parasites. Molecular Biology and Evolution, 35(2), 287–298. https://doi.org/10.1093/molbev/msx258










