#!/bin/sh

####################################################
## Admixture proportions using ANGSD and ngsAdmix ##
####################################################

# This script will be used to estimate admixture proportions of our 23 genomes of wild and domestic canids following 
# Matteo Fumagalli's tutorial: https://github.com/mfumagalli/ngsTools/blob/master/TUTORIAL.md#basic-filtering-using-angsd
# We will first call genotype likelihoods in BEAGLE format and then calculate Admixture proportions. 

###################################
## Variable and PATHs definition ##
###################################
# PATH to a list of the 23 autosomal genomes in .bam resulting preprocessing pipeline
BAM='/path/to/bamfolder/list'
# PATH where our program executables are stored (can also be defined as variables)
PROG='/path/to/programs/folder'
# PATH with the fasta reference genome (in this case domestic dog CanFam3.1 - Lindblad-Toh et al., (2005), 
# but can be also Lycaon pictus - Campana et al., (2016)
REF='/path/to/ref/CanFam3.1'
# PATH to output folder where our output will be stored:
OUT='/path/to/ANGSD/output/3.ngsAdmix/'
# PATH with a set of genomic regions at 10-kb distances upwards and downwards from genes. 
SITES='/path/to/neutral/regions/'

######################################
## Estimating Admixture proportions ##
######################################

# We estimated a lower limit of coverage per genome of 5X (minIndDepth 5). 
# Max depth of all genomes combined was an estimation of the cumulative mean coverage of all 23 genomes multiplied by two.
# For some analyses, we are interested in having the neutral regions trimmed; for some analyses not. 
# We observed that running the whole autosomal genome would require more than 5 days (real time), so instead we subdivided
# the program in 38 instances -for 38 chromosomes- and ran them in parallel. The program took max. 16 hours (for chromosome 1).

##########################
# 1st step. We estimate genotype likelihoods in BEAGLE format (repeat for every autosome).

$PROG/angsd/angsd -P 16 -bam $BAM/23genomes.chr01.filelist -ref $REF/chr01.fasta -out $OUT/chr01.canfam -r chr01 \
	-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
	-minMapQ 20 -minQ 20 -minInd 23 -setMinDepth 115 -setMaxDepth 742 -doCounts 1 \
	-GL 1 -doMajorMinor 4 -doMaf 1 -skipTriallelic 1 -doGlf 2 -sites $SITES -SNP_pval 1e-6 

#########################
# 2nd step. Run ngsAdmix in beagle data

# First we need to unzip all .beagle.gz files and concatenate them. 
# We simply take all data except for the first row (which is descriptive of the columns) and concatenate them with command cat

# Then we can run ngsAdmix. 

for i in {2..23}; do 
	$PROG/ngsTools/ngsAdmix/NGSadmix -likes $OUT/autosomes.filtered.beagle -K $i -minMaf 0.05 -seed 1 \
		-o $OUT/ngsAdmix.autosomes.K$i.run 2>$OUT/run.$i.K.log
done

# The final step is plotting. We have modified the R script provided by the developers so that we could have more color categories.

for i in {2..23}; do
	Rscript $PROG/ngsTools/ngsAdmix/Scripts/plotAdmix.R -i $OUT/ngsAdmix.autosomes.K$i.run.qopt \
		-o $OUT/ngsAdmix.autosomes.$i.run.pdf
done

################
## REFERENCES ##
################

# Lindblad-Toh, K., et al. (2005). Genome sequence, comparative analysis and haplotype structure of the domestic dog. Nature, 438(7069), 803–819. https://doi.org/10.1038/nature04338
# Campana, M. G., et al. (2016). Genome sequence, population history, and pelage genetics of the endangered African wild dog (Lycaon pictus). BMC Genomics, 17(1), 1–10. https://doi.org/10.1186/s12864-016-3368-9
# Liu, Y. H., et al. (2018). Whole-Genome sequencing of African dogs provides insights into adaptations against tropical parasites. Molecular Biology and Evolution, 35(2), 287–298. https://doi.org/10.1093/molbev/msx258



