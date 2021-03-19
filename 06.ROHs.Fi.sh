#!/bin/sh

############################################################################
## Runs of Homozygosity (ROHs) and genomewide Inbreeding coefficient (Fi) ##
############################################################################

# We are calculating heterozygosity through measure of genomewide homozygosity. This can be calculated through two ways: 
# 1. Calculating directly genomewide heterozygosity (genotype likelihood-based or SNP-based)
# 2. Calculating homozygosity through Runs of Homozygosity (ROHs - genotype likelihood-based or SNP-based)

###################################
## Variable and PATHs definition ##
###################################

# PATH where our program executables are stored (can also be defined as variables)
PROG='/path/to/programs/folder'
# PATH to the .bam files. 
IN='/path/to/bamfolder'
# PATH to the PLINK-format files (output of 02.e) 
SNP='/path/to/bamfolder/individual/SNP'
# PATH with the fasta reference genome CanFam3.1 - Lindblad-Toh et al., (2005) 
REF='/path/to/ref/CanFam3.1'
# PATH to output folder where our ngsPSMC output will be stored:
OUT='/path/to/ANGSD/output/10.Fi/'

#####################################################
## 1. Calculating directly genomewide homozygosity ## 
#####################################################

######################
## a. ngsF analysis ##
######################

# The first method requires the use of ngsF from the ngsTools software (https://github.com/mfumagalli/ngsTools). 
# A tutorial on how to use ngsF can be found here: https://github.com/mfumagalli/ngsTools/blob/master/TUTORIAL.md#inbreeding

# First of all, we call genotype likelihoods in binary format (-doGlf 3). This format is accepted by ngsF. 
# We are using two populations: afr_north (composed of Senegal, west Morocco, east Morocco, Algeria and Egypt)
# and afr_east (Ethiopia and Kenya). Maximum depth will be 2X the cumulative mean depth of every genome involved. 

$PROG/angsd/angsd -P 10 -bam $IN/afr_north.filelist -ref $REF -out $OUT/afr_north.glf3 \
	-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
	-minMapQ 20 -minQ 20 -minInd 5 -setMinDepth 5 -setMaxDepth 143 -doCounts 1 \
	-GL 1 -doGlf 3 -minIndDepth 5 -doMajorMinor 1 -doMaf 1 -skipTriallelic 1 -SNP_pval 1e-3

# With ngsF, we estimate individual genomewide inbreeding coefficients (Fi) with 20 iterations. 
# afr_north.indF will be the final output.

nsites=`zcat $OUT/afr_north.glf3.mafs.gz | tail -n+2 | wc -l` 

zcat $OUT/afr_north.glf3.glf.gz | $PROG/ngsTools/ngsF/ngsF --n_ind 5 --n_sites $nsites --glf - --out $OUT/afr_north.approx_indF --approx_EM --init_values u --n_threads 4

zcat $OUT/afr_north.glf3.glf.gz | $PROG/ngsTools/ngsF/ngsF.sh --n_ind 5 --n_sites $nsites --glf - --out $OUT/afr_north.20it.indF --init_values $OUT/afr_north.approx_indF.pars --n_threads 4

#######################
## b. PLINK analysis ##
#######################

# With PLINK, we make use of PLINK-format SNP files generated in script 02.e. We directly called SNPs in two populations, 
# afr_north and afr_east. Output will be afr_north.tfam and afr_north.tped.

$PROG/plink/plink --tfile $SNP/afr_north --het --dog --allow-no-sex --out $OUT/afr_north.het

cat afr_north.het.het | awk '{print $6}'

##############################################
## 2. Calculating homozygosity through ROHs ## 
##############################################

##########################################################
## a. Genotype-likelihood based ROH calculation (ROHan) ##
##########################################################

# We used a novel program called ROHan (Renaud et al., 2019), specifically designed to detect ROHs in low coverage data. 
# The program (and a tutorial) can be found here: https://github.com/grenaud/rohan
# ROHan requires a file with the sequencing error rates in an Illumina platform ($error) and an estimation of rates of local 
# heterozygosity (rohmu). Also, defining a desired window size (5 Mb, for example) is desirable. 
# ROHan will estimate: local rates of heterozygosity, local MCMC and HMM posterior per window according to minimum, mid and maximum 
# estimates of heterozygosity, average coverage and estimation of fraction of the genome in ROH and heterozygosity outside ROHs.

$PROG/rohan/rohan --auto $OUT/afr_wolf.algeria.rohan --err $error -t 10 --rohmu 5e-5 -o $OUT/afr_wolf.algeria.autosomes --size 500000 $ref $bam

#############################################
## b. SNP-based ROH calculation with PLINK ##
#############################################

# First, PLINK-format SNP files were generated per population as in step 02.e. 
# SNPs were pruned for linkage disequilibrium as in Sams & Boyko, (2019). 

# We generate a bfile and filter for SNPs in close linkage disequilibrium in 200-bp windows with a step size of 100 bp and R2=0.9.
$PROG/plink/plink --tfile $SNP/afr_east.autosomes --allow-no-sex --dog --make-bed --recode --out $OUT/afr_east.autosomes

$PROG/plink/plink --bfile $OUT/afr_east.autosomes --allow-no-sex --dog --make-bed --missing-genotype N --indep-pairwise 200 100 0.90 --maf 0.05 --out $OUT/afr_east.autosomes.prune

# After having pruned, we use the list of SNPs pruned to filter the original PLINK format file. 
$PROG/plink/plink --file $SNP/afr_east.autosomes --allow-no-sex --dog --extract $OUT/afr_east.autosomes.prune.prune.in --make-bed --out $OUT/afr_east.autosomes.pruned

# Finally, we estimate Runs of Homozygosity as defined in Sams & Boyko, (2019). 
$PROG/plink/plink --bfile $OUT/afr_east.autosomes.pruned --homozyg --homozyg-snp 41 --homozyg-kb 500 --homozyg-window-snp 41 --homozyg-window-missing 0 --homozyg-window-threshold 0.05 --homozyg-window-het 0 --homozyg-density 5000 --homozyg-gap 1000 --dog --out $OUT/afr_east.autosomes.pruned.samsboyko

################
## REFERENCES ##
################

# Renaud, G., et al. (2019). Joint Estimates of Heterozygosity and Runs of Homozygosity for Modern and Ancient Samples. Genetics, 212(July), 587–614. https://doi.org/10.1534/genetics.119.302057

# Sams, A. J., & Boyko, A. R. (2019). Fine-scale resolution of runs of homozygosity reveal patterns of inbreeding and substantial overlap with recessive disease genotypes in domestic dogs. G3: Genes, Genomes, Genetics, 9(1), 117–123. https://doi.org/10.1534/g3.118.200836









