#!/bin/sh

###############################################
## Site Frequency Spectrum (SFS) calculation ##
###############################################

# Since we have only few genomes with low coverage, we are not able to run some powerful analyses - like IBD-based analyses.
# So we will use thetas and neutrality tests as an indirect estimation of recent population changes.
# The Site Frequency Spectrum (SFS) is the distribution of allele frequencies of a number of loci in a population or in an 
# individual. It is very useful to calculate genomewide heterozygosity, Fst and thetas. 

# Methods can be found in https://github.com/mfumagalli/ngsTools/blob/master/TUTORIAL.md , http://popgen.dk/angsd/index.php/SFS_Estimation and the theory behind it is presented here https://pubmed.ncbi.nlm.nih.gov/22911679/

###################################
## Variable and PATHs definition ##
###################################
# We will run an unfolded SFS calculation using the Lycaon pictus (Campana et al., 2016) genome as ancestral reference. 
# The domestic dog (CanFam3.1) genome will be the reference "modern" reference. 

# PATH where our program executables are stored (can also be defined as variables)
PROG='/path/to/programs/folder'
# PATH to the .bam files. We can either run SFS per individual (useful to calculate Fst and individual heterozygosity), or per 
# population.
IND='/path/to/bamfolder/individual'
POP='/path/to/bamfolder/list'
# PATH with the fasta reference genome CanFam3.1 - Lindblad-Toh et al., (2005) 
REF='/path/to/ref/CanFam3.1'
# PATH with the fasta reference genome Lycaon pictus - Campana et al., (2016)
ANC='/path/to/ref/L.pictus'
# PATH to output folder where our SFS output will be stored:
OUT='/path/to/ANGSD/output/7.SFS/'
# PATH with a set of genomic regions at 10-kb distances upwards and downwards from genes. 
SITES='/path/to/neutral/regions/'

###############################################
## Site Frequency Spectrum (SFS) calculation ##
###############################################

# The first step consists in computing posterior probabilities of Sample Allele Frequency. Both for small populations and 
# individuals, we can run SFS throughout the whole autosomal genome. It takes less than 15 hours (for a population of 3 individuals)

# Example of a population SAF file being called. For Thetas and Fst calculations we filter out genic regions. Minimum depth is 
# defined by sites represented at least 5 times per genome (or 15 in total) and maximum depth is 2.5X the cumulative mean genome 
# depth.The list here points to whole autosomes of the three coyotes (Mexico, Midwest and California - see paper).

$PROG/angsd/angsd -P 8 -bam $IN/coyote.autosomes.bamlist -anc $ANC -ref $REF -out $OUT/coyote.autosomes.wogenic \
	-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
	-minMapQ 20 -minQ 20 -minInd 3 -setMinDepth 15 -setMaxDepth 119 -doCounts 1 \
	-GL 1 -doSaf 1 -sites $SITES

# Example of SAF file being called for an individual (Algerian African golden wolf).

$PROG/angsd/angsd -P 8 -bam $IN/afr_wolf.algeria.autosomes.list -anc $ANC -ref $REF -out $OUT/afr_wolf.algeria.autosomes.wogenic \
	-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
	-minMapQ 20 -minQ 20 -minInd 1 -setMinDepth 5 -setMaxDepth 71 -doCounts 1 \
	-GL 1 -doSaf 1 -sites $SITES

# Second step: We use the .saf files from previous step to call individual SFS with the program realSFS from the ANGSD package.
# It can be found here: http://www.popgen.dk/angsd/index.php/RealSFS
# Attention! realSFS consumes a lot of RAM memory. Depending on the size of your genome and to avoid crashes, you might want to 
# increase the memory settings first. realSFS takes a rather short time for a genome (15 minutes).

$PROG/angsd/misc/realSFS $OUT/coyote.autosomes.wogenic.saf.idx -P 16 > $OUT/coyote.autosomes.wogenic.sfs

$PROG/angsd/misc/realSFS $OUT/afr_wolf.algeria.autosomes.wogenic.saf.idx -P 16 > $OUT/afr_wolf.algeria.autosomes.wogenic.sfs

###############################
## Genomewide heterozygosity ##
###############################

# Genomewide heterozygosity can be simply calculated dividing the number of singletons by the total number of SFS values. 
# For example, if our individual SFS file is:

echo $OUT/afr_wolf.algeria.autosomes.wogenic.sfs
350688535	994011		998982946

# Then the genomewide heterozygosity is: 994011/(350688535+994011+998982946), calculated as in Gopalakrishnan et al., (2018)
# Since we have some medium and low genome coverages, we want to know if there is a bias at the calculation of SFS. 
# We downsampled the Kenyan genome (24X) to coverages pf 7X, 9X, 11.2X and 15X with samtools view -bs and repeated the analysis, 
# correcting the genomewide heterozygosity wherever needed. 

######################################
## Fst between pairs of individuals ##
######################################

# For Fst calculations, we will compute the 2DSFS between pairs of .saf files. 
# Example: calculating Fst between afr_wolf.algeria and other African golden wolves:

for name in afr_wolf.egypt afr_wolf.ethiopia afr_wolf.kenya afr_wolf.moroccopub afr_wolf.moroccounpub afr_wolf.senegal; do
	$PROG/angsd/misc/realSFS -P 16 $OUT/afr_wolf.algeria.autosomes.wogenic.saf.idx $OUT/$name.autosomes.wogenic.saf.idx > $OUT/afr_wolf.algeria.$name.2d.sfs
done

# This 2DSFS file will be used to estimate Fst in a 50kb sliding windows scan with realSFS as in:
# https://github.com/mfumagalli/ngsTools/blob/master/TUTORIAL.md 

for name in afr_wolf.egypt afr_wolf.ethiopia afr_wolf.kenya afr_wolf.moroccopub afr_wolf.moroccounpub afr_wolf.senegal 
do
	$PROG/angsd/misc/realSFS fst index -P 16 $OUT/afr_wolf.algeria.autosomes.wogenic.saf.idx $OUT/$name.autosomes.wogenic.saf.idx -sfs $OUT/afr_wolf.algeria.$name.2d.sfs -fstout $OUT/afr_wolf.algeria.$name.wogenic -whichFST 1
	$PROG/angsd/misc/realSFS fst print $OUT/afr_wolf.algeria.$name.wogenic.fst.idx | awk -v OFS='\t' '{sum3+=$3; sum4+=$4} END {print sum3/sum4}' > $OUT/afr_wolf.algeria.$name.gnome.fst
done

#################################
## Thetas and neutrality tests ##
#################################

# Thetas can be caltulated using ANGSD and realSFS as in https://github.com/mfumagalli/ngsTools/blob/master/TUTORIAL.md
# The way to do it has been described here: http://www.popgen.dk/angsd/index.php/Thetas,Tajima,Neutrality_tests

# First: we compute the allele frequency posterior probabilities and statistics (-doThetas) using SFS as prior information (-pest)
# Since we want to know population changes, we will calculate thetas for populations, not individuals. 

#1. Estimate pestPG
$PROG/angsd/angsd -P 10 -bam $IN/AFR_EAST.autosomes.bamlist -anc $ANC -ref $REF -out $OUT/AFR_EAST.autosomes.wogenic \
	-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
	-minMapQ 20 -minQ 20 -minInd 4 -setMinDepth 5 -setMaxDepth 106 -doCounts 1 \
	-GL 1 -doSaf 1 -doThetas 1 -pest $OUT/afr_east.autosomes.wogenic.sfs

#2. Print thetas
$PROG/angsd/misc/thetaStat print $OUT/AFR_EAST.autosomes.wogenic.thetas.idx > $OUT/AFR_EAST.autosomes.wogenic.thetas.print

#3. Check window-based thetas
$PROG/angsd/misc/thetaStat do_stat $OUT/AFR_EAST.autosomes.wogenic.thetas.idx
$PROG/angsd/misc/thetaStat do_stat $OUT/AFR_EAST.autosomes.wogenic.thetas.idx -win 50000 -step 50000 -outnames $OUT/AFR_EAST.autosomes.wogenic.thetas.out

# We then estimate mean and standard deviation of the genomic windows.
awk '{for(i=1;i<=NF;i++) {sum[i] += 10**$i; sumsq[i] += (10**$i)^2}} END {for (i=1;i<=NF;i++) {printf "%f %f \n", sum[i]/NR, sqrt((sumsq[i]-sum[i]^2/NR)/NR)}}' $OUT/AFR_EAST.autosomes.wogenic.thetas.print >> $OUT/AFR_EAST.autosomes.wogenic.thetas.meanstd.logscaled

# After knowing the mean and standard deviation of the number of sites per genome window, we can estimate a confidence interval of # 95% (mean +/- 2 standard deviations). In AFR_EAST example, it is:

# For thetas
cat AFR_EAST.autosomes.wogenic.thetas.out | awk '{for(i=4;i<=8;i++) {if ($14 > 30297) {theta[i]=($i/$14); sum[i]+=theta[i]; sumsq[i]+=(theta[i]^2)}}} END {for(i=4;i<=8;i++) {printf "theta/site= mean %f stdev %f \n", sum[i]/44010, sqrt((sumsq[i]-sum[i]^2/44010)/44010)}}'

# For neutrality tests
cat AFR_EAST.autosomes.wogenic.thetas.out | awk '{for(i=9;i<=13;i++) {if ($14 > 30297) {sum[i]+=$i; sumsq[i]+=($i^2)}}} END {for(i=9;i<=13;i++) {printf "neut.test_wgs= mean %f stdev %f \n", sum[i]/44010, sqrt((sumsq[i]-sum[i]^2/44010)/44010)}}'

################
## REFERENCES ##
################

# Gopalakrishnan, S., et al. (2018). Interspecific Gene Flow Shaped the Evolution of the Genus Canis. Current Biology, 28(21), 3441-3449.e5. https://doi.org/10.1016/j.cub.2018.08.041







