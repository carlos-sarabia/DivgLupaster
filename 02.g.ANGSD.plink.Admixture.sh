#!/bin/sh

#########################
## SNP-based Admixture ##
#########################

# We used the PLINK format output of step 02.e to run a PCA with Admixture - http://dalexander.github.io/admixture/download.html

###################################
## Variable and PATHs definition ##
###################################
# PATH where our program executables are stored (can also be defined as variables)
PROG='/path/to/programs/folder'
# PATH to folder where our input (output of 02.d) is stored:
IN='/path/to/ANGSD/output/4.PLINK/'
# PATH to output folder where our output will be stored:
OUT='/path/to/ANGSD/output/6.Admixture/'

#########################
## SNP-based Admixture ##
#########################

# We are using the output of step 02.d in PLINK-format bfile. We have a set of neutral SNPs in Hardy-Weinberg equilibrium.

# First we run Admixture
for K in {05..23};do 
	$PROG/Admixture/admixture --cv $IN/autosomes.canfam.pruned.bed $K -j20 | tee $OUT/admix_results/log${K}.autosomes.canfam.pruned.txt; done

# Then, we gather the log-likelihood results to draw a best-fit K graph (in R).
grep 'CV' $OUT/admix_results/log${K}.autosomes.canfam.pruned.txt >> $OUT/admix_results/admixCVerror.autosomes.canfam.pruned.txt

################
## REFERENCES ##
################

# Alexander, D. H., Novembre, J., & Lange, K. (2009). Fast model-based estimation of ancestry in unrelated individuals. Genome Research, 19(9), 1655–1664. https://doi.org/10.1101/gr.094052.109
