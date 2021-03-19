#!/bin/sh

###################
## SNP-based PCA ##
###################

# We used the PLINK format output of step 02.e to run a PCA with flashPCA.

###################################
## Variable and PATHs definition ##
###################################
# PATH where our program executables are stored (can also be defined as variables)
PROG='/path/to/programs/folder'
# PATH to folder where our input (output of 02.d) is stored:
IN='/path/to/ANGSD/output/4.PLINK/'
# PATH to output folder where our output will be stored:
OUT='/path/to/ANGSD/output/5.flashPCA/'

###################
## SNP-based PCA ##
###################

# We are using the output of step 02.d in PLINK-format bfile. We have a set of neutral SNPs in Hardy-Weinberg equilibrium.
# To run flashPCA, we need to run this program at the flashPCA folder.
$PROG/flashPCA/flashpca --bfile $IN/autosomes.canfam.pruned --numthreads 8
	
# We send files to storage folder
cp p*txt $OUT
cp e*txt $OUT

# We change names so that they do not overlap
mv $OUT/pve.txt $OUT/pve.canfam.autosomes.pruned.txt
mv $OUT/pcs.txt $OUT/pcs.canfam.autosomes.pruned.txt
mv $OUT/eigenvalues.txt $OUT/eigenvalues.canfam.autosomes.pruned.txt
mv $OUT/eigenvectors.txt $OUT/eigenvectors.canfam.autosomes.pruned.txt

# Once there, we need to check the pcs file and edit it:

# First, we paste 1st two columns (with species names and individual names) to the table.
paste cols.txt pcs.canfam.autosomes.pruned.txt > pasted.pcs.canfam.autosomes.pruned.txt
			
# Second, we set all multiple spaces as one
sed -n 's/ \+/ /gp' pasted.pcs.canfam.autosomes.pruned.txt > sspaced.pcs.canfam.autosomes.pruned.txt
			
# Third, we concatenate the full header to the file.
cat header.txt sspaced.pcs.canfam.autosomes.pruned.txt > complete.pcs.canfam.autosomes.pruned.txt

# Fourth, we change all spaces for tabs
sed -i 's/\s/\t/g' complete.pcs.canfam.autosomes.pruned.txt 

# Finally, we plot the PCA in R. 

################
## REFERENCES ##
################

# Abraham, G., & Inouye, M. (2014). Fast principal component analysis of large-scale genome-wide data. PLoS ONE, 9(4), 1–5. https://doi.org/10.1371/journal.pone.0093766
