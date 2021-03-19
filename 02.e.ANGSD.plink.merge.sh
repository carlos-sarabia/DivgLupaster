#!/bin/sh

##############################################
## Merge PLINK output and convert to bfiles ##
##############################################
# We use the PLINK-format output of step 02.b (section: SNP discovery) and merge all chromosomal SNP files.
# Also, we filter out: SNPs under linkage disequilibrium (LD), SNPs not in Hardy-Weinberg Equilibrium (HWE), 
# SNPs with a minor allele frequency (MAF) less than 5% and missing genotypes. 

###################################
## Variable and PATHs definition ##
###################################
# PATH where our program executables are stored (can also be defined as variables)
PROG='/path/to/programs/folder'
# PATH to folder where our input (the output of step 02.b) is stored:
IN='/path/to/ANGSD/output/1.angsd/'
# PATH to folder where our output will be stored:
OUT='/path/to/ANGSD/output/4.PLINK/'

###############################
## SNP merging and filtering ##
###############################

# The output of step 02.b will come in PLINK transposed file format (tfiles). So, we need to convert them to bfile format.
# Since PLINK can also work with lists of data, we will store the file names in several lists before the next step.

for i in {01..38}; do
  	$PROG/plink/plink  --tfile $IN/chr$i.canfam --allow-no-sex --dog --make-bed --out $OUT/chr$i.canfam
        echo $OUT/chr$i.canfam.bed >> $OUT/bed.canfam.list
        echo $OUT/chr$i.canfam.bim >> $OUT/bim.canfam.list
        echo $OUT/chr$i.canfam.fam >> $OUT/fam.canfam.list
done

paste $OUT/bed.canfam.list $OUT/bim.canfam.list $OUT/fam.canfam.list > $OUT/dataset.canfam.txt

$PROG/plink/plink --bfile $OUT/chr01.canfam --merge-list $OUT/dataset.canfam.txt --dog --make-bed --out $OUT/autosomes.canfam

# We now filter SNPs in LD, HWE, MAF and missing genotypes. We handle canid data, so we set up option --dog.
# We filter SNPs in linkage disequilibrium in windows of 50 SNPs with a step size of 5 SNPs and a R2 of 0.5.

$PROG/plink/plink --bfile $OUT/autosomes.canfam --allow-no-sex --dog --make-bed --missing-genotype N --indep-pairwise 50 5 0.5 --hwe 0.001 --maf 0.05 --out $OUT/autosomes.canfam.pruned

# More details on this data type and how to work with them in https://zzz.bwh.harvard.edu/plink/data.shtml






