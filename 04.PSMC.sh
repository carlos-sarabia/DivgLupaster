#!/bin/sh

#######################################################
## Pairwise Sequentially Markovian Coalescent (PSMC) ##
#######################################################

# We are running PSMC to estimate deep demographic history of African golden wolves - divided per lineage. 
# PSMC is a software developed by Li and Durbin, (2011), which is based in the calculation of coalescence events between two 
# alleles at every locus. To achieve this, we must call genotypes with the bcftools software: 
# http://samtools.github.io/bcftools/bcftools.html

###################################
## Variable and PATHs definition ##
###################################

# PATH where our program executables are stored (can also be defined as variables)
PROG='/path/to/programs/folder'
# PATH to the .bam files. 
IN='/path/to/bamfolder/individual'
# PATH with the fasta reference genome CanFam3.1 - Lindblad-Toh et al., (2005) 
REF='/path/to/ref/CanFam3.1'
# PATH to output folder where our PSMC output will be stored:
OUT='/path/to/ANGSD/output/8.PSMC/'
# PATH with a set of genomic regions at 10-kb distances upwards and downwards from genes. 
SITES='/path/to/neutral/regions/'

###################
## PSMC analysis ##
###################

# To run a PSMC analysis, we will follow the tutorial in https://github.com/lh3/psmc
# Calling genotypes in low coverage data might be tricky for low coverage samples as tested in Nadachowska-Brzyska et al., (2013)
# We will correct our PSMC plots using the False Negative Rate (FNR) option of PSMC. 

# First. We call a consensus sequence considering alleles with a minimum of 4X and a maximum of twice the average mean coverage. 
# Using the example of the Algerian African golden wolf:

$PROG/bcftools mpileup -C50 -f $REF -Ou $IN/afr_wolf.algeria.autosomes.bam | $PROG/bcftools call -c | $PROG/bcftools/bin/vcfutils.pl vcf2fq -d 4 -D 23 | gzip > $OUT/afr_wolf.algeria.fq.gz

# We have now a whole genome diploid consensus sequence of the African golden wolf. 
# We use this consensus sequence to run PSMC. We have used the settings of Freedman et al., (2014) with wild canids:
# 64 time intervals divided in 1 time segments of 6 intervals and 58 individual time intervals - "1*6+58*1"
# Also, we previously downsampled the Kenyan (24X) genome to the same coverage of the Algerian individual and visually corrected
# the plot to calculate the FNR.

$PROG/psmc/utils/fq2psmcfa -q20 $OUT/afr_wolf.algeria.fq.gz > $OUT/afr_wolf.algeria.psmcfa

# We bootstrap 50 times the PSMC analysis. 
for i in {01..50}; 
do $PROG/psmc/psmc -N20 -t10 -r6.3291 -b -p "1*6+58*1" -o $OUT/afr_wolf.algeria.round.$i.psmc $OUT/afr_wolf.algeria.split.psmcfa
done

$PROG/psmc/psmc -N20 -t10 -r6.3291 -p "1*6+58*1" -o $OUT/afr_wolf.algeria.unsplit.psmc $OUT/afr_wolf.algeria.psmcfa

# All psmc files are concatenated for the final plot.
cat $OUT/afr_wolf.algeria.unsplit.psmc $OUT/afr_wolf.algeria.round.*.psmc > $OUT/afr_wolf.algeria.combined.psmc

# Finally, we generate the plot:
# We used a generation time of 3 years (Freedman et al., 2014) and a mutation rate of 4.5*10-9 (Koch et al., 2019)
# -N is the FNR correction. 

$PROG/psmc/psmc_plot.pl -P "right bottom" -Y16 -f Helvetica,12 -u 4.5e-09 -g 3 -w 2 -N 0.21 "algeria.corrected" algeria $OUT/afr_wolf.algeria.combined.50.psmc

################
## REFERENCES ##
################

# Freedman, A. H., et al.(2014). Genome Sequencing Highlights the Dynamic Early History of Dogs. PLoS Genetics, 10(1). https://doi.org/10.1371/journal.pgen.1004016

# Koch, E. M., et al. (2019). De Novo Mutation Rate Estimation in Wolves of Known Pedigree. Molecular Biology and Evolution, 36(11), 2536–2547. https://doi.org/10.1093/molbev/msz159

# Nadachowska-Brzyska, K., et al. (2013). Demographic Divergence History of Pied Flycatcher and Collared Flycatcher Inferred from Whole-Genome Re-sequencing Data. PLoS Genetics, 9(11). https://doi.org/10.1371/journal.pgen.1003942





