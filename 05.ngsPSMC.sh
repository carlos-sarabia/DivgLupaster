#!/bin/sh

####################################################################################
## Genotype likelihood-based Pairwise Sequentially Markovian Coalescent (ngsPSMC) ##
####################################################################################

# We are running ngsPSMC to estimate more recent demographic history of African golden wolves - divided per lineage. 
# ngsPSMC is a software developed by Shchur, Korneliussen and Nielsen (2017), that uses the same principles as PSMC, but 
# using inferred genotype likelihoods as input, which decreases the uncertainty of called genotypes in low coverage genomes.
# The program can be found as a part of the ANGSD package in: https://github.com/ANGSD/ngsPSMC

###################################
## Variable and PATHs definition ##
###################################

# PATH where our program executables are stored (can also be defined as variables)
PROG='/path/to/programs/folder'
# PATH to the .bam files. 
IN='/path/to/bamfolder/individual'
# PATH with the fasta reference genome CanFam3.1 - Lindblad-Toh et al., (2005) 
REF='/path/to/ref/CanFam3.1'
# PATH to output folder where our ngsPSMC output will be stored:
OUT='/path/to/ANGSD/output/9.ngsPSMC/'

######################
## ngsPSMC analysis ##
######################

# To run a ngsPSMC analysis, we will follow the tutorial in https://github.com/ANGSD/ngsPSMC
# The first thing to take into consideration is that ngsPSMC is still in beta and some parameters cannot be optimized.
# Therefore, we have worked with it using parameters as close as possible to the tutorial provided, but adapting several values to
# our model animal (African golden wolves). 

# First step. We estimate genotype likelihoods with ANGSD.

$PROG/angsd/angsd -P 14 -i $IN/afr_wolf.algeria.autosomes.bam -ref $REF -dopsmc 1 -out $OUT/afr_wolf.algeria.psmc -gl 1 \
	-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
	-minInd 1 -minIndDepth 5 -minq 20 -minmapq 20

# The output of previous step will serve as input for ngsPSMC. Parameters changed from the tutorial:
# We have estimated individual theta as in script 0.3.SFS.Het.Fst.thetas.sh
# Genomewide rho was estimated from a recombination map from dogs (Auton et al., 2013)
# Time segment parameters were estimated as in PSMC (see script 04.PSMC.sh).
# -init parameters represents initial population effective size, relative to Ne=10000. 
# After some attempts, we observed that the ngsPSMC plot started from 70K years ago. We used the calculated population size 70kyr 
# ago from the PSMC plot (see script 04.PSMC.sh), which was Ne=8800

$PROG/ngspsmc/ngspsmc $OUT/afr_wolf.algeria.autosomes.psmc.psmc.idx -p "1*6+58*1" -dospline 0 -nthreads 38 \
	-nIter 50 -init 0.88 -theta 0.000665 -rho 0.00218063 > $OUT/afr_wolf.algeria.autosomes.thetarho.ngspsmc

# Third: plot the ngsPSMC graph.
# The program requires a .txt file with setting units to be defined. Original parameters were: 
# Initial Ne=8800; mutation rate=4.5e-09 (Koch et al., 2019); binsize=100 (default); generation time=3 years (Freedman et al., 2014)

python3.5 $PROG/ngspsmc/utils/psmc_plot.py -psmc algeria $OUT/afr_wolf.algeria.autosomes.thetarho.ngspsmc 0 0 \
	--funits $PROG/ngspsmc/utils/setunits.algeria.txt --fout afr_wolf.algeria.autosomes.ci.thetarho.ngspsmc.pdf \
	--maxY 1.7

# This command was repeated by changing --maxY value from 1.7 to 0.2 (x10000) to increase the zoom over changes in Ne between 100 
# and 10000 years ago

################
## REFERENCES ##
################

# Auton, A., et al. (2013). Genetic Recombination Is Targeted towards Gene Promoter Regions in Dogs. PLoS Genetics, 9(12). https://doi.org/10.1371/journal.pgen.1003984

# Shchur, V., Korneliussen, T. S., & Nielsen, R. (2017). ngsPSMC: genotype likelihood-based PSMC for analysis of low coverage NGS. Retrieved from https://github.com/ANGSD/ngsPSMC







