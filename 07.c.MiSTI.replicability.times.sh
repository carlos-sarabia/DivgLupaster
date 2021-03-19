#!/bin/sh

#######################################################################################
## Estimating divergence with PSMC-based Migration and Split Time Inference (MiSTI)  ##
################# Replicability test: b. Definition of time segments ##################
#######################################################################################

# MiSTI is a novel program developed by Vladimir Shchur (https://github.com/vlshchur/MiSTI) that is able to estimate migration rates over time and split time between two populations. 
# To this end, the program requires:
# A PSMC file of each genome to be compared. This was done in step 04.PSMC.sh
# A joint site frequency spectrum (2DSFS) file calculated from both genomes (script 03.SFS.Het.Fst.thetas.sh). 
# (Optional) Migration bands: defined segments of time where we estimate different rates of migration. It can be optimized. 

# In this script, we are running the same command between the Algerian and Egyptian golden wolf genomes. 
# We defined different time segments, using both humid/dry time segments and paired (even and uneven) time segments. 
# Results are in Table S8.

###################################
## Variable and PATHs definition ##
###################################

# PATH where our program executables are stored (can also be defined as variables)
PROG='/path/to/programs/folder'
# PATH to the .bam files. 
IN='/path/to/bamfolder/individual'
# PATH with the fasta reference genome CanFam3.1 - Lindblad-Toh et al., (2005) 
REF='/path/to/ref/CanFam3.1'
# PATH to output folder where our PSMC input is stored:
psmc='/path/to/ANGSD/output/8.PSMC/PSMC'
# PATH to output folder where our 2DSFS input is stored:
dsfs='/path/to/ANGSD/output/7.SFS/2DSFS'

####################
## MiSTI analysis ##
####################

# We used the same files at every time (afr_wolf.algeria and afr_wolf.egypt) and repeated the previous steps from scripts 7.a and 7.b.
# Heterozygosity loss is taken into account with command --hetloss

# humid/dry
time parallel --header : -j 10 $misti -uf --bsSize 10 --hetloss 0.21472 0.09886 --funits $PROG/MiSTI/setunits.canid.txt $in/psmc/afr_wolf.algeria.normalX.wgenic.psmc $in/psmc/afr_wolf.egypt.normalX.wgenic.psmc $in/misfs/algeria.egypt.mi.sfs {per} -o $in/mi/algeria.egypt.opt.mi -mi 1 0 2 00.0 1 -mi 2 0 2 00.0 1 -mi 1 3 6 00.0 1 -mi 2 3 6 00.0 1 -mi 1 7 13 00.0 1 -mi 2 7 13 00.0 1 -mi 1 14 17 00.0 1 -mi 2 14 17 00.0 1 -mi 1 18 19 00.0 1 -mi 2 18 19 00.0 1 -mi 1 20 22 00.0 1 -mi 2 20 22 00.0 1 -mi 1 23 29 00.0 1 -mi 2 23 29 00.0 1 -mi 1 30 32 00.0 1 -mi 2 30 32 00.0 1 -mi 1 33 34 00.0 1 -mi 2 33 34 00.0 1 -mi 1 35 36 00.0 1 -mi 2 35 36 00.0 1 -mi 1 37 38 00.0 1 -mi 2 37 38 00.0 1 -mi 1 39 40 00.0 1 -mi 2 39 40 00.0 1 -mi 1 41 44 00.0 1 -mi 2 41 44 00.0 1 >> algeria.egypt.opt.out ::: per 3 4 5 6 7 8 9 10 11

# paired even
time parallel --header : -j 20 $misti -uf --bsSize 10 --hetloss 0.21472 0.09886 --funits $PROG/MiSTI/setunits.canid.txt  $in/psmc/afr_wolf.algeria.normalX.wgenic.psmc $in/psmc/afr_wolf.egypt.normalX.wgenic.psmc $in/misfs/algeria.egypt.mi.sfs {per} -o $in/mi/algeria.egypt.opt.mi -mi 1 0 3 00.0 1 -mi 2 0 3 00.0 1 -mi 1 4 5 00.0 1 -mi 2 4 5 00.0 1 -mi 1 6 7 00.0 1 -mi 2 6 7 00.0 1 -mi 1 8 9 00.0 1 -mi 2 8 9 00.0 1 -mi 1 10 11 00.0 1 -mi 2 10 11 00.0 1 -mi 1 12 13 00.0 1 -mi 2 12 13 00.0 1 -mi 1 14 15 00.0 1 -mi 2 14 15 00.0 1 -mi 1 16 17 00.0 1 -mi 2 16 17 00.0 1 -mi 1 18 19 00.0 1 -mi 2 18 19 00.0 1 -mi 1 20 21 00.0 1 -mi 2 20 21 00.0 1 -mi 1 22 23 00.0 1 -mi 2 22 23 00.0 1 -mi 1 24 25 00.0 1 -mi 2 24 25 00.0 1 -mi 1 26 27 00.0 1 -mi 2 26 27 00.0 1 -mi 1 28 29 00.0 1 -mi 2 28 29 00.0 1 -mi 1 30 31 00.0 1 -mi 2 30 31 00.0 1 -mi 1 32 33 00.0 1 -mi 2 32 33 00.0 1 -mi 1 34 35 00.0 1 -mi 2 34 35 00.0 1 -mi 1 36 37 00.0 1 -mi 2 36 37 00.0 1 -mi 1 38 39 00.0 1 -mi 2 38 39 00.0 1 -mi 1 40 41 00.0 1 -mi 2 40 41 00.0 1 -mi 1 42 43 00.0 1 -mi 2 42 43 00.0 1 -mi 1 44 45 00.0 1 -mi 2 44 45 00.0 1 >> algeria.egypt.paired.even.out ::: per 3 4 5 6 7 8 9 10 11

# paired uneven
time parallel --header : -j 20 $misti -uf --bsSize 10 --hetloss 0.21472 0.09886 --funits $PROG/MiSTI/setunits.canid.txt  $in/psmc/afr_wolf.algeria.normalX.wgenic.psmc $in/psmc/afr_wolf.egypt.normalX.wgenic.psmc $in/misfs/algeria.egypt.mi.sfs {per} -o $in/mi/algeria.egypt.opt.paired.uneven.mi -mi 1 0 2 00.0 1 -mi 2 0 2 00.0 1 -mi 1 3 4 00.0 1 -mi 2 3 4 00.0 1 -mi 1 5 6 00.0 1 -mi 2 5 6 00.0 1 -mi 1 7 8 00.0 1 -mi 2 7 8 00.0 1 -mi 1 9 10 00.0 1 -mi 2 9 10 00.0 1 -mi 1 11 12 00.0 1 -mi 2 11 12 00.0 1 -mi 1 13 14 00.0 1 -mi 2 13 14 00.0 1 -mi 1 15 16 00.0 1 -mi 2 15 16 00.0 1 -mi 1 17 18 00.0 1 -mi 2 17 18 00.0 1 -mi 1 19 20 00.0 1 -mi 2 19 20 00.0 1 -mi 1 21 22 00.0 1 -mi 2 21 22 00.0 1 -mi 1 23 24 00.0 1 -mi 2 23 24 00.0 1 -mi 1 25 26 00.0 1 -mi 2 25 26 00.0 1 -mi 1 27 28 00.0 1 -mi 2 27 28 00.0 1 -mi 1 29 30 00.0 1 -mi 2 29 30 00.0 1 -mi 1 31 32 00.0 1 -mi 2 31 32 00.0 1 -mi 1 33 34 00.0 1 -mi 2 33 34 00.0 1 -mi 1 35 36 00.0 1 -mi 2 35 36 00.0 1 -mi 1 37 38 00.0 1 -mi 2 37 38 00.0 1 -mi 1 39 40 00.0 1 -mi 2 39 40 00.0 1 -mi 1 42 43 00.0 1 -mi 2 42 43 00.0 1 -mi 1 44 45 00.0 1 -mi 2 44 45 00.0 1 >> algeria.egypt.paired.uneven.out ::: per 3 4 5 6 7 8 9 10 11


####################################################
## Calculating the split time and migration rates ##
####################################################

# The output algeria.egypt.opt.out shows the different paralellized runs of MiSTI and the log-likelihood values. 
# This command will show the different log likelihoods. We will choose the one closest to 0. 
grep "llh =" algeria.kenya.opt.out

# Calculating the split time is not straightforward. Normally, a list of split times will look like this:
algeria.kenya				
splitT	time		Time/10000	Llh /100000		llh
10	14014.906	1.4014906	-2.4498167494272	-244981.67494272
11	16454.32606	1.6454326	-2.4141074259614	-241410.74259614
12	17095.6004	1.70956004	-2.40589688625513	-240589.688625513
13	20439.20546	2.043920546	-2.36995939498768	-236995.939498768
14	20568.79202	2.056879202	-2.36876460523784	-236876.460523784
15	24068.03546	2.406803546	-2.35175114554301	-235175.114554301
16	25009.44351	2.500944351	-2.34804439531577	-234804.439531577
17	28006.614	2.8006614	-2.34162141334184	-234162.141334184
18	29802.81405	2.9802814	-2.33895852365319	-233895.852365319
19	32281.0112	3.22810112	-2.3387775032911	-233877.75032911
20	34976.85229	3.49768522	-2.34017491043334	-234017.491043334
21	36920.16933	3.69201693	-2.34303327439572	-234303.327439572
22	40561.62956	4.05616295	-2.35082403745713	-235082.403745713

# We repeat the steps of scripts 7.a and 7.b to estimate divergence time, adjusting a polynomial curve between time and llh (we have included time/10000 and llh/100000 here for the sake of simplification) and comparing 95%, 99% and 99.9% of the curve maximum between the three estimates.
# We also added an estimate of curve maximum with humid/dry periods and all time steps between 1 and 44 (150kyr ago). 

################
## REFERENCES ##
################

# Freedman, A. H., et al.(2014). Genome Sequencing Highlights the Dynamic Early History of Dogs. PLoS Genetics, 10(1). https://doi.org/10.1371/journal.pgen.1004016

# Koch, E. M., et al. (2019). De Novo Mutation Rate Estimation in Wolves of Known Pedigree. Molecular Biology and Evolution, 36(11), 2536–2547. https://doi.org/10.1093/molbev/msz159

