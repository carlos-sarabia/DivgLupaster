#!/bin/sh

#######################################################################################
## Estimating divergence with PSMC-based Migration and Split Time Inference (MiSTI)  ##
############### Replicability test: a. High and low genomic coverages #################
#######################################################################################

# MiSTI is a novel program developed by Vladimir Shchur (https://github.com/vlshchur/MiSTI) that is able to estimate migration rates over time and split time between two populations. 
# To this end, the program requires:
# A PSMC file of each genome to be compared. This was done in step 04.PSMC.sh
# A joint site frequency spectrum (2DSFS) file calculated from both genomes (script 03.SFS.Het.Fst.thetas.sh). 
# (Optional) Migration bands: defined segments of time where we estimate different rates of migration. It can be optimized. 

# In this script, we are running the same command between the Algerian and Kenyan African golden wolf genomes. 
# Kenyan AGW (24X) has been downsampled to 15X, 11.2X, 9X and 7X.
# Results are in Table S7.

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

###
# 1st step: convert ANGSD 2DSFS file -> MiSTI 2DSFS file
###
# For each coverage, we convert the 2DSFS file from step 03.SFS.Het.Fst.thetas.sh to MiSTI format.
# Most of these calculations would normally take minutes or seconds and can be run in any home PC.
for cov in 24X 15X 11.2X 9X 7X; do
$PROG/misti/utils/ANGSDSFS.py $dsfs/afr_wolf.algeria.afr_wolf.kenya.$cov.sfs afr_wolf.algeria afr_wolf.kenya.$cov > algeria.kenya.mi.$cov.sfs; done

###
# 2nd step: calculate the list of time segments between the Algerian and Kenyan AGW PSMC files
###
# We define a set of units for canids: mutation rate 4.5e-09 (Koch et al., 2019), generation time 3 years (Freedman et al., 2014)
for cov in 24X 15X 11.2X 9X 7X; do
$PROG/misti/utils/calc_time.py --funits $PROG/misti/misti/setunits.canid.txt $psmc/afr_wolf.algeria.psmc $psmc/afr_wolf.kenya.$cov.psmc > timescale.algeria.kenya.$cov.txt; done

###
# 3rd step: Define a number of common time steps to run between the different coverages. 
###
# We have observed that definition of time segments does not affect the divergence time estimation 
#(see Supplementary File: "Replicability tests: results b. Definition of time segments"). # However, we want our time segments to be as close to each other as possible. 
# Our time steps are like this:

time step	algeria.kenya.7X	algeria.kenya.9X	algeria.kenya.11.2X	algeria.kenya.15X	algeria.kenya.24X
0	0		0		0		0		0
1 	 1237		1487		1551		1776		1919
2 	 2205		2205		2205		2205		2205
3 	 2579		3100		3234		3703		4003
4 	 4037		4585		4585		4585		4585
5 	 4585		4852		5060		5793		6263
6 	 5618		6752		7042		7154		7154
7 	 7154		7154		7154		8060		8715
8 	 7334		8814		9192		9928		9928
9 	 9197		9928		9928		10519		11375
10 	 9928		11052		11525		12921		12921
...
# In this table we see how several time steps coincide (0, 2205, 4585, 7154, 9928), which are the time steps extracted from the Algerian .psmc file. We will use these steps as beginning of every time segment. 

###
# 4th step: Run MiSTI. 
###
# We have calculated the heterozygosity loss per genome due to coverage in step 03.SFS.Het.Fst.thetas.sh. 
# We are taking into account genomic heterozygosity loss using the --hetloss command.
# We will be using GNU Parallel (Tange, 2018), which can paralellize running softwares in home PCs. 
# Our runs:

# Algeria.Kenya.24X
time parallel --header : -j 20 $misti -uf --bsSize 10 --hetloss 0.21472 0.0 --funits $MiSTI/setunits.canid.txt $in/psmc/afr_wolf.algeria.normalX.wgenic.psmc $in/psmc/afr_wolf.kenya.24X.wgenic.psmc $in/misfs/algeria.kenya.24X.mi.sfs {per} -o $in/mi/algeria.kenya.24X.opt.extended.mi -mi 1 0 3 00.0 1 -mi 2 0 3 00.0 1 -mi 1 4 9 00.0 1 -mi 2 4 9 00.0 1 -mi 1 10 35 00.0 1 -mi 2 10 35 00.0 1 -mi 1 36 37 00.0 1 -mi 2 36 37 00.0 1 -mi 1 38 39 00.0 1 -mi 2 38 39 00.0 1 -mi 1 40 41 00.0 1 -mi 2 40 41 00.0 1 -mi 1 42 43 00.0 1 -mi 2 42 43 00.0 1 -mi 1 44 45 00.0 1 -mi 2 44 45 00.0 1 -mi 1 46 50 00.0 1 -mi 2 46 50 00.0 1 ::: per 15 16 17 18 19 20 21 22 23 24 25 >> output.test.a/algeria.kenya.24X.extended.out

# Algeria.Kenya.15X
time parallel --header : -j 20 $misti -uf --bsSize 10 --hetloss 0.21472 0.09886 --funits $MiSTI/setunits.canid.txt $in/psmc/afr_wolf.algeria.normalX.wgenic.psmc $in/psmc/afr_wolf.kenya.15X.wgenic.psmc $in/misfs/algeria.kenya.15X.mi.sfs {per} -o $in/mi/algeria.kenya.15X.opt.extended.mi -mi 1 0 3 00.0 1 -mi 2 0 3 00.0 1 -mi 1 4 9 00.0 1 -mi 2 4 9 00.0 1 -mi 1 10 36 00.0 1 -mi 2 10 36 00.0 1 -mi 1 37 38 00.0 1 -mi 2 37 38 00.0 1 -mi 1 39 40 00.0 1 -mi 2 39 40 00.0 1 -mi 1 41 42 00.0 1 -mi 2 41 42 00.0 1 -mi 1 43 44 00.0 1 -mi 2 43 44 00.0 1 -mi 1 45 46 00.0 1 -mi 2 45 46 00.0 1 -mi 1 47 51 00.0 1 -mi 2 47 51 00.0 1 ::: per 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 >> output.test.a/algeria.kenya.15X.extended.out

# Algeria.Kenya.11.2X
time parallel --header : -j 20 $misti -uf --bsSize 10 --hetloss 0.21472 0.18911 --funits $MiSTI/setunits.canid.txt $in/psmc/afr_wolf.algeria.normalX.wgenic.psmc $in/psmc/afr_wolf.kenya.11.2X.wgenic.psmc $in/misfs/algeria.kenya.11.2X.mi.sfs {per} -o $in/mi/algeria.kenya.11.2X.opt.extended.mi -mi 1 0 3 00.0 1 -mi 2 0 3 00.0 1 -mi 1 4 10 00.0 1 -mi 2 4 10 00.0 1 -mi 1 11 37 00.0 1 -mi 2 11 37 00.0 1 -mi 1 38 39 00.0 1 -mi 2 38 39 00.0 1 -mi 1 40 41 00.0 1 -mi 2 40 41 00.0 1 -mi 1 42 43 00.0 1 -mi 2 42 43 00.0 1 -mi 1 44 45 00.0 1 -mi 2 44 45 00.0 1 -mi 1 46 47 00.0 1 -mi 2 46 47 00.0 1 -mi 1 48 52 00.0 1 -mi 2 48 52 00.0 1 ::: per 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 >> output.test.a/algeria.kenya.11.2X.extended.out

# Algeria.Kenya.9X
time parallel --header : -j 20 $misti -uf --bsSize 10 --hetloss 0.21472 0.21472 --funits $MiSTI/setunits.canid.txt $in/psmc/afr_wolf.algeria.normalX.wgenic.psmc $in/psmc/afr_wolf.kenya.9X.wgenic.psmc $in/misfs/algeria.kenya.9X.mi.sfs {per} -o $in/mi/algeria.kenya.9X.opt.extended.mi -mi 1 0 3 00.0 1 -mi 2 0 3 00.0 1 -mi 1 4 10 00.0 1 -mi 2 4 10 00.0 1 -mi 1 11 38 00.0 1 -mi 2 11 38 00.0 1 -mi 1 39 40 00.0 1 -mi 2 39 40 00.0 1 -mi 1 41 42 00.0 1 -mi 2 41 42 00.0 1 -mi 1 43 44 00.0 1 -mi 2 43 44 00.0 1 -mi 1 45 46 00.0 1 -mi 2 45 46 00.0 1 -mi 1 47 48 00.0 1 -mi 2 47 48 00.0 1 -mi 1 49 53 00.0 1 -mi 2 49 53 00.0 1 ::: per 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 >> output.test.a/algeria.kenya.9X.extended.out

# Algeria.Kenya.7X
time parallel --header : -j 20 $misti -uf --bsSize 10 --hetloss 0.21472 0.23025 --funits $MiSTI/setunits.canid.txt $in/psmc/afr_wolf.algeria.normalX.wgenic.psmc $in/psmc/afr_wolf.kenya.7X.wgenic.psmc $in/misfs/algeria.kenya.7X.mi.sfs {per} -o $in/mi/algeria.kenya.7X.opt.extended.mi -mi 1 0 4 00.0 1 -mi 2 0 4 00.0 1 -mi 1 5 11 00.0 1 -mi 2 5 11 00.0 1 -mi 1 12 40 00.0 1 -mi 2 12 40 00.0 1 -mi 1 41 42 00.0 1 -mi 2 41 42 00.0 1 -mi 1 43 44 00.0 1 -mi 2 43 44 00.0 1 -mi 1 45 46 00.0 1 -mi 2 45 46 00.0 1 -mi 1 47 48 00.0 1 -mi 2 47 48 00.0 1 -mi 1 49 50 00.0 1 -mi 2 49 50 00.0 1 -mi 1 51 55 00.0 1 -mi 2 51 53 00.0 1 ::: per 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 >> output.test.a/algeria.kenya.7X.extended.out

# Every time segment defined starts in a common time step found across all five time scale files.  
# The values in real time of these segments will be given by files timescale.algeria.kenya.$cov.txt; with cov being 24X, 15X, 11.2X, 9X and 7X.

####################################################
## Calculating the split time and migration rates ##
####################################################

# The output algeria.kenya.$cov.extended.out shows the different paralellized runs of MiSTI and the log-likelihood values. 
# This command will show the different log likelihoods. We will choose the one closest to 0. 
grep "llh =" algeria.kenya.opt.out

# Calculating the split time is not straightforward. Normally, a list of split times will look like this:
algeria.kenya.24X				
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

# The split time must be between 29802-32281 years ago, but it is impossible to calculate it like this. 
# We plot a scatterplot and adjust a polynomial curve between time and llh (we have included time/10000 and llh/100000 here for the sake of simplification). 
# It is the same step done as in script 07.a. We extract the polynomial curve and R2 
# We also compare the upper 95%, 99% and 99.9% of the curve maximum between the five estimates.

################
## REFERENCES ##
################

# Freedman, A. H., et al.(2014). Genome Sequencing Highlights the Dynamic Early History of Dogs. PLoS Genetics, 10(1). https://doi.org/10.1371/journal.pgen.1004016

# Koch, E. M., et al. (2019). De Novo Mutation Rate Estimation in Wolves of Known Pedigree. Molecular Biology and Evolution, 36(11), 2536–2547. https://doi.org/10.1093/molbev/msz159

