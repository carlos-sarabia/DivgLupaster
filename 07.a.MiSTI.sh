#!/bin/sh

#######################################################################################
## Estimating divergence with PSMC-based Migration and Split Time Inference (MiSTI)  ##
#######################################################################################

# MiSTI is a novel program developed by Vladimir Shchur (https://github.com/vlshchur/MiSTI) that is able to estimate migration 
# rates over time and split time between two populations. 
# To this end, the program requires:
# A PSMC file of each genome to be compared. This was done in step 04.PSMC.sh
# A joint site frequency spectrum (2DSFS) file calculated from both genomes (script 03.SFS.Het.Fst.thetas.sh). 
# (Optional) Migration bands: defined segments of time where we estimate different rates of migration. It can be optimized. 

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

# We convert the 2DSFS file from step 03.SFS.Het.Fst.thetas.sh to MiSTI format. Attention! This 2DSFS MUST include genic regions.
# With this example, we will calculate the split time and migration rates between the Algerian and Kenyan African golden wolves.
# Most of these calculations would normally take minutes or seconds and can be run in any home PC.

$PROG/misti/utils/ANGSDSFS.py $dsfs/afr_wolf.algeria.afr_wolf.kenya.wgenic.sfs afr_wolf.algeria afr_wolf.kenya > algeria.kenya.mi.sfs

# We calculate the list of time segments between the Algerian and Kenyan AGW PSMC files. 
# We define a set of units for canids: mutation rate 4.5e-09 (Koch et al., 2019), generation time 3 years (Freedman et al., 2014)

$PROG/misti/utils/calc_time.py --funits $PROG/misti/misti/setunits.canid.txt $psmc/afr_wolf.algeria.psmc $psmc/afr_wolf.kenya.psmc > timescale.algeria.kenya.txt

# We have calculated the heterozygosity loss per genome due to coverage in step 03.SFS.Het.Fst.thetas.sh.
# Also, the file timescale.algeria.kenya.txt will have a number of defined time segments with certain time points in the past. 
# We have based our knowledge in humid/dry periods in Sahara from the literature (see paper) to define a number of time segments 
# under which there could be more or less migration. To avoid biases, we let the program optimize the migration rates by itself.
# These time segments could have (or not) affected our lineages since they might not have been in the same place. For this reason, 
# we repeated the analyses three times being dynamic with time segments and received the same results for time split.
# We will be using GNU Parallel (Tange, 2018), which can paralellize running softwares in home PCs. 
# An example of one of the runs: 

time parallel --header : -j 20 $misti -uf --bsSize 10 --hetloss 0.21472 0.0 --funits $PROG/misti/setunits.canid.txt $psmc/afr_wolf.algeria.psmc $psmc/afr_wolf.kenya.psmc algeria.kenya.mi.sfs {per} -o algeria.kenya.opt.mi -mi 1 0 3 00.0 1 -mi 2 0 3 00.0 1 -mi 1 4 9 00.0 1 -mi 2 4 9 00.0 1 -mi 1 10 35 00.0 1 -mi 2 10 35 00.0 1 -mi 1 36 37 00.0 1 -mi 2 36 37 00.0 1 -mi 1 38 39 00.0 1 -mi 2 38 39 00.0 1 -mi 1 40 41 00.0 1 -mi 2 40 41 00.0 1 -mi 1 42 43 00.0 1 -mi 2 42 43 00.0 1 -mi 1 44 45 00.0 1 -mi 2 44 45 00.0 1 -mi 1 46 50 00.0 1 -mi 2 46 50 00.0 1 >> output.modelized/algeria.kenya.opt.out ::: per 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50

# With this command, we are defining several time segments: 0-3, 4-9, 10-35, 36-37, 38-39, 40-41, 42-43, 44-45, 46-50. 
# The values in real time of these segments will be given by the file timescale.algeria.kenya.txt

####################################################
## Calculating the split time and migration rates ##
####################################################

# The output algeria.kenya.opt.out shows the different paralellized runs of MiSTI and the log-likelihood values. 
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

# The split time must be between 29802-32281 years ago, but it is impossible to calculate it like this. However, we can plot a 
# scatterplot and adjust a polynomial curve between time and llh (we have included time/10000 and llh/100000 here for the sake of 
# simplification). 
# To do so, we have plotted llh/100000 vs time/10000 in an Excel table and have adjusted a polynomial curve of degree 5. R2 was 
# 0.9945 and the equation 9.556E-04x⁵ – 2.763E-02x⁴ + 2.98E-01x³ – 1.483x² + 3.319x -4.942.
# Finally, to calculate the split time, we have used the Newton-Raphson approach with the help of a website to calculate the points 
# in which the derivative of the polynomial reaches 0 (https:www.symbolab.com). We have taken the 0.1%, 1% and 5% of upper points 
# between 0 and 100000 years ago to estimate a confidence interval. 
# These time points and intervals were plotted with data of d¹⁸O from the NGRIP with glacial/interglacial times (see paper).

# To estimate migration rates, the lowest value of llh was taken and migration rates were observed per time segment. 

################
## REFERENCES ##
################

# Freedman, A. H., et al.(2014). Genome Sequencing Highlights the Dynamic Early History of Dogs. PLoS Genetics, 10(1). https://doi.org/10.1371/journal.pgen.1004016

# Koch, E. M., et al. (2019). De Novo Mutation Rate Estimation in Wolves of Known Pedigree. Molecular Biology and Evolution, 36(11), 2536–2547. https://doi.org/10.1093/molbev/msz159

