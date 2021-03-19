#!/bin/sh

#######################################################################################
## Estimating divergence with PSMC-based Migration and Split Time Inference (MiSTI)  ##
############### Replicability test: c. Cavalli-Sforza (1969) equation #################
#######################################################################################

# In this script, we are estimating Watterson's theta and genomewide Fst to apply at the Cavalli-Sforza's  equation (1969) and compare with MiSTI results. 
# The equation is (equation 4 from Supplementary File: Replicability of results in divergence time estimation. C. Comparison to Cavalli-Sforza estimation.")
# t = (- log (1-ḞST))*2*(θW/(4*μ)): 
# where:
# ḞST is a mean estimate of Fst (with genic regions)
# θW is genomewide Watterson's theta (with genic regions)
# μ is mutation rate per site and generation (we will use μ=4.0*10^-9 as in Skoglund et al., (2015) and 
# μ=4.5*10^-9 as in Frantz et al., (2016) and Koch et al., (2019))
# t will be the estimation of divergence time. 

# We will take into consideration the wide standard deviation of this estimate (Nielsen et al., 1998)

#################################
## Estimate Watterson's theta ###
#################################

# 1st step: we will estimate genomewide θW
# To this end, we will follow the same steps as in script 03.SFS.Het.Fst.thetas.sh (section: thetas). 
# Thetas will be computed per each pair of genomes. We will treat each pair as a separate population to
# estimate pseudo- population effective sizes. The example will be "algeria - egypt"
# Genic regions will be included.

#1. Estimate pestPG
$PROG/angsd/angsd -P 10 -bam $IN/algeria.egypt.autosomes.bamlist -anc $ANC -ref $REF -out $OUT/algeria.egypt.autosomes \
	-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
	-minMapQ 20 -minQ 20 -minInd 2 -setMinDepth 5 -setMaxDepth 56 -doCounts 1 \
	-GL 1 -doSaf 1 -doThetas 1 -pest $OUT/algeria.egypt.autosomes.sfs

#2. Print thetas
$PROG/angsd/misc/thetaStat print $OUT/algeria.egypt.autosomes.wgenic.thetas.idx > $OUT/algeria.egypt.autosomes.wgenic.thetas.print

#3. Check window-based thetas
$PROG/angsd/misc/thetaStat do_stat $OUT/algeria.egypt.autosomes.wogenic.thetas.idx
$PROG/angsd/misc/thetaStat do_stat $OUT/algeria.egypt.autosomes.wogenic.thetas.idx -win 50000 -step 50000 -outnames $OUT/algeria.egypt.autosomes.wogenic.thetas.out

# We will then plot the distribution of sites per window and eliminate those windows where we have less than 5% of the curve distribution. 
# In bash:
cut -f2-14 afr_wolf.algeria.afr_wolf.egypt.wgenic.thetas.ww.out | tail -n+2 > algeria.egypt.thetas

head algeria.egypt.thetas
chr01	225000	38.040209	40.647157	11.641463	60.203793	50.425475	0.715373	1.933220	1.853722	-1.556514	1.146807	40910

# The desired column is 13. 
# In R:

>col = read.table("algeria.egypt.test", as.is=T)
>hist(col$V13, breaks=100)
>quantile(col$V13, probs = c(0.01,0.05,0.99, 0.995)); 
      1%       5%      99%    99.5% 
17300.10 38125.60 48447.88 48627.97 

# We repeat this command for every pair of genomes and extracted those windows who had more than the lowest 5% of data. 
# After this step, we estimate mean and standard deviation in R using commands: mean(col$V13); stdev (col$V13)

##############################
## Estimate genomewide Fst ###
##############################

# 2nd step: we will estimate genomewide Fst
# To this end, we will follow the same steps as in script 03.SFS.Het.Fst.thetas.sh (section: Fst). 
# Fst will be computed per each pair of genomes. The example will be "algeria - egypt"
# Genic regions will be included.

# We first compute 2DSFS files using SAF files.

$PROG/angsd/misc/realSFS -P 16 $OUT/afr_wolf.algeria.autosomes.wgenic.saf.idx $OUT/afr_wolf.egypt.autosomes.wgenic.saf.idx > $OUT/afr_wolf.algeria.afr_wolf.egypt.2d.sfs

# This 2DSFS file will be used to estimate Fst in a 50kb sliding windows scan with realSFS as in:
# https://github.com/mfumagalli/ngsTools/blob/master/TUTORIAL.md 

$PROG/angsd/misc/realSFS fst index -P 16 $OUT/afr_wolf.algeria.autosomes.wgenic.saf.idx $OUT/afr_wolf.egypt.autosomes.wgenic.saf.idx -sfs $OUT/afr_wolf.algeria.afr_wolf.egypt.2d.sfs -fstout $OUT/afr_wolf.algeria.afr_wolf.egypt.wgenic -whichFST 1

$PROG/angsd/misc/realSFS fst stats2 afr_wolf.algeria.afr_wolf.egypt.wgenic.fst.idx -win 50000 -step 50000 -whichFST 1 > afr_wolf.algeria.afr_wolf.egypt.wgenic.fst.ww.out

# .ww.out files will be extracted and the same process as in previous step will be done:
cut -f2-5 afr_wolf.algeria.afr_wolf.egypt.wgenic.fst.ww.out | tail -n+2 > algeria.egypt.fst.file

# Wanted column is #3. As in previous step, we plot the curve and estimate lowest 5% in R. 
>col = read.table("algeria.egypt.fst.file", as.is=T)
>hist(col$V3, breaks=100)
>quantile(col$V3, probs = c(0.01,0.05,0.99, 0.995)); mean (col$V3); stdev (pbs$V3); 
      1%       5%      99%    99.5% 
15224.32 20375.10 45778.67 46414.33

# We repeat this command for every pair of genomes and extracted those windows who had more than the lowest 5% of data. 

########################################################################
## Estimating divergence times with Cavalli-Sforza's (1969) equation ###
########################################################################

# We used mean estimations of Fst (using 95% top significant data under the curve), and mean Watterson's theta with its standard deviation in an Excel file.
# Cavalli-Sforza (1969) and Nielsen et al. (1998) use ḞST without a standard deviation. This is because Fst can take extremely low and high values at a genomic scale. 
# This happens with genes under selection, telomeric and centromeric regions, regions of low mappability and high repeatibility (SINEs, LINEs) and regions with high homozygosity (ROHs). 
# Therefore, we just use the mean value of Fst (estimated) and take into account stdev of Watterson's theta for the table. 
# μ is mutation rate per site and generation (we will use μ=4.0*10^-9 as in Skoglund et al., (2015) and 
# μ=4.5*10^-9 as in Frantz et al., (2016) and Koch et al., (2019))

# Results are in Table S9. 

################
## REFERENCES ##
################

# Cavalli-Sforza, L. L. (1969). Human diversity. Proc. 12th Int.Congr. Genet. 2:405-416

# Frantz, L.A.F., Mullin, V.E., Pionnier-Capitan, M., Lebrasseur, O., Ollivier, M., Perri, A., Linderholm, A., Mattiangeli, V., Teasdale, M.D., Dimopoulos, E.A. (2016). Genomic and archaeological evidence suggest a dual origin of domestic dogs. Science 352(6290):1228–1231.

# Freedman, A. H., et al.(2014). Genome Sequencing Highlights the Dynamic Early History of Dogs. PLoS Genetics, 10(1). https://doi.org/10.1371/journal.pgen.1004016

# Koch, E. M., et al. (2019). De Novo Mutation Rate Estimation in Wolves of Known Pedigree. Molecular Biology and Evolution, 36(11), 2536–2547. https://doi.org/10.1093/molbev/msz159

# Nielsen, R., Mountain, J.L., Huelsenbeck, J.P., Slatkin, M. (1998). Maximum-Likelihood estimation of population divergence times and population phylogeny in models without mutation. Evolution, 52(3). 1998. pp. 669-677

# Skoglund P, Ersmark E, Palkopoulou E, Dalen L. 2015. Ancient wolf genomereveals an early divergence of domestic dog ancestors and admixture into high-latitude breeds. Curr Biol. 25(11):1515–1519.

