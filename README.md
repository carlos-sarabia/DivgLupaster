# DivgLupaster

This repository contains all code files released in Sarabia, C.; vonHoldt, B.; Larrasoaña, J.C.; Uríos, V.; Leonard, J.A. Pleistocene climate fluctuations drove demographic history of African golden wolves (Canis lupaster). Molecular Ecology, 00:1-20. doi:10.1111/mec.15784

####

Files 00. to 01. refer to read pre-processing from raw fastq files to mapped, sorted and indel realigned .bam files.

Files 02. describe a series of commands with ANGSD/NGSadmix and exploratory analyses: calculation of a distribution of qscores, calling genotype likelihoods, PCA with genotype likelihoods and/or called SNPs, an Admixture-like plot with genotype likelihoods and/or called SNPs. 

File 03. describes a series of calculations for summary statistics: SFS discovery, genomewide heterozygosity corrected by individual depth of coverage, genomewide Fst and thetas (Watterson, Pairwise thetas - nucleotide diversity, Fu & Li's theta, Fay's H theta) and neutrality tests (Tajima's D, Fu&Li F's, Fu&Li's D, Fay's H, Zeng's E).

File 04. describes the generation of a PSMC plot generation with a bootstrapping of 50 iterations.

File 05. describes the generation of a ngsPSMC plot using parameters as close as possible to the PSMC. 

File 06. is a comparison of methods to estimate individual and population's inbreeding coefficients (Fi), some of which are based in Runs of Homozygosity (ROH) estimations

Files 07. explain the use of MiSTi to calculate divergences between lineages using SFS and PSMC data corrected by coverage. Files 07.b-d test the replicability of the MiSTi test. 

####
