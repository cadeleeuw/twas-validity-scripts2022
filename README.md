# twas-validity-scripts2022

This repository contains the simulation and analysis scripts for the paper "Reconsidering the validity of transcriptome-wide association studies" ([preprint](https://www.biorxiv.org/content/10.1101/2021.08.15.456414v1)). Content by directory is as follows:

**main:** scripts for real data analysis
- run_bivar.r: runs LAVA-rG and LAVA-TWAS
- run_fusion.bs: runs FUSION

**lava:** scripts for parts of LAVA analysis external to LAVA R package
- lava_twas.r: LAVA-TWAS code
- lava_approx.r: normal approximation for LAVA-rG bivariate analysis

**simulations:** scripts for generation and analysis of simulated data
- generate.r: generates input for all the analysis (both raw E and Y, and corresponding SNP summary statistics)
	- there are 10 different data blocks (from chromosomes 1-10), results are aggregated over thes blocks
	- each call generates 100 iterations for the specified phenotype sample size and number of SNPs
	- simulated values are generated for the full range of (relevant) heritability and expression sample size values
- run_comm.r: run CoMM analysis for a specific pair of conditions at a time 
	- uses raw E and Y
- run_lava.r: run LAVA-TWAS analysis for all pairs of the list of input conditions (for specific number of SNPs in the block)
	- all analyses use the N=1000 genotype block as reference data, regardless of the sample sizes used when generating the data
	- uses summary statistics for both E and Y
- fusion_weights.bs: fits prediction models (BLUP, BSLMM, Elastic Net, LASSO) on the simulated expression
	- to make simulations computationally more feasible, cross-validation is turned off, and true r2 values used in the simulation process are provided to FUSION
	- uses raw E
- fusion_analysis.bs: runs TWAS analysis with fitted prediction models, for all combinations of expression and outcome conditions for given outcome sample size
	- all analyses use the N=1000 genotype block as reference data, regardless of the sample sizes used when generating the data
	- uses summary statistics for Y
- fusion (directory): contains three R helper scripts used in fusion_weights.bs and fusion_analysis.bs

