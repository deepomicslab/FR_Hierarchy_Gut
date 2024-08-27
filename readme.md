# Script description

## Package requirement  

PySEAT may have conflict with numpy version. We recommand: numpy = 1.22.4 and pyseat = 0.0.1.3

To reproduce the result, first please download and unzip data.zip to /data directory.

If you want to start the analysis from GCN, please run /script_procedure/step1_compute_distance first to compute distance matrix of GCN. This may take some time. To save time, you can directly use sp_d.tsv in /data directory which is preproduced by step1_compute_distance.

Prior structure is generated by scripts in /script_GCN_d3/GCN_tree, please run this script before FMT, NSCLC, Anti_analysis which depend on the prior sturcture.

## input: following files (can be found at data/)  

**&ensp;metadata.tsv**(metadata from GutMeta, disease should be in the header for phenotype comparison)  
**&ensp;abd.tsv** (header: species, index: sample name)  
read by abd_profile.input_profile()  
**&ensp;GCN.tsv** (header: KO, index: species)  
read by GCN.input_GCN()  

## <span id="abd">script of taxonomy abundance difference (script_abundance_check)</span>  

Check the abundance difference for each taxon (including NAFLD 16s OTU).  

## script of completeness (script_completeness)  

1. completeness  
Compute the module completeness of each taxon (including NAFLD 16s OTU).

2. test_diff  
Test completeness enrichment based one the GCN prior GCN structure and **should be used after GCN_tree script in script_GCN_d3**

## script of prior GCN structure (script_GCN_d3)  

1. <span id="tree">GCN_tree</span>  
Make prior GCN tree structure.  

2. SE_diff / NFR_diff  
Check SE/nFR difference of disease and health group.  

3. distribution_se  
Plot SE distribution for disease and health group.  

## script of FMT (FMT)  

[**GCN_tree result is required**](#tree)  

Scripts related to two FMT dataset analysis.

1. analysis_se  
Mutiple regression on SE value, days after FMT and fraction at each cluster/super-cluster.  

2. analysis_nfr  
Mutiple regression on nFR, days after FMT and fraction at each cluster/super-cluster.

## script of Antibiotic treatment (script_Anti_analysis)

[**GCN_tree result is required**](#tree)  

1. analysis_se/analysis_nfr  
Check SE/nFR difference of control and exposed group at each clsuter/super-cluter.  

2. analysis_se_exposed/analysis_nfr_exposed  
Check SE/nFR difference of six participants exhibited a bloom of the opportunistic pathogen Enterobacter cloacae complex at the E7 timepoint in exposed group and control group at each clsuter/super-cluter.  

3. merge  
Merge and plot the difference test result of nFR and SE in control and exposed group.

4. merge_exposed  
Merge and plot the difference test result of nFR and SE of the six samples and control group.

5. boxplot
Draw boxplot for SE at each cluster/super-cluster.

## script of lCFR procedure （script_procedure）  

1. step0_NAFLD  
An example of comparing keyston clusters of taxa on NAFLD dataset.  
  
2. step1_compute_distance  
An example of computing taxa distance and KO distance from GCN.
  
3. step2_cluster_analysis  
An example of analyzing keystone cluster and keystone taxon for metagenomics abundance profiles in cMD by constructing posterior structure.  

4. step3_count_support  
An example of checking valid keystone-taxon enterotype with more than one network of size larger than 10 supporting.

5. utils  

    **a. log_effect**  
    An example of computing lCFR and showing distribution of lCFR values and CFR values without log and normalization.

    **b. nestedness_experiment**  
    An example to test the nestedness compared with NULL experiments of lCFR.

    **c. evaluation**  
    An example to evaluate the feature of GCN.  

6. draw_*  
Scripts used to plot the result of previous step.  

## script of NSCLC (script_NSCLC)  

1. SE  
Test difference of SE between response group and non-response group at each cluster/super-cluster and compute FR S score for each sample.

2. sig_SE  
Test difference of SE between response group and non-response group at SIG1/SIG2 clsuter raised in original study and compute S score for each sample.

3. distribution
Plot SE distribution for response group and non-response group.  

4. combination  
Compute combined S score for each sample.

5. The r script
Used to produce the analysis in original study and is provided by https://github.com/valerioiebba/TOPOSCORE/tree/main.  

## plot keystone graph (script_keystone_graph)  

[**GCN_tree result is required**](#tree)  
[**abundance difference result is required**](#abd)  

1. run.ipynb  
Plot keystone result.

## script of eigenspecies analysis (script_eigen_graph_preservation)  

[**GCN_tree result is required**](#tree)  

1. run.ipynb  
Find eigen species and plot the result.
