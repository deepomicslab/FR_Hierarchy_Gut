# Script description

## Package requirement  

PySEAT may have conflict with numpy version. We recommand: numpy = 1.22.4 and pyseat = 0.0.1.3

## input: following files (can be found at data/)  

**&ensp;metadata.tsv**(metadata from GutMeta, disease should be in the header for phenotype comparison)  
**&ensp;abd.tsv** (header: species, index: sample name)  
read by abd_profile.input_profile()  
**&ensp;GCN.tsv** (header: KO, index: species)  
read by GCN.input_GCN()  

## script of taxonomy abundance difference (script_abundance_check)  

Check the abundance difference for each taxon (including NAFLD 16s OTU).  

## script of completeness (script_completeness)  

Compute the module completeness of each taxon (including NAFLD 16s OTU).

## script of prior GCN structure (script_GCN_d3)  

1. GCN_tree  
Make prior GCN tree structure.  

2. SE_diff / NFR_diff  
Check SE/nFR difference of disease and health group.  

3. tree_annotation  
Make annotation file for itol to show different cohort on the prior tree.  

## script of FMT (FMT)  

**GCN_tree result is required**  

1. analysis_se  
Mutiple regression on SE value, days after FMT and fraction at each cluster/super-cluster.  

2. analysis_nfr  
Mutiple regression on nFR, days after FMT and fraction at each cluster/super-cluster.

## script of Antibiotic treatment (script_Anti)

**GCN_tree result is required**  

1. analysis_se/ analysis_nfr  
Check SE/nFR difference of control and exposed group at each clsuter/super-cluter.  

2. merge  
Merge and plot the result of nFR and SE in control and exposed group.

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
