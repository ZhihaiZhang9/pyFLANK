Overview

pyFLANK is an open-source and automated Python implementation which detects FST outliers using a null distribution inferred from quasi-independent loci inspired by the R package OutFLANK(https://doi.org/10.1086/682949). Our tool integrates three approaches to identify loci obeying a null distribution: graph neural network (GNN) inference, linkage disequilibrium (LD)-based inference, and user-defined input. Because pyFLANK uses GNN-based inference of quasi-independent loci, it yields a more accurate null model with less need for user parameter input. 

FST calculation is based on Weir and Cockerham (1984).

Key Features

Graph-based representation of local dependency context of loci and their dependence structure
GNN-based null distribution inference
Compatible with standard FST-based workflows
Designed to complement, not replace, LD pruning and clumping

Date Format

pyFLANK requires two datasets ("Import" model requires prior nuetural loci dataset in VCF format):

1. VCF (v4) file

Only diploid biallelic data is supported for current version.
Compressed (vcf.gz) or not compressed vcf files are both supported.
Phased or unphased genotype are both supported.

2. Population information

A population file is also required as following format:

id          pop
sample1     pop1
sample2     pop1
sample3     pop2
sample4     pop3



Usage


1. Usage for the "Import" model
#"Import" model requires a prior neutral variants to import.
#Usage:

python pyFLANK.py -vcf sim1a.vcf.gz -pop pop1.txt -nim Import -neu which_pruned.vcf.gz -o sim1a_import

"-neu" option used to import the prior neutral variants in VCF format as "-vcf".

2. Usage for the "LD" model
#using "LD" model to infer the neutral variants
python pyFLANK.py -vcf sim1a.vcf.gz -pop pop1.txt -nim LD [-ldwin 1000 -ldcutoff 0.1 -m Bonferroni] -o sim1a_ld



3. Usage for the "GNN" model
#using "GNN" model to infer the neutral variants

python pyFLANK.py -vcf sim1a.vcf.gz -pop pop1.txt -nim GNN [-gnnwin 1000 -loscutoff 0.1 -m Bonferroni]  -o sim1a_gnn

-loscutoff option used to set the loss function threshold to terminate the deep learning epoches.


Interpretation of Results
Users should interpret pyFLANK results conservatively: the method aims to improve calibration and reduce redundant outlier clustering, rather than dramatically increasing power relative to established pipelines.


This repository is actively maintained, and bug reports and feature requests are welcome via GitHub Issues.

Citation

"pyFLANK, a graph neural network based null distribution inference model for FST outlier detection" (under review)
