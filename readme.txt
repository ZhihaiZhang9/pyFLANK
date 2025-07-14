#usage for the "Import" model

python pyFLANK.py -vcf sim1a.vcf.gz -pop pop1.txt -nim Import -neu which_pruned.vcf.gz -o sim1a_Import

#usage for the "LD_prune" model

python pyFLANK.py -vcf sim1a.vcf.gz -pop pop1.txt -nim LD [-ldwin 1000 -ldcutoff 0.1 -m Bonferroni] -o sim1a_LD

#usage for the "GNN" model

python pyFLANK.py -vcf sim1a.vcf.gz -pop pop1.txt -nim GNN [-gnnwin 1000 -loscutoff 0.1 -m Bonferroni]  -o sim1a_gnn