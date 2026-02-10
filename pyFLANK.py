#!/usr/bin/python3

import os, sys, argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


import logging

'''
import torch
import torch.nn.functional as F
from torch_geometric.nn import GCNConv, GATConv
from torch_geometric.data import Data
from torch_geometric.utils import to_networkx
import networkx as nx
from scipy.stats import f, ks_2samp
from scipy.spatial.distance import squareform
'''

current_dir = os.path.dirname(__file__)
lib_path = os.path.join(current_dir, 'lib')
sys.path.append(lib_path)

#load build libraries
import vcfReadAndParse
import Fst_diploids
import parameter_estimate_on_neutral_FST_distribution
import visualization
import ld
import gnn



# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')


parser = argparse.ArgumentParser(
    usage='./%(prog)s -vcf vcfFile -pop populationFile -nim ["Import", "LD", "GNN"] [-neu NeutralFile] -o outputfile',
    description = ''' This software implements the method developed by Whitlock and Lotterhos (2015) to use 
    likelihood on a trimmed distribution of FST values to infer the distribution of FST for neutral markers.''',
    
    epilog = """Written by Zhihai Zhang (Sep,04,2024)""",
    )

parser.add_argument(
    '-vcf',
    '--vcffile',
    required = True,
    help = 'input vcf file',
    metavar=''
)


parser.add_argument(
    '-pop',
    '--population',
    required = True,
    help = 'population file',
    metavar=''
)


parser.add_argument(
    '-nim',
    '--neutral_infer_method',
    required = True,
    choices = ['Import', 'LD', 'GNN'],
    help = 'input neutral file, available options: "Import", "LD", "GNN"',
    metavar=''
)

parser.add_argument(
    '-neu',
    '--neutral_file',
    required = False,
    help = 'input neutral file, which is required when set "-nim Import"',
    metavar=''
)


parser.add_argument(
    '-ldwin',
    '--ld_window_bp',
    required = False,
    default= 100000,
    type = int,
    help = 'input ld window distance (bp), which is required when set "-nim LD", default is 100kb',
    metavar=''
)

parser.add_argument(
    '-ldcutoff',
    '--ld_cutoff',
    required = False,
    default= 0.1,
    type = float,
    help = 'input ld cutoff for ld prune, which is required when set "-nim LD", default is 0.1',
    metavar=''
)

parser.add_argument(
    '-gnnwin',
    '--gnn_window_bp',
    required = False,
    default= 200000,
    type = int,
    help = 'input gnn window distance (bp), which is required when set "-nim GNN", default is 200kb',
    metavar=''
)

parser.add_argument(
    '-loscutoff',
    '--loss_threshold',
    required = False,
    default= 0.1,
    type = float,
    help = 'input ld cutoff for ld prune, which is required when set "-nim LD", default is 0.1',
    metavar=''
)


'''
parser.add_argument(
    '-num_epo',
    '--number_of_epoch',
    required = False,
    default = 20,
    type = int,
    help = 'number of epoch to train the model',
    metavar=''
)
'''


parser.add_argument(
    '-m',
    '--multiTestCorrectionMethod',
    required = False,
    default= 'bonferroni',
    help = "Multiple test correction method, choose 'bonferroni' or 'fdr', the default is 'bonferroni'.",
    metavar=''
)



parser.add_argument(
    '-o',
    '--outputFilePrefix',
    required = True,
    help = 'outputfile prefix',
    metavar=''
)


args = parser.parse_args()


#excute

#set the input and output files
vcf_file_path = args.vcffile


pop = pd.read_csv(args.population, sep='\\s+', usecols = [1], skiprows = 1, header = None)

merged_output_file_path = args.outputFilePrefix

#read and extract the genotype
position_info, genotype_info = vcfReadAndParse.read_vcf_file(vcf_file_path)

#format conversion
binary_data = vcfReadAndParse.process_genotype(genotype_info)

G_matrix = pd.DataFrame(binary_data).fillna(0).astype(int)


#calculate FST
fst1 = []
for i in range(len(G_matrix)):

    fst1_line = Fst_diploids.get_fsts_diploids(pop, G_matrix.iloc[i].T)
    chrom, pos = position_info[i]
    fst1.append({"CHROM": chrom, "POS": pos, **fst1_line})
    
fst1_df = pd.DataFrame(fst1)



#check the He and Fst

visualization.scatter_plot(f"{args.outputFilePrefix}_He_FST_check_WholeProfile", fst1_df, 'He', 'FST', None, None, 100, None, None, None)
plt.close()



visualization.scatter_plot(f"{args.outputFilePrefix}_He_FST_check_ChromWise", fst1_df, 'He', 'FST', 'CHROM', 'CHROM', 100, legend_title = 'Chromosome', lengend_loc = 'center right', highlight=None)
plt.close()

#check the FST and FSTNoCorr

visualization.qq_plot(f"{args.outputFilePrefix}_FST_FSTNoCorrcheck_WholeProfile", fst1_df, 'FST', 'FSTNoCorr', None, None, 100, None, None)
plt.close()

visualization.qq_plot(f"{args.outputFilePrefix}_FST_FSTNoCorrcheck_ChromWise", fst1_df, 'FST', 'FSTNoCorr', 'CHROM', 'CHROM', 100, legend_title = 'Chromosome', lengend_loc = 'center right')
plt.close()



if args.neutral_infer_method == 'Import':
    
    if args.neutral_file:
        null_vcf_file_path = args.neutral_file
        if null_vcf_file_path is not None:

            null_position_info, _ = vcfReadAndParse.read_vcf_file(null_vcf_file_path)

            null_pos_df = pd.DataFrame(null_position_info, columns = ['CHROM', 'POS'])
            
            fst1_df['POS'] = fst1_df['POS'].astype(int)
            fst1_df['CHROM'] = fst1_df['CHROM'].astype(str)
            
            null_pos_df['POS'] = null_pos_df['POS'].astype(int)
            null_pos_df['CHROM'] = null_pos_df['CHROM'].astype(str)

            null_fst_df = fst1_df.merge(null_pos_df, on = ['CHROM', 'POS'], how = 'inner')
        else:
            raise ValueError("Error: 'The 'neutral_file' argument must be set when 'neutral_infer_method' is 'Import'.")
    else:
        raise ValueError("Error: 'The 'neutral_file' argument must be set when 'neutral_infer_method' is 'Import'.")

if args.neutral_infer_method == 'LD':
    
    genotype_matrix = G_matrix.to_numpy()
    
    prune_index = ld.LD_model(genotype_matrix, position_info, ld_window_bp = args.ld_window_bp, ld_threshold = args.ld_cutoff)
    
    fst1_df['POS'] = fst1_df['POS'].astype(int)
    fst1_df['CHROM'] = fst1_df['CHROM'].astype(str)
    
    prune_index['position'] = prune_index['position'].astype(int)
    prune_index['chr'] = prune_index['chr'].astype(str)
    
    
    null_fst_df = fst1_df.merge(prune_index, left_on=['CHROM', 'POS'], right_on=['chr','position'], how = 'inner')
    null_fst_df = null_fst_df[null_fst_df["FSTNoCorr"] < 0.1][fst1_df.columns]
    
    if null_fst_df.shape[0] < 30:
        print("The number of loci used for calibrating the null distribution of FST is too small, please adjust the parameters.")
        sys.exit()
    else:
        print(f"The number of loci used for calibrating the null distribution of FST is {null_fst_df.shape[0]}")
        
        
if args.neutral_infer_method == 'GNN':
    
    
    #using deep learning GNN model to infer the null fst df
    genotype_matrix = G_matrix.to_numpy()
    
    independent_index = gnn.GNN_model(genotype_matrix, position_info, gnn_window_bp = args.gnn_window_bp, loss_threshold = args.loss_threshold)
    
    fst1_df['POS'] = fst1_df['POS'].astype(int)
    fst1_df['CHROM'] = fst1_df['CHROM'].astype(str)
    
    independent_index['POS'] = independent_index['POS'].astype(int)
    independent_index['CHROM'] = independent_index['CHROM'].astype(str)
    
    null_fst_df = fst1_df.merge(independent_index, on = ['CHROM', 'POS'], how = 'inner')
    
        
    null_fst_df = null_fst_df[null_fst_df["FSTNoCorr"] < 0.1]
    
    if null_fst_df.shape[0] < 30:
        print("The number of loci used for calibrating the null distribution of FST is too small, please adjust the parameters.")
        sys.exit()
    else:
        print(f"The number of loci used for calibrating the null distribution of FST is {null_fst_df.shape[0]}")
    


##initialize calculation dataframe
#estimate the parameters on the neutral FST distribution

NumPopSampled = pop.iloc[:,0].nunique()

null_para = parameter_estimate_on_neutral_FST_distribution.NullCalibration(null_fst_df, left_trim_fraction=0.05, right_trim_fraction=0.05, Hmin=0.1, number_of_samples=NumPopSampled, qthreshold=0.05)

if null_para is None:
    sys.exit()


#visualize FST distrubtion of neutral sites and p value histogram
visualization.OutFLANKResultsPlotter(null_para, withOutliers=True, NoCorr=True, Hmin=0.1, binwidth=0.001, Zoom=False, RightZoomFraction=0.05, titletext = args.outputFilePrefix)
plt.close()

visualization.OutFLANKResultsPlotter(null_para, withOutliers=True, NoCorr=True, Hmin=0.1, binwidth=0.001, Zoom=True, RightZoomFraction=0.05, titletext = args.outputFilePrefix)
plt.close()

visualization.histgram_plot(null_para['results']['pvaluesRightTail'], 10, titletext = f"{args.outputFilePrefix}_Histogram_of_pvaluesRightTail", x_lab = 'pvaluesRightTail')
plt.close()



#using estimated neutral mean FST and df to calculate P-values for all loci
fst_final = parameter_estimate_on_neutral_FST_distribution.p_outlier_finder_chisq_no_corr(fst1_df, null_para['FSTNoCorrbar'], null_para['dfInferred'], qthreshold=0.05, Hmin=0.1, mulTestCorrect = args.multiTestCorrectionMethod)

fst_final_filtered_outlier = fst_final[fst_final['OutlierFlag'] == True]
fst_final_filtered_He = fst_final[fst_final['He'] > 0.1]

#FST treatment
#convert to numeric
fst_final_filtered_He['FST'] = pd.to_numeric(fst_final_filtered_He['FST'], errors='coerce')

#replace +/- inf with NaN
fst_final_filtered_He['FST'].replace([np.inf, -np.inf], np.nan, inplace=True)

#set negative FST to zero
fst_final_filtered_He.loc[fst_final_filtered_He['FST'] < 0, 'FST'] = 0

#drop rows with NaN / None values in 'FST' column
fst_final_filtered_He = fst_final_filtered_He.dropna(subset=['FST'])



fst_final_filtered_He_outlier = fst_final_filtered_He[fst_final_filtered_He['OutlierFlag'] == True]

#visulation

visualization.scatter_plot(f'{args.outputFilePrefix}_He_FST_check_final_WholeProfile', fst_final_filtered_He, 'He', 'FST', 'CHROM', 'CHROM', 100, legend_title = 'Chromosome', lengend_loc = 'center right', highlight = None)
plt.close()

visualization.scatter_plot(f'{args.outputFilePrefix}_He_FST_check_final_ChromWise', fst_final_filtered_He, 'He', 'FST', 'CHROM', 'CHROM', 100, legend_title = 'Chromosome', lengend_loc = 'center right', highlight = fst_final_filtered_He_outlier)
plt.close()

visualization.manhattan_plot(fst_final_filtered_He[['CHROM', 'POS', 'FST', 'qvalues']], 'CHROM', 'POS', 'FST', 'qvalues', threshold=0.05, title=f"{args.outputFilePrefix}_FST_manhattan", LOG = False)
plt.close()

visualization.manhattan_plot(fst_final_filtered_He[['CHROM', 'POS', 'pvaluesRightTail', 'qvalues']], 'CHROM', 'POS', 'pvaluesRightTail', 'qvalues', threshold=0.05, title=f"{args.outputFilePrefix}_pvaluesRightTail_manhattan", LOG = True)
plt.close()


#save
fst_output1 = fst_final[['CHROM', 'POS', 'He', 'FST', 'T1', 'T2', 'FSTNoCorr', 'T1NoCorr', 'T2NoCorr', 'meanAlleleFreq', 'pvalues', 'pvaluesRightTail', 'qvalues', 'OutlierFlag']]
fst_output1.to_csv(f"{merged_output_file_path}_final.csv", sep= '\t', index = False)

fst_output2 = fst_final_filtered_He[['CHROM', 'POS', 'He', 'FST', 'T1', 'T2', 'FSTNoCorr', 'T1NoCorr', 'T2NoCorr', 'meanAlleleFreq', 'pvalues', 'pvaluesRightTail', 'qvalues', 'OutlierFlag']]

fst_output2.to_csv(f"{merged_output_file_path}_filtered_He.csv", sep= '\t', index = False)

fst_output3 = fst_final_filtered_He[['CHROM', 'POS', 'He', 'FST', 'T1', 'T2', 'FSTNoCorr', 'T1NoCorr', 'T2NoCorr', 'meanAlleleFreq', 'pvalues', 'pvaluesRightTail', 'qvalues', 'OutlierFlag']]

fst_output3.to_csv(f"{merged_output_file_path}_filtered_Outlier.csv", sep= '\t', index = False)

fst_output4 = fst_final_filtered_He_outlier[['CHROM', 'POS', 'He', 'FST', 'T1', 'T2', 'FSTNoCorr', 'T1NoCorr', 'T2NoCorr', 'meanAlleleFreq', 'pvalues', 'pvaluesRightTail', 'qvalues', 'OutlierFlag']]

fst_output4.to_csv(f"{merged_output_file_path}_filtered_He_Outlier.csv", sep= '\t', index = False)
