#Fst_Diploid

import numpy as np
import pandas as pd

def wc_fst_diploids_2alleles(sample_mat):
    """
    Calculates FST with and without correction for local sample sizes, for diploid biallelic data.
    Based on Weir and Cockerham (1984).
    
    :param sample_mat: A numpy array where each row represents a population, 
                       with columns representing homozygotes for one allele, heterozygotes, and homozygotes for the other allele.
    :return: A dictionary containing He, FST, T1, T2, FSTNoCorr, T1NoCorr, T2NoCorr, and meanAlleleFreq.
    """
    
    sample_sizes = np.sum(sample_mat, axis=1)
    n_ave = np.mean(sample_sizes)
    n_pops = sample_mat.shape[0]
    r = n_pops
    n_c = (n_pops * n_ave - np.sum(sample_sizes**2) / (n_pops * n_ave)) / (n_pops - 1)
    p_freqs = (sample_mat[:, 0] + sample_mat[:, 1] / 2) / sample_sizes
    p_ave = np.sum(sample_sizes * p_freqs) / (n_ave * n_pops)

    s2 = np.sum(sample_sizes * (p_freqs - p_ave)**2) / ((n_pops - 1) * n_ave)
    
    if s2 == 0:
        #return 0
        #optimize the above command:
        return {key: np.nan for key in ["He", "FST", "T1", "T2", "FSTNoCorr", "T1NoCorr", "T2NoCorr", "meanAlleleFreq"]}

    h_freqs = sample_mat[:, 1] / sample_sizes
    h_ave = np.sum(sample_sizes * h_freqs) / (n_ave * n_pops)

    a = n_ave / n_c * (s2 - 1 / (n_ave - 1) * (p_ave * (1 - p_ave) - (r - 1) / r * s2 - 1 / 4 * h_ave))
    b = n_ave / (n_ave - 1) * (p_ave * (1 - p_ave) - (r - 1) / r * s2 - (2 * n_ave - 1) / (4 * n_ave) * h_ave)
    c = 1 / 2 * h_ave

    a_no_corr = n_ave / n_c * s2
    b_no_corr = p_ave * (1 - p_ave) - (r - 1) / r * s2 - 2 * n_ave / (4 * n_ave) * h_ave
    c_no_corr = 1 / 2 * h_ave

    He = 1 - np.sum([p_ave**2, (1 - p_ave)**2])

    FST = a / (a + b + c)
    FSTNoCorr = a_no_corr / (a_no_corr + b_no_corr + c_no_corr)

    return {
        "He": He,
        "FST": FST,
        "T1": a,
        "T2": a + b + c,
        "FSTNoCorr": FSTNoCorr,
        "T1NoCorr": a_no_corr,
        "T2NoCorr": a_no_corr + b_no_corr + c_no_corr,
        "meanAlleleFreq": p_ave
    }
    




def get_fsts_diploids(pop_name_list, snp_data_column):
    """
    Calculates FST etc. from a single locus from a column of individual data.
    
    :param pop_name_list: List of population names.
    :param snp_data_column: Column of individual SNP data (0, 1, 2 or 9 for missing data).
    :return: Dictionary with FST calculations.
    """
    
    
    popnames = np.array(pop_name_list) #numpy.ndarray   (1056, 1)
    
    # Removing missing data
    valid_indices = np.where(snp_data_column != 9)
    #valid_indices = np.unique(np.where(snp_data_column != 9))
    pop_name_temp = popnames[valid_indices] #numpy.ndarray
    pop_name_temp_flat = pop_name_temp.flatten()
    snp_data_temp = snp_data_column[valid_indices[0]] #snp_data_column: pandas.core.frame.DataFrame; snp_data_temp: pandas.core.frame.DataFrame
    snp_data_temp_flat = snp_data_temp.values.flatten()
    
    #het_counts = pd.crosstab(pop_name_temp, snp_data_temp).fillna(0).values
    #het_counts = pd.crosstab(pd.Series(pop_name_temp.flatten()), pd.Series(snp_data_temp.values.flatten())).fillna(0).values
    
    #het_counts = pd.crosstab(pd.Series(pop_name_temp_flat), pd.Series(snp_data_temp_flat)).fillna(0).values
    
    #optimize above command:
    genotypes = pd.Series(snp_data_temp_flat, dtype="category")
    genotypes = genotypes.cat.set_categories([0, 1, 2])
    pop_series = pd.Series(pop_name_temp_flat)
    
    het_counts_df = pd.crosstab(pop_series, genotypes).reindex(columns=[0, 1, 2], fill_value=0)
    het_counts = het_counts_df.values
    
    

    # Case: all individuals are genetically identical at this locus
    if het_counts.shape[1] == 1:
        return {key: np.nan for key in ["He", "FST", "T1", "T2", "FSTNoCorr", "T1NoCorr", "T2NoCorr", "meanAlleleFreq"]}
    
    #fill the missed genotype
    
    if het_counts.shape[1] == 2:
        if "".join(map(str, het_counts.shape)) == "01":
            het_counts = np.hstack([het_counts, np.zeros((het_counts.shape[0], 1))])
        elif "".join(map(str, het_counts.shape)) == "12":
            het_counts = np.hstack([np.zeros((het_counts.shape[0], 1)), het_counts])
        elif "".join(map(str, het_counts.shape)) == "02":
            het_counts = np.hstack([het_counts[:, 0:1], np.zeros((het_counts.shape[0], 1)), het_counts[:, 1:2]])
            
    return wc_fst_diploids_2alleles(het_counts)