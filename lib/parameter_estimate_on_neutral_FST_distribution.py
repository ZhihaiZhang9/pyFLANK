import numpy as np
import pandas as pd
from scipy.stats import chi2
from statsmodels.stats.multitest import fdrcorrection
from statsmodels.stats.multitest import multipletests
from scipy.optimize import minimize
from scipy.special import gammaincc, gamma



# Mean FstNoCorr calculation from DataList
def fstBarCalculatorNoCorr(DataList):
    # Calculates mean FstNoCorr from the dataframe, using sum(T1NoCorr) / sum(T2NoCorr)
    return np.sum(DataList['T1NoCorr'][~DataList['OutlierFlag']]) / np.sum(DataList['T2NoCorr'][~DataList['OutlierFlag']])

# Mean Fst calculation from DataList
def fstBarCalculator(DataList):
    # Calculates mean Fst from the dataframe, using sum(T1) / sum(T2)
    return np.sum(DataList['T1'][~DataList['OutlierFlag']]) / np.sum(DataList['T2'][~DataList['OutlierFlag']])

'''
# Effective number of samples maximum likelihood estimate
def effective_number_samples_mle(FstVect, Fstbar, NumberOfSamples, SmallestFstInTrimmedList, LargestFstInTrimmedList):
    """
    This function finds the maximum likelihood value of the effective number of samples 
    for a given list of Fst values.
    
    FstVect should already be purged of NaN values and loci with too low heterozygosity or MAF.
    """
    
    #test
    
    FstVect = working_dataframe[putative_neutral_list_temp]['FSTNoCorr']
    Fstbar = Fstbar_no_corr_temp
    NumberOfSamples = number_of_samples
    SmallestFstInTrimmedList = low_trim_point
    LargestFstInTrimmedList = high_trim_point
    
    
    # Sort the Fst vector
    sorted_fst = np.sort(FstVect)
    
    # Minimum Fst considered in the trimmed data
    low_trim_point = max(Fstbar / 100, SmallestFstInTrimmedList)
    
    # Filter the Fst vector within the range
    trimmed_fst_vect = FstVect[(FstVect >= low_trim_point) & (FstVect <= LargestFstInTrimmedList)]
    
    def local_nll_all_data(df_inferred):
        def local_nll_one_locus(Fst):
            return neg_ll_df_fst_trim(Fst, df_inferred, Fstbar, low_trim_point, LargestFstInTrimmedList)
        
        return np.sum([local_nll_one_locus(Fst) for Fst in trimmed_fst_vect])
    
    # Optimization using L-BFGS-B method
    result = minimize(local_nll_all_data, NumberOfSamples, bounds=[(2, None)], method="L-BFGS-B")
    
    return result.x[0]

# Incomplete gamma function
def incomplete_gamma_function(a, z):
    """
    This is equivalent to Mathematica's Gamma[a, z] according to 
    http://r.789695.n4.nabble.com/Incomplete-Gamma-function-td833545.html
    """
    return gammaincc(a, z) * gamma(a)

# Negative log likelihood for Fst with trimming
def neg_ll_df_fst_trim(Fst, df_inferred, Fstbar, low_trim_point, high_trim_point):
    """
    Finds contribution to the negative log likelihood of a given locus' Fst for a given df_inferred.
    This is based on chi-squared distribution of neutral Fst.
    """
    #test
    Fst = 0.023219
    df_inferred = df_inferred_temp
    Fstbar = Fstbar_no_corr_temp
    
    
    df = df_inferred
    term1 = df * Fst / (2 * Fstbar)
    term2 = df * Fstbar * np.log(2) / (2 * Fstbar)
    term3 = -df * Fstbar * np.log(df) / (2 * Fstbar)
    term4 = -(df - 2) * Fstbar * np.log(Fst) / (2 * Fstbar)
    term5 = df * Fstbar * np.log(Fstbar) / (2 * Fstbar)
    term6 = 2 * Fstbar * np.log(-incomplete_gamma_function(df / 2, df * high_trim_point / (2 * Fstbar)) + incomplete_gamma_function(df / 2, df * low_trim_point / (2 * Fstbar)))
    
    return term1 + term2 + term3 + term4 + term5 + term6
'''


def effective_number_samples_mle(FstVect, Fstbar, NumberOfSamples, SmallestFstInTrimmedList, LargestFstInTrimmedList):
    """
    This function finds the maximum likelihood value of the effective number of samples 
    for a given list of Fst values.
    
    FstVect should already be purged of NaN values and loci with too low heterozygosity or MAF.
    """
    
    # Minimum Fst considered in the trimmed data
    low_trim_point = max(Fstbar / 100, SmallestFstInTrimmedList)
    
    # Filter the Fst vector within the range
    trimmed_fst_vect = FstVect[(FstVect >= low_trim_point) & (FstVect <= LargestFstInTrimmedList)]
    
    # Nested function for calculating the negative log likelihood for a single locus
    def neg_ll_df_fst_trim(Fst, df_inferred):
        df = df_inferred
        
        '''
        term1 = df * Fst / (2 * Fstbar)
        term2 = df * Fstbar * np.log(2) / (2 * Fstbar)
        term3 = -df * Fstbar * np.log(df) / (2 * Fstbar)
        term4 = -(df - 2) * Fstbar * np.log(Fst) / (2 * Fstbar)
        term5 = df * Fstbar * np.log(Fstbar) / (2 * Fstbar)
        
        # Handle potential invalid values for the incomplete gamma function
        low_gamma = incomplete_gamma_function(df / 2, df * low_trim_point / (2 * Fstbar))
        high_gamma = incomplete_gamma_function(df / 2, df * LargestFstInTrimmedList / (2 * Fstbar))

        # Ensure that values are valid for logarithm
        if high_gamma <= 0 or low_gamma <= 0:
            return np.inf  # Return a large value to indicate an invalid likelihood
        
        term6 = 2 * Fstbar * np.log(high_gamma - low_gamma)
        '''
        
        #return term1 + term2 + term3 + term4 + term5 + term6
        return 1/(2*Fstbar)*(df * Fst +df * Fstbar * np.log(2) - df * Fstbar * np.log(df)-(df-2)*Fstbar * np.log(Fst)+df * Fstbar * np.log(Fstbar) + 2*Fstbar * np.log(-incomplete_gamma_function(df/2,df*LargestFstInTrimmedList/(2*Fstbar))+incomplete_gamma_function(df/2,df * low_trim_point/(2*Fstbar))))

    # Total negative log likelihood for all data
    def local_nll_all_data(df_inferred):
        return np.sum([neg_ll_df_fst_trim(Fst, df_inferred) for Fst in trimmed_fst_vect])

    # Optimization using L-BFGS-B method
    result = minimize(local_nll_all_data, NumberOfSamples, bounds=[(2, None)], method="L-BFGS-B")
    
    return result.x[0]

def incomplete_gamma_function(a, z):
    """
    This is equivalent to Mathematica's Gamma[a, z].
    """
    return gammaincc(a, z) * gamma(a)


def p_two_sided_from_chisq(x, df):
    """
    Takes a value x, finds the two-sided p-value for comparison to a chi-square distribution
    with df degrees of freedom.
    """
    p_one_sided = chi2.cdf(x, df)
    p_two_sided = np.where(p_one_sided > 0.5, (1 - p_one_sided) * 2, p_one_sided * 2)
    return p_two_sided

# Chi-squared p-value function for neutrality
def p_outlier_finder_chisq_no_corr(data_list, Fstbar, df_inferred, qthreshold=0.05, Hmin=0.1, mulTestCorrect = 'bonferroni'):
    #Finds outliers based on chi-squared distribution
    #Takes given values of dfInferred and Fstbar, and returns a list of p-values and q-values for all loci based on chi-square.
    #Assumes that the DataList input has a column called $FSTNoCorr and that empty columns exist for $qvalues and $OutlierFlag 
    
    #Divide DataList into 3 lists:  DataListGood has $FST>0; DataListNeg has cases where $FST <=0; and
    #   DataListNA has cases where $FST is NA.
    #DataListNeg is necessary to keep separate here because these cases do not have meaningful results with the chi-square approach;
    #   however, they do carry information.
    
    #test
    '''
    data_list = fst1_df
    Fstbar = null_para['FSTNoCorrbar']
    df_inferred = null_para['dfInferred']
    qthreshold=0.05
    Hmin=0.1
    '''
    
    if 'indexOrder' not in data_list:
        data_list['indexOrder'] = range(len(data_list))
    
    keepers = (data_list['FSTNoCorr'] > 0) & (data_list['He'] >= Hmin)
    data_list_good = data_list[keepers].copy()
    data_list_others = data_list[~keepers].copy()

    # Calculate chi-squared p-values
    p_list = p_two_sided_from_chisq(data_list_good['FSTNoCorr']* df_inferred / Fstbar,df_inferred)
    p_list_right_tail = 1 - chi2.cdf(data_list_good['FSTNoCorr'] * df_inferred / Fstbar, df_inferred)
    
    mulTestCorrect = mulTestCorrect.lower()
    
    if mulTestCorrect == "bonferroni":
        # Adjust p-values using Bonferroni correction
        reject, qtemp, _, _ = multipletests(p_list_right_tail, alpha = 0.05, method = 'bonferroni')
    elif mulTestCorrect == 'fdr':
        # Adjust p-values using FDR
        reject, qtemp = fdrcorrection(p_list_right_tail, alpha=qthreshold, method='indep')
    #else:
    #    raise ValueError("Invalid multiple test correction method. Choose 'bonferroni' or 'fdr'.")
    
    data_list_good['pvalues'] = p_list
    data_list_good['pvaluesRightTail'] = p_list_right_tail
    data_list_good['qvalues'] = qtemp
    data_list_good['OutlierFlag'] = data_list_good['qvalues'] < qthreshold

    # Combining the good and bad loci back
    data_list_combined = pd.concat([data_list_good, data_list_others], axis=0).sort_values(by='indexOrder')
    return data_list_combined




def NullCalibration(fst_dataframe, left_trim_fraction=0.05, right_trim_fraction=0.05, Hmin=0.1, number_of_samples=39, qthreshold=0.05):
    
    #NumberOfSamples is the number of populations sampled
    
    #initialize
    len_df = len(fst_dataframe['FSTNoCorr'])
    
    fst_dataframe.loc[:, 'indexOrder'] = np.arange(len_df)
    fst_dataframe.loc[:, 'GoodH'] = np.where(fst_dataframe['He'] < Hmin, "lowH", "goodH")
    #fst_dataframe.loc[:, 'OutlierFlag'] = np.where(pd.isna(fst_dataframe['FSTNoCorr']), np.nan, 'False')
    fst_dataframe.loc[:, 'OutlierFlag'] = fst_dataframe['FSTNoCorr'].isna()
    fst_dataframe.loc[:, 'qvalues'] = np.full(len_df, np.nan)
    fst_dataframe.loc[:, 'pvalues'] = np.full(len_df, np.nan)
    fst_dataframe.loc[:, 'pvaluesRightTail'] = np.full(len_df, np.nan)
    
    

    #filter and sort
    nonkeepers = (fst_dataframe['FSTNoCorr'].isna()) | (fst_dataframe['He'] < Hmin)
    working_dataframe = fst_dataframe[~nonkeepers].copy().sort_values(by='FSTNoCorr')
    
    stored_data_frame_na = fst_dataframe[nonkeepers]

    n_loci_total = len(working_dataframe)
    
    #check if there are enough loci to proceed
    
    #####added on 02042025
    if n_loci_total == 0:
        print("No valid FST values after filtering, Exiting function.")
        return None
    
    #compute trim indices
    smallest_keeper = int(np.ceil(n_loci_total * left_trim_fraction))
    largest_keeper = int(np.floor(n_loci_total * (1 - right_trim_fraction)))
    
    #####added on 02042025
    #ensure indices are within valid bounds
    smallest_keeper = max(0, min(smallest_keeper, n_loci_total - 1))
    largest_keeper = max(0, min(largest_keeper, n_loci_total -1))
    
    
    low_trim_point = working_dataframe['FSTNoCorr'].iloc[smallest_keeper]
    high_trim_point = working_dataframe['FSTNoCorr'].iloc[largest_keeper]
    
    

    #check trim point
    if low_trim_point < 0:
        raise ValueError("The smallest FST in the trimmed set must be > 0. Please use a larger LeftTrimFraction.")
    if high_trim_point >= 1:
        raise ValueError("The largest FST in the trimmed set must be < 1. Please use a larger RightTrimFraction.")

    #initialize the parameters of iterization
    putative_neutral_list_temp = working_dataframe['FSTNoCorr'] > 0
    old_outlier_flag = np.zeros(n_loci_total, dtype=bool)
    count = 0

    while True:
        count += 1
        if count > 20:
            print("Exceeded iteration maximum.")
            break

        Fstbar_no_corr_temp = fstBarCalculatorNoCorr(working_dataframe[putative_neutral_list_temp])
        
        
        #df_inferred_temp = effective_number_samples_mle(working_dataframe['FSTNoCorr'].loc[putative_neutral_list_temp], Fstbar_no_corr_temp, number_of_samples, low_trim_point,high_trim_point)
        df_inferred_temp = effective_number_samples_mle(working_dataframe[putative_neutral_list_temp]['FSTNoCorr'], Fstbar_no_corr_temp, number_of_samples, low_trim_point,high_trim_point)

        working_dataframe = p_outlier_finder_chisq_no_corr(working_dataframe, Fstbar_no_corr_temp, df_inferred_temp, qthreshold, Hmin)
        
        #working_dataframe.to_csv('working_dataframe.csv')
        
        '''
        #calculate p and q
        p_list_right_tail = 1 - chi2.cdf(working_dataframe['FSTNoCorr'] * df_inferred_temp / Fstbar_no_corr_temp, df_inferred_temp)
        _, qtemp = fdrcorrection(p_list_right_tail, alpha=qthreshold, method='indep')

        #working_dataframe['pvaluesRightTail'] = p_list_right_tail
        working_dataframe.loc[:,'pvaluesRightTail'] = p_list_right_tail
        #working_dataframe['qvalues'] = qtemp
        working_dataframe.loc[:, 'qvalues'] = qtemp
        #working_dataframe['OutlierFlag'] = working_dataframe['qvalues'] < qthreshold
        working_dataframe.loc[:, 'OutlierFlag'] = working_dataframe['qvalues'] < qthreshold
        '''

        # outlier list and neutral loci
        #putative_neutral_list_temp = ~working_dataframe['OutlierFlag']
        putative_neutral_list_temp = working_dataframe['OutlierFlag'] == False
        if putative_neutral_list_temp.sum() == 0:
            print("No loci in neutral list...")
            return None

        if np.array_equal(old_outlier_flag, working_dataframe['OutlierFlag']):
            break

        old_outlier_flag = working_dataframe['OutlierFlag']

    number_low_fst_outliers = sum(working_dataframe['OutlierFlag'] & (working_dataframe['FSTNoCorr'] < low_trim_point))
    number_high_fst_outliers = sum(working_dataframe['OutlierFlag'] & (working_dataframe['FSTNoCorr'] > high_trim_point))

    #FSTbar = working_dataframe.loc[putative_neutral_list_temp, 'FSTNoCorr'].mean()
    FSTbar = fstBarCalculator(working_dataframe[putative_neutral_list_temp])
    
    results_dataframe = pd.concat([working_dataframe, stored_data_frame_na], ignore_index = True)
    results_dataframe = results_dataframe.sort_values(by='indexOrder', ascending = True).reset_index(drop=True)
    
    
    return {
        'FSTbar': FSTbar,
        'FSTNoCorrbar': Fstbar_no_corr_temp,
        'dfInferred': df_inferred_temp,
        'numberLowFstOutliers': number_low_fst_outliers,
        'numberHighFstOutliers': number_high_fst_outliers,
        'results': results_dataframe
    }