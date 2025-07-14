import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import chi2
import numpy as np
import pandas as pd
import re

#extract numeric part for proper sorting of chromosomes
def extract_chr_numer(chr_str):
    if isinstance(chr_str, str):
        if 'X' in chr_str: return 1000
        if 'Y' in chr_str: return 1001
        if 'MT' in chr_str or 'M' in chr_str: return 1002
        match = re.search(r'\d+', str(chr_str))
        return int(match.group()) if match else float('inf')
    return 9999

def scatter_plot(plt_title, plotdf, x_info, y_info, diff_points, diff_points_style, marker_size, legend_title = None, lengend_loc = None, highlight=None):
    
    #test
    '''
    plt_title = 'He_FST_check'
    plotdf = fst1_df
    x_info = 'He'
    y_info = 'FST' 
    diff_points = None 
    diff_points_style = None
    marker_size = 100
    legend_title = None
    lengend_loc = None
    highlight = None
    
    '''
    
    
    #set the figure size
    plt.figure(figsize = (10, 6))
    
    #create the scatter plot
    sns.scatterplot(data = plotdf, x = x_info, y = y_info, hue = diff_points, style = diff_points_style, s = marker_size)
    
    # Highlight specific points
    if highlight is not None and not highlight.empty:
        #plt.scatter(x[highlight_indices], y[highlight_indices], color='red', s=100, label='Highlighted Points', edgecolor='black')
        #sns.scatterplot(data = highlight, x = x_info, y = y_info, hue = diff_points, style = diff_points_style, edgecolor = 'red', facecolor = 'none', marker = 'o', s = marker_size + 20)
        #plt.scatter(highlight[x_info], highlight[y_info], facecolors = 'none', edgecolors = 'red', marker = 'o', s = marker_size + 10, label = 'Highlight Points')
        plt.scatter(highlight[x_info], highlight[y_info], facecolors = 'none', edgecolors = 'red', s = marker_size + 10, label = 'Highlight Points')
    
    
    #add labels and titles
    plt.title(plt_title, fontsize = 16)
    plt.xlabel(x_info, fontsize = 14)
    plt.ylabel(y_info, fontsize = 14)
    if legend_title is not None:
        plt.legend(title = legend_title, bbox_to_anchor = (1.5,0.5),
                    loc = lengend_loc, 
                    borderaxespad = 0, ncol = 2)
    #show the plot
    #plt.grid(Ture)
    #plt.show()
    
    #adjust layout to prevent overlap with legend
    plt.tight_layout()
    #plt.show()
    
    #save the plot
    plt.savefig(plt_title, dpi = 300)
    
    #plt.close()


def qq_plot(plt_title, df, x_info, y_info, diff_points, diff_points_style, marker_size, legend_title = None, lengend_loc = None):
    
    #test
    '''
    plt_title = 'FST_FSTNoCorrcheck_WholeProfile'
    df = fst1_df
    x_info = 'FST'
    y_info = 'FSTNoCorr'
    diff_points = None
    diff_points_style = None
    marker_size = 100
    '''
    
    
    scatter_plot(plt_title, df, x_info, y_info, diff_points, diff_points_style, marker_size, legend_title, lengend_loc)
    
    #add the diagonal line y = x
    
    maxY = max(df[x_info].max(), df[y_info].max())
    
    #plt.plot([0, maxY + 0.1], [0, maxY + 0.1], color = 'red', linestyle = '--', linewidth = 2)
    plt.plot([0, maxY], [0, maxY], color = 'red', linestyle = '--', linewidth = 2)
    
    #save
    plt.savefig(plt_title, dpi = 300)
    
    #plt.close()



def FstDistPlotter(df, FSTlist, FSTbar, binwidth=0.001, titletext = None):
    
    #test
    '''
    df = null_para['dfInferred']
    FSTlist = null_para['results']['FSTNoCorr']
    FSTbar = null_para['FSTbar']
    binwidth=0.005
    titletext=None
    '''
    
    xPlotUpperBound = np.ceil(max(FSTlist) * 100) / 100
    breakslist = np.arange(0, xPlotUpperBound + binwidth, binwidth)
    
    y = np.array([chi2.cdf(((i - 0.5) * binwidth) / FSTbar * df, df=df) - chi2.cdf((((i - 1.5) * binwidth)) / FSTbar * df, df=df) for i in range(len(breakslist))])
    y = len(FSTlist) * y
    
    plt.hist(FSTlist, bins=breakslist, color="darkgoldenrod", edgecolor="black")
    plt.plot(breakslist, y, color="darkblue", lw=3)
    plt.xlabel("Fst")
    plt.title(titletext)
    #plt.show()
    
    #save
    #plt.savefig(titletext, dpi = 300)
    
    #plt.close()


def FstDistPlotterZoom(df, FSTlist, FSTbar, binwidth=0.005, titletext = None, RightZoomFraction=0.1):
    
    FSTlistNoNA = FSTlist[~np.isnan(FSTlist)]
    
    xPlotUpperBound = np.ceil(max(FSTlistNoNA) * 100) / 100
    xPlotLowerBound = np.floor(np.quantile(FSTlistNoNA, 1 - RightZoomFraction) * 100) / 100
    flist = FSTlistNoNA[FSTlistNoNA > xPlotLowerBound]
    
    breakslist = np.arange(xPlotLowerBound, xPlotUpperBound, binwidth)
    
    y = np.array([chi2.cdf((xPlotLowerBound + (i - 0.5) * binwidth) / FSTbar * df, df=df) - chi2.cdf((xPlotLowerBound + (i - 1.5) * binwidth) / FSTbar * df, df=df) for i in range(len(breakslist))])
    y = len(FSTlistNoNA) * y
    
    plt.hist(flist, bins=breakslist, color="darkgoldenrod", edgecolor="black")
    plt.plot(breakslist, y, color="darkblue", lw=3)
    plt.xlabel("Fst")
    plt.title(titletext)
    #plt.show()
    
    #save
    #plt.savefig(titletext, dpi = 300)
    
    #plt.close()
    

'''
def FstDistPlotterAddBadCurve(df, FSTlist, FSTbar, binwidth=0.005, RightZoomFraction=0.99):
    
    #test
    df = null_para['dfInferred']
    FSTlist = null_para['results']['FST']
    FSTbar = null_para['FSTbar']
    binwidth=0.001
    RightZoomFraction=0.99
    
    
    FSTlistNoNA = FSTlist[~np.isnan(FSTlist)]
    
    xPlotUpperBound = np.ceil(max(FSTlistNoNA) * 100) / 100
    xPlotLowerBound = np.floor(np.quantile(FSTlistNoNA, 1 - RightZoomFraction) * 100) / 100
    
    breakslist = np.arange(xPlotLowerBound, xPlotUpperBound, binwidth)
    
    y = np.array([chi2.cdf((xPlotLowerBound + (i - 0.5) * binwidth) / FSTbar * df, df=df) - chi2.cdf((xPlotLowerBound + (i - 1.5) * binwidth) / FSTbar * df, df=df) for i in range(len(breakslist))])
    y = len(FSTlistNoNA) * y
    
    plt.plot(breakslist, y, color="red", lw=3)
    plt.show()


def OutFLANKBadCurvePlotter(badDF, OFoutput, withOutliers=True, NoCorr=True, Hmin=0.1, binwidth=0.005, Zoom=False, RightZoomFraction=0.99, titletext=None):
    data = OFoutput['results'][OFoutput['results']['He'] > Hmin]
    
    if NoCorr:
        flist = data['FSTNoCorr']
        fbar = sum(data['T1NoCorr']) / sum(data['T2NoCorr'])
    else:
        flist = data['FST']
        fbar = OFoutput['FSTbar']
    
    flist = flist[~np.isnan(flist)]
    keeperlist = data[~data['OutlierFlag']].index
    
    if not withOutliers:
        flist = flist.loc[keeperlist]
    
    FstDistPlotterAddBadCurve(badDF, flist, fbar, binwidth, RightZoomFraction)
'''    

def OutFLANKResultsPlotter(OFoutput, withOutliers=True, NoCorr=True, Hmin=0.1, binwidth=0.005, Zoom=False, RightZoomFraction=0.05, titletext=None):

    #test
    '''
    OFoutput = null_para
    withOutliers=True
    NoCorr=True
    Hmin=0.1
    binwidth=0.001
    Zoom=False
    RightZoomFraction=0.05 
    titletext=None
    '''

    data = OFoutput['results'][OFoutput['results']['He'] > Hmin]
    
    if NoCorr:
        flist = data['FSTNoCorr']
        fbar = sum(data['T1NoCorr']) / sum(data['T2NoCorr'])
        titletext = f"{titletext} Fst without sample size correction"
    else:
        flist = data['FST']
        fbar = OFoutput['FSTbar']
        titletext = f"{titletext} Fst with sample size correction"
    
    flist = flist[~np.isnan(flist)]
    keeperlist = data[~data['OutlierFlag']].index
    
    if not withOutliers:
        flist = flist.loc[keeperlist]
    
    if Zoom:
        FstDistPlotterZoom(OFoutput['dfInferred'], flist, fbar, binwidth, titletext, RightZoomFraction)
        plt.savefig(f"{titletext} (Right tail zoomed)", dpi = 300)
    else:
        FstDistPlotter(OFoutput['dfInferred'], flist, fbar, binwidth, titletext)
        plt.savefig(f"{titletext}", dpi = 300)



    
def histgram_plot(data, binwd, titletext = None, x_lab = None):
        
    #test
    '''
    data = null_para['results']['pvaluesRightTail']
    binwd = 30
    titletext = 'histogram of pvaluesRightTail'
    x = 'pvaluesRightTail'
    '''
    
    plt.hist(data, bins = binwd, color = 'blue', edgecolor = 'black')
    plt.title(titletext)
    plt.xlabel(x_lab)
    plt.ylabel('Frequency')
        
    #save
    plt.savefig(titletext, dpi = 300)
    
    #plt.close()
        

def manhattan_plot(plotdf, chr_col, bp_col, val_col, qvalues, threshold=None, title=None, LOG = False):

    #test
    '''
    plotdf = fst_final_filtered_He[['CHROM', 'POS', 'FST', 'qvalues']]
    chr_col = 'CHROM'
    bp_col = 'POS'
    val_col = 'FST'
    qvalues = 'qvalues'
    threshold = 0.05
    LOG = True
    '''
    
    plotdf[bp_col] = pd.to_numeric(plotdf[bp_col], errors= 'coerce')
    # Ensure the data is sorted by chromosome and base-pair position
    plotdf = plotdf.sort_values([chr_col, bp_col])

    # Convert chromosome column to categorical (useful for sorting chromosomes numerically)
    #plotdf[chr_col] = plotdf[chr_col].astype('category')
    
    #sort chromosomes by numeric part
    unique_chroms = sorted(plotdf[chr_col].unique(), key = extract_chr_numer)
    #apply ordered categorical with the sorted chromosome list
    plotdf[chr_col] = pd.Categorical(plotdf[chr_col], categories = unique_chroms, ordered = True)
    

    # Create a list of colors for alternating chromosomes
    colors = ['#1f77b4', '#ff7f0e']

    # Create a new figure for the plot
    plt.figure(figsize=(12, 6))

    # Track the base-pair positions for chromosome ticks
    x_labels = []
    x_ticks = []

    # Plot the data, chromosome by chromosome
    last_x = 0
    
    if LOG:
        ## Calculate the -log10 of the p-values for plotting
        plotdf[val_col] = np.where(plotdf[val_col] > 0 , plotdf[val_col], np.nan)
        plotdf['-log10p'] = -np.log10(plotdf[val_col])
        for i, (chromosome, group) in enumerate(plotdf.groupby(chr_col, observed = True)):
            group['ind'] = range(last_x, last_x + len(group))
            plt.scatter(group['ind'], group['-log10p'], c=colors[i % len(colors)], s=10, label=f'Chr {chromosome}')
            if threshold:
                #log_threshold = -np.log10(threshold)
                #find all lines with qvalue <= threshold:
                filtered = plotdf[plotdf[qvalues] <= threshold]
                #find the max target value as the final threshold
                corrected_threshold = filtered[val_col].max()
                log_threshold = -np.log10(corrected_threshold)
                
                sig = group[group['-log10p'] > log_threshold]
                plt.scatter(sig['ind'], sig['-log10p'], color='red', s=10)
                plt.axhline(log_threshold, color='red', linestyle='dashed', label='Threshold' if threshold else None)
            last_x += len(group)
            x_ticks.append((group['ind'].iloc[-1] + group['ind'].iloc[0]) / 2)
            x_labels.append(chromosome)
            
        plt.ylabel(f'-log{val_col}')
    else:
        for i, (chromosome, group) in enumerate(plotdf.groupby(chr_col, observed = True)):
            group['ind'] = range(last_x, last_x + len(group))
            plt.scatter(group['ind'], group[val_col], c=colors[i % len(colors)], s=10, label=f'Chr {chromosome}')
            if threshold:
                #find all lines with qvalue <= threshold:
                filtered = plotdf[plotdf[qvalues] <= threshold]
                #find the max target value as the final threshold
                corrected_threshold = filtered[val_col].min()
                
                #sig = group[group[val_col] > threshold]
                sig = group[group[val_col] > corrected_threshold]
                plt.scatter(sig['ind'], sig[val_col], color='red', s=10)
                #plt.axhline(threshold, color='red', linestyle='dashed', label='Threshold' if threshold else None)
                plt.axhline(corrected_threshold, color='red', linestyle='dashed', label='Threshold' if threshold else None)   
            last_x += len(group)
            x_ticks.append((group['ind'].iloc[-1] + group['ind'].iloc[0]) / 2)
            x_labels.append(chromosome)
            
        plt.ylabel(val_col)

    # Highlight points above the significance threshold, if provided
    
    '''
    if threshold:
        sig = df[df['-log10p'] > -np.log10(threshold)]
        plt.scatter(sig['ind'], sig['-log10p'], color='red', s=10)
    '''
    

    # Add labels and title
    #plt.axhline(-np.log10(threshold), color='red', linestyle='dashed', label='Threshold' if threshold else None)
    
    plt.xlabel('Chromosome')
    plt.xticks(ticks=x_ticks, labels=x_labels)
    plt.title(title if title else 'Manhattan Plot')

    # Display legend
    #plt.legend()

    # Show the plot
    #plt.tight_layout()
    #plt.show()
    
    #save
    
    plt.savefig(title, dpi = 300)
    
    #plt.close()
    
    