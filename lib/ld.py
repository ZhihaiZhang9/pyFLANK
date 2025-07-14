import allel
import numpy as np
import pandas as pd
from scipy.spatial.distance import squareform



def LD_model(genotype_matrix, position_info, ld_window_bp = 100000, ld_threshold = 0.2):
    
    # Ensure the genotype data is a numpy array
    if isinstance(genotype_matrix, pd.DataFrame):
        genotype_matrix = genotype_matrix.values
    
    position_info_array = np.array(position_info, dtype=object)
    
    chr = position_info_array[:, 0].astype(str)
    
    
    #genotype_matrix_keep_list = []
    prune_all_indices = []
    
    for c in np.unique(chr):
        positions = position_info_array[chr == c][:, 1].astype(int)
        
        genotype_c = genotype_matrix[chr == c]
        prune_list = []
        i = 0
        while i < len(positions):
            
            pos = positions[i]
            
            print(f'Processing variant position: {pos}')
            
            
            #define the ld window range
            win_start = pos
            win_end = pos + ld_window_bp
            win_mask = (positions >= win_start) & (positions <= win_end)
            
            genotype_win = genotype_c[win_mask, :] 
            
            
            idxs_in_window = np.where(win_mask)[0] # get original indices of window
            
            max_var_count = np.sum(win_mask)
            print(f'Containing {max_var_count} variations within {ld_window_bp/1000:.0f}-KB span')
            
            if genotype_win.shape[0] > 1:
                r = allel.rogers_huff_r(genotype_win)
                #r = allel.rogers_huff_d(genotype_data)
                #r2 = r ** 2
                r2 = squareform(r ** 2)  #the diagonal will fill with '0'
                #np.fill_diagonal(r2,1)
                
                prune = np.all(r2 < ld_threshold, axis=1)
                
                selected_idxs = idxs_in_window[prune]
                
                #keep_list.append(keep)
                prune_list.extend(selected_idxs)
                
                i += max_var_count
            
            else:
                
                i += 1
            
        positions_prune = positions[prune_list]
        
        
        prune_df = pd.DataFrame({
            'chr': c,
            'position': positions_prune
        })
        
        prune_all_indices.append(prune_df)
        
    #prune_all = np.concatenate(prune_all_indices)
    
    #prune_all_df = pd.DataFrame(prune_all)
    
    prune_all_df = pd.concat(prune_all_indices, ignore_index=True)
    
    return prune_all_df