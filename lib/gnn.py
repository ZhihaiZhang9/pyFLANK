import torch
import torch.nn.functional as F
from torch_geometric.nn import GCNConv, GATConv
from torch_geometric.data import Data
from torch_geometric.utils import to_networkx
import numpy as np
import pandas as pd
import networkx as nx
from scipy.stats import f, ks_2samp
from scipy.spatial.distance import squareform
import allel

#from sklearn.metrics import mutual_info_score
#from itertools import combinations

#define the GNN model

class GNNModel(torch.nn.Module):
    def __init__(self, in_channels, out_channels):
        super(GNNModel, self).__init__()
        self.conv1 = GCNConv(in_channels, 16) #Graph Convolution Layer
        self.conv2 = GCNConv(16, out_channels) #Second Layer for output embedding
        
    def forward(self, x, edge_index):
        x = self.conv1(x, edge_index)
        x = F.relu(x)     #rectified linear unit activation function
        x = self.conv2(x, edge_index)
        return x



def calculate_distance(embeddings):
    independent_nodes = []
    independent_nodes.append(embeddings[0])  #selecte the first node
    distances_list = []
    for emb in embeddings[1:]:
        #calculate minimum Euclidean distance of emb to current node
        distances = [np.linalg.norm(emb - node) for node in independent_nodes]
        #print(distances)
        #if min(distances) > threshold:
        #    independent_nodes.append(emb)
        distances_list.append(distances)
    #return independent_nodes
    return distances_list


def GNN_model(genotype_matrix, position_info, gnn_window_bp = 200000, loss_threshold = 0.1):    
    


    #define nodes and edges based on LD (or correlation threashold)

    cor_threashold = 0.2

    
    
    # Ensure the genotype data is a numpy array
    if isinstance(genotype_matrix, pd.DataFrame):
        genotype_matrix = genotype_matrix.values
    
    position_info_array = np.array(position_info, dtype=object)
    
    
    chr = position_info_array[:, 0].astype(str)
    
    
    #genotype_matrix_keep_list = []
    
    
    independent_all_indices = []
    
    for c in np.unique(chr):
        
        #print(c)
        
        #c = '1'
        #c='Chr01'
        
        positions = position_info_array[chr == c][:, 1].astype(int)
        
        genotype_c = genotype_matrix[chr == c]
        embeddings_list = []
        i = 0
        while i < len(positions):
            
            pos = positions[i]
            #num_variants, _ = genotype_matrix.shape
            
            
            print(f'Processing variant position: {pos}')
            
            
            
            #convert window bp to variant count
            
            
            #define the ld window range
            win_start = pos
            win_end = pos + gnn_window_bp
            win_mask = (positions >= win_start) & (positions <= win_end)
            
            #genotype_win = genotype_matrix[win_mask, :] 
            genotype_win = genotype_c[win_mask, :]
            
            #idxs_in_window = np.where(win_mask)[0] # get original indices of window
            
            max_var_count = np.sum(win_mask)
            print(f'Containing {max_var_count} variations within {gnn_window_bp/1000:.0f}-KB span')
            
            
            
            if genotype_win.shape[0] > 1:
                
                correlation_matrix = np.corrcoef(genotype_win) ** 2

                edge_indices = np.argwhere(np.triu(correlation_matrix, k = 1) > cor_threashold)

                edges = edge_indices.tolist()
                
                #convert edges to PyTorch tensor format
                if len(edges) == 0:
                    #i += max_var_count
                    #continue
                    num_nodes = genotype_win.shape[0]
                    edge_index = torch.arange(0, num_nodes, dtype = torch.long).repeat(2,1)
                
                else:
                
                    edge_index = torch.tensor(edges, dtype = torch.long).t().contiguous()


                #create node features ()
                node_features = torch.tensor(genotype_win.mean(axis = 1).reshape(-1, 1), dtype = torch.float)

                #define the graph data structure

                data = Data(x = node_features, edge_index = edge_index)
                
                #instantiate model, optimizer, and loss function

                model = GNNModel(in_channels = 1, out_channels = 1)
                optimizer = torch.optim.Adam(model.parameters(), lr = 0.01) #adaptive moment estimation, to handle sparse gradients and avoid getting stuck in local minima. model.parameters() passes the parameters, weights and biases of the model to the optimizer, so it can adjust them during training to minimize the loss.

                #train loop
                    
                epoch = 0
                while True:
                    model.train()
                    optimizer.zero_grad()
                    out = model(data.x, data.edge_index)
                    
                    
                    #Laplacian-style smoothness regularization
                    edge_weights = torch.pow(out[data.edge_index[0]] - out[data.edge_index[1]], 2).sum()
                    loss = edge_weights
                    loss.backward() #conputes the gradients of the loss with respect to each model parameter.
                    optimizer.step()
                    
                    #if epoch % 100 == 0:
                    #    print(f'Epoch {epoch}, Loss: {loss.item()}')
                    print(f'Epoch {epoch + 1}, Loss: {loss.item()}')
                    
                        
                    if loss.item() < loss_threshold:
                        break
                    
                    epoch += 1
                
                #Node selection based on embedding to find independent subset

                embeddings = model(data.x, data.edge_index).detach().numpy().flatten()
                embeddings_list.append(embeddings)
                
                i += max_var_count
                            
            else:
                
                i += 1
            
        embeddings_all_np = np.concatenate(embeddings_list)
        distances = calculate_distance(embeddings_all_np)
        threshold = np.percentile(distances, 85)
        independent_nodes = [i for i in range(len(distances)) if distances[i] > threshold]
        
        independent_positions = positions[independent_nodes]
        
        independent_df = pd.DataFrame({
        'CHROM': c,
        'POS': independent_positions
        })
        
        independent_all_indices.append(independent_df)
    
    independent_all_df = pd.concat(independent_all_indices, ignore_index=True)
    
    return independent_all_df