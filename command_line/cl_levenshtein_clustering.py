from ExpoSeq.settings.collecting_all_arguments import ExpoSeqArgs, prep_args
from ExpoSeq.plots.levenshtein_clustering import PrepareData, LevenshteinClustering
import pandas as pd
import networkx as nx
import numpy as np 

def call_args():
    Args = ExpoSeqArgs()
    Args.add_save_csv()
    Args.add_sequencing_report()
    Args.add_samples()
    Args.add_batch_size()
    Args.add_levenshtein_distance()
    Args.add_region_plots()
    Args.add_binding_data()
    Args.add_antigen_names()
    Args.add_dimension()
    return Args

if __name__ == "__main__":
    args = call_args()
    parser = prep_args(args)
    sequencing_report = pd.read_csv(parser.sequencing_report)
    if parser.binding_data is not None:
        bind_data = pd.read_csv(parser.binding_data)
    else: 
        bind_data = None
    # plot prepare
    PrepData = PrepareData()
    
    LSCluster = LevenshteinClustering(sequencing_report, 
                                      parser.samples, 
                                      batch_size = parser.batch_size, 
                                      region_string = parser.region_plots, 
                                      max_ld = parser.levenshtein_distance, 
                                      binding_data=bind_data, 
                                      antigens = parser.antigen_names) # get graph object and associated report with data
    G = LSCluster.G
    assert LSCluster.report.shape[0] > 0, "No sequences found for the given sample"
    
    pos = nx.spring_layout(G, dim = parser.dimension, seed = 42) 
    
    nodes = np.array([pos[v] for v in sorted(G)]) # extract nodes from graph object
    
    edges = np.array([(pos[u], pos[v]) for u, v in G.edges()]) # extract edges from graph object
    
    print(edges)
    if parser.dimension == 2:
        data = pd.DataFrame(nodes, columns = ["x", "y"]) # this is the x and y coordinate of the nodes
    else:
        data = pd.DataFrame(nodes, columns = ["x", "y", "z"])
    # if there is no matching edge the corresponding values for matching x and y are nan
    data["matching_x"] = np.nan # this is the x coordinate where the edge should end
    data["matching_y"] = np.nan # this is the y coordinate where the edge should end
    if parser.dimension == 3:
        data["matching_z"] = np.nan # this is the z coordinate where the edge should end
        
    data["nodesize"] = LSCluster.nodesize # nodesize which should be used
    
    data["Sequence"] = LSCluster.report[parser.region_plots].to_list() # add sequences
    
    data["ClusterID"] = data["Sequence"].map(LSCluster.partition) # this is the cluster id of the sequence
    
    for edge_pair in edges:
        if parser.dimension == 2:
            condition = (data["x"] == edge_pair[0, 0]) & (data["y"] == edge_pair[0, 1])
        else:
            condition = (data["x"] == edge_pair[0, 0]) & (data["y"] == edge_pair[0, 1]) & (data["z"] == edge_pair[0, 2])
        if condition.any():
            index = data[condition].index[0]
            data.loc[index, "matching_x"] = edge_pair[0, 0]
            data.loc[index, "matching_y"] = edge_pair[0, 1]
            if parser.dimension == 3:
                data.loc[index, "matching_z"] = edge_pair[0, 2]
        
    data.to_csv(
        parser.save_csv
    )
    
