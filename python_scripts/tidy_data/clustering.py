#import imnet
from scipy.sparse import csr_matrix
from networkx.convert_matrix import from_scipy_sparse_matrix
from numpy import array, empty
import matplotlib.pyplot as plt
import editdistance
## this script generates too many figures and should be run one library at a time
## If i am tu use these graphs
## Should nodesize be equal to identical CDR3Hs?

# imnet process_strings.distance matrix from https://github.com/rokroskar/imnet/blob/master/imnet/process_strings.py

#def cleaning(sample_name, report):
 #   sample = report[report["Experiment"] == sample_name]
  #  aa = list(sample["aaSeqCDR3"])
   # mat_arr = array(list(imnet.process_strings.distance_matrix(aa, min_ld=0, max_ld=2)))
    #mat = csr_matrix((mat_arr[:, 2], (mat_arr[:, 0], mat_arr[:, 1])), (len(aa), len(aa)))
    #G = from_scipy_sparse_matrix(mat)
    # nx.draw_networkx(G_raph, node_size=10, with_labels=False)
    #degree_sequence = sorted([d for n, d in G.degree()], reverse=True)
   # dmax = max(degree_sequence)
    #return G, degree_sequence

def cleaning(sample_name, report):
    sample = report[report["Experiment"] == sample_name]
    aa = list(sample["aaSeqCDR3"])
    num_strings = len(aa)
    # Create an empty distance matrix
    distance_matrix = empty((num_strings, num_strings))
    # Fill the distance matrix with the Levenshtein distances between the strings,
    # setting a maximum distance of 2 and a minimum distance of 0
    max_distance = 2
    min_distance = 0
    for i in range(num_strings):
        for j in range(num_strings):
            distance = editdistance.distance(aa[i], aa[j])
            distance_matrix[i, j] = min(max(distance, min_distance), max_distance)

    mat = csr_matrix((distance_matrix[:, 2], (distance_matrix[:, 0], distance_matrix[:, 1])), (len(aa), len(aa)))
    G = from_scipy_sparse_matrix(mat)
    # plt.figure(figsize=(12,12))
    # nx.draw_networkx(G_raph, node_size=10, with_labels=False)
    degree_sequence = sorted([d for n, d in G.degree()], reverse=True)
   # dmax = max(degree_sequence)
    return G, degree_sequence

