
from scipy.sparse import csr_matrix
from networkx.convert_matrix import from_scipy_sparse_matrix
import numpy as np
import editdistance
## this script generates too many figures and should be run one library at a time
## If i am tu use these graphs
## Should nodesize be equal to identical CDR3Hs?

# imnet process_strings.distance matrix from https://github.com/rokroskar/imnet/blob/master/imnet/process_strings.py


def cleaning(report):
    aa = list(report["aaSeqCDR3"])
    num_strings = len(aa)
    string_index = [*range(0, num_strings, 1)]
    max_ld = 2
    min_ld = 0
    collect_coordinates = np.array([0, 0, 0]).reshape(1,3)
    for i in range(len(string_index)):
        string1 = aa[i]
        index_1 = i
        for index in range(index_1):
            string2 = aa[index]
            distance = editdistance.distance(string1, string2)
            if distance>= min_ld and distance <= max_ld:
                collect_coordinates = np.concatenate((collect_coordinates,
                                                      np.array([index_1, index, distance]).reshape(1, 3)))

    mat = csr_matrix((collect_coordinates[:, 2],
                      (collect_coordinates[:, 0], collect_coordinates[:, 1])),
                     (len(aa), len(aa)))
    G = from_scipy_sparse_matrix(mat)
    degree_sequence = sorted([d for n, d in G.degree()], reverse=True)
    return G, degree_sequence

