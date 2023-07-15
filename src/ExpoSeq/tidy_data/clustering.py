from scipy.sparse import csr_matrix
import numpy as np
import editdistance
import networkx as nx

def cleaning(report, max_ld, min_ld):
    aa = list(report["aaSeqCDR3"])
    num_strings = len(aa)
    string_index = [*range(0, num_strings, 1)]
    max_ld = max_ld
    min_ld = min_ld
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
    G = nx.convert_matrix.from_scipy_sparse_matrix(mat)
    degree_sequence = sorted([d for n, d in G.degree()], reverse=True)
    return G, degree_sequence

