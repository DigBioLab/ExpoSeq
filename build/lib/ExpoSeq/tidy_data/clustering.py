import numpy as np
import editdistance
import networkx as nx

def cleaning(report, max_ld, min_ld, region_string):
    G = nx.Graph()
    aa = list(report[region_string])
    num_strings = len(aa)
    string_index = [*range(0, num_strings, 1)]
    max_ld = max_ld
    min_ld = min_ld
    for seq in aa:
        G.add_node(seq)
    for i in range(len(aa)):
        string1 = aa[i]
        index_1 = i
        for index in range(index_1):
            string2 = aa[index]
            distance = editdistance.distance(string1, string2)
            if distance>= min_ld and distance <= max_ld:
                G.add_edge(string1, string2)
    return G
