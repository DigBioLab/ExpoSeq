
import numpy as np
import editdistance


def create_distance_matrix(aa, max_dist=20):  # max_dist is a large finite value
    num_sequences = len(aa)
    distance_matrix = np.full((num_sequences, num_sequences), max_dist)  # initialize with max_dist

    for i in range(num_sequences):
        for j in range(i+1, num_sequences):
            lev_distance = editdistance.distance(aa[i], aa[j])
            if lev_distance < 2:
                distance_matrix[i, j] = lev_distance
                distance_matrix[j, i] = lev_distance  # The distance matrix is symmetric

    return distance_matrix

def get_clustered_sequences(aa):
    clustered_sequences = set()

    for i in range(len(aa)):
        for j in range(i+1, len(aa)):
            if editdistance.distance(aa[i], aa[j]) < 2:
                clustered_sequences.add(aa[i])
                clustered_sequences.add(aa[j])

    return list(clustered_sequences)
