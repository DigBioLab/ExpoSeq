import pandas as pd
import editdistance
import numpy as np

def create_distance_matrix(aa):
    num_sequences = len(aa)
    distance_matrix = np.zeros((num_sequences, num_sequences))  # initialize with zeros

    for i in range(num_sequences):
        for j in range(i+1, num_sequences):
            lev_distance = editdistance.distance(aa[i], aa[j])
            distance_matrix[i, j] = lev_distance
            distance_matrix[j, i] = lev_distance  # The distance matrix is symmetric

    return distance_matrix

def get_clustered_sequences(aa, max_cluster_dist):
    clustered_sequences = set()

    for i in range(len(aa)):
        for j in range(i+1, len(aa)):
            if editdistance.distance(aa[i], aa[j]) < max_cluster_dist:
                clustered_sequences.add(aa[i])
                clustered_sequences.add(aa[j])

    return list(clustered_sequences)


def label_bind_seqs(mix, region_string, aa_clustered, antigens):
    values = mix.loc[mix[region_string].isin(aa_clustered), antigens]
    #values = pd.DataFrame(values, columns = ["binding"])
    filtered_df = values[antigens][values[antigens] > 1].dropna(how='all')
    indices = filtered_df.index.tolist()
    key_sequences = []
    key_values = []
    for i in indices:
        key_sequences.append(mix.iloc[i][region_string])
        max_value = max(mix.iloc[i][column] for column in antigens)
        key_values.append(max_value)
    seq_val = pd.DataFrame([key_sequences, key_values]).T
    seq_val.columns = ["binding_seqs", "binding_values"]
    seq_val = seq_val.sort_values(by = "binding_values", ascending = False)
    seq_val.reset_index(inplace = True)
    binding_seqs = seq_val.binding_seqs.to_list()
    for label in range(len(aa_clustered)):
        label_str = aa_clustered[label]
        if label_str in binding_seqs:
            pos_label = binding_seqs.index(label_str)
            label_str = "(" + str(pos_label + 1) + ") " + label_str
            binding_seqs[pos_label] = label_str
        aa_clustered[label] = label_str
    return aa_clustered, binding_seqs, seq_val

def sort_binding_seqs(binding_seqs, seq_val, ascending=True):

    if ascending:
        zipped = list(zip(binding_seqs, seq_val.binding_values))

        # Sort by binding_values
        zipped_sorted = sorted(zipped, key=lambda x: x[1])

        # Unzip the sorted list
        binding_seqs_sorted, binding_values_sorted = zip(*zipped_sorted)
    else:
        binding_seqs_sorted = binding_seqs
        binding_values_sorted = seq_val.binding_values.to_list()
    return binding_seqs_sorted, binding_values_sorted