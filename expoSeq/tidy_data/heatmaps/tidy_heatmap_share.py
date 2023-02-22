import numpy as np
from .bool_sequences_matrix import find_seq_matches
import pandas as pd

def heatmap_share(sequencing_report, protein, specific_experiments):
    heatmap_axis, unique_sequences, unique_experiments = find_seq_matches(sequencing_report,
                                                                          protein,
                                                                          specific_experiments)
    heatmap_absolute = np.zeros([heatmap_axis, heatmap_axis])


    for index_sample_one in range(heatmap_axis):
        d1 = unique_sequences.iloc[:, index_sample_one]
        shape_first_sample = sequencing_report[sequencing_report["Experiment"] == unique_experiments[index_sample_one]].shape[0]
        for index_sample_two in range(heatmap_axis):
            shape_second_sample = \
                sequencing_report[sequencing_report["Experiment"] == unique_experiments[index_sample_two]].shape[0]

            d2 = unique_sequences.iloc[:, index_sample_two]
            inner_join = np.where(np.logical_and(d1 == 1, d2 == 1))
            identity_both = inner_join[0].shape[0]
            lib_matches = d1 * d2
            lib_matches = lib_matches.sum()
            a = lib_matches
            b = d1.sum() - lib_matches
            c = d2.sum() - lib_matches
            heatmap_absolute[index_sample_one, index_sample_two] = float(lib_matches / ((shape_second_sample + shape_first_sample)/2))
            heatmap_absolute[index_sample_two, index_sample_one] = float(lib_matches / ((shape_second_sample + shape_first_sample)/2))
    matrix = pd.DataFrame(heatmap_absolute,
                          index=list(unique_experiments),
                          columns=list(unique_experiments))
    return matrix, unique_sequences, unique_experiments