import numpy as np
import pandas as pd
from .bool_sequences_matrix import find_seq_matches


def heatmap_sorensen(sequencing_report, protein, specific_experiments):
    heatmap_axis, unique_sequences, unique_experiments = find_seq_matches(sequencing_report,
                                                                          protein,
                                                                          specific_experiments)
    heatmap_absolute_sorensen = np.zeros([heatmap_axis, heatmap_axis])
    for index_sample_one in range(heatmap_axis):
        d1 = unique_sequences.iloc[:, index_sample_one]
        sequencing_report[sequencing_report["Experiment"] == unique_experiments[index_sample_one]].shape[0]
        for index_sample_two in range(heatmap_axis):
            sequencing_report[sequencing_report["Experiment"] == unique_experiments[index_sample_two]].shape[0]
            d2 = unique_sequences.iloc[:, index_sample_two]
            lib_matches = d1 * d2
            lib_matches = lib_matches.sum()
            b = d1.sum() - lib_matches
            c = d2.sum() - lib_matches
            heatmap_absolute_sorensen[index_sample_one, index_sample_two] = 2 * lib_matches / (2 * lib_matches + b + c)
            heatmap_absolute_sorensen[index_sample_two, index_sample_one] = 2 * lib_matches / (2 * lib_matches + b + c)
    matrix = pd.DataFrame(heatmap_absolute_sorensen,
                          index = list(unique_experiments),
                          columns = list(unique_experiments))
    return matrix, unique_sequences, unique_experiments