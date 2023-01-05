import numpy as np
import seaborn as sns
from python_scripts.tidy_data.heatmaps.tidy_morosita_horn import cleaning_data
#from python_scripts.main import local_pattern_more_digits, grouped_filenames
import pandas as pd

def cleaning_jaccard(sequencing_report, strand_column):
    unique_sequences = pd.DataFrame(sequencing_report.nSeqCDR3.unique())
    unique_sequences.rename(columns={0: strand_column},
                            inplace=True)
    for experiment in unique_experiments:
        local_data = sequencing_report[sequencing_report["Experiment"] == experiment][strand_column]
        unique_sequences = unique_sequences.merge(local_data, on=strand_column, how='left', indicator=True)
        unique_sequences[experiment] = unique_sequences.pop('_merge').eq("both")
    unique_sequences.drop(strand_column,
                        inplace=True,
                        axis=1)
    unique_sequences = unique_sequences.astype("int")
    heatmap_axis = len(unique_experiments)
    heatmap_absolute_jaccard = np.zeros([heatmap_axis, heatmap_axis])
    for index_sample_one in range(heatmap_axis):
        d1 = unique_sequences.iloc[:, index_sample_one]
        for index_sample_two in range(heatmap_axis):
            d2 = unique_sequences.iloc[:, index_sample_two]
            lib_matches = d1 * d2
            lib_matches = lib_matches.sum()
            a = lib_matches
            b = d1.sum() - lib_matches
            c = d2.sum() - lib_matches
            heatmap_absolute_jaccard[index_sample_one, index_sample_two] = lib_matches/(lib_matches+b+c)
            heatmap_absolute_jaccard[index_sample_two, index_sample_one] = lib_matches / (lib_matches + b + c)
    matrix = pd.DataFrame(heatmap_absolute_jaccard,
                          index = list(unique_experiments),
                          columns = list(unique_experiments))

    return matrix, unique_sequences, unique_experiments


def heatmap_sorensen(sequencing_report, protein):
    unique_experiments = sequencing_report["Experiment"].unique()
    if protein == True:
        strand_column = "aaSeqCDR3"
    else:
        strand_column = "nSeqCDR3"
    unique_sequences = pd.DataFrame(sequencing_report.nSeqCDR3.unique())
    unique_sequences.rename(columns={0: strand_column},
                            inplace=True)
    for experiment in unique_experiments:
        local_data = sequencing_report[sequencing_report["Experiment"] == experiment][strand_column]
        unique_sequences = unique_sequences.merge(local_data, on=strand_column, how='left', indicator=True)
        unique_sequences[experiment] = unique_sequences.pop('_merge').eq("both")
    unique_sequences.drop(strand_column,
                          inplace=True,
                          axis=1)
    unique_sequences = unique_sequences.astype("int")
    heatmap_axis = len(unique_experiments)
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


def heatmap_direct(sequencing_report, strand_column ):
    unique_sequences = pd.DataFrame(sequencing_report.nSeqCDR3.unique())
    unique_sequences.rename(columns={0: strand_column},
                            inplace=True)
    for experiment in unique_experiments:
        local_data = sequencing_report[sequencing_report["Experiment"] == experiment][strand_column]
        unique_sequences = unique_sequences.merge(local_data, on=strand_column, how='left', indicator=True)
        unique_sequences[experiment] = unique_sequences.pop('_merge').eq("both")

    unique_sequences.drop(strand_column,
                          inplace=True,
                          axis=1)

    unique_sequences = unique_sequences.astype("int")
    heatmap_axis = len(unique_experiments)
    heatmap_absolute_sorensen = np.zeros([heatmap_axis, heatmap_axis])

    for index_sample_one in range(heatmap_axis):
        d1 = unique_sequences.iloc[:, index_sample_one]
        start = 0
        shape_first_sample = \
        sequencing_report[sequencing_report["Experiment"] == unique_experiments[index_sample_one]].shape[0]
        for index_sample_two in range(heatmap_axis):
            shape_second_sample = \
            sequencing_report[sequencing_report["Experiment"] == unique_experiments[index_sample_two]].shape[0]
            d2 = unique_sequences.iloc[:, index_sample_two]
            lib_matches = d1 * d2
            lib_matches = lib_matches.sum()
            a = lib_matches
            b = d1.sum() - lib_matches
            c = d2.sum() - lib_matches
            heatmap_absolute[index_sample_one, index_sample_two] = float(lib_matches / shape_second_sample)
            heatmap_absolute[index_sample_two, index_sample_one] = float(lib_matches / shape_second_sample)
    matrix = pd.DataFrame(heatmap_absolute_sorensen,
                          index = list(unique_experiments),
                          columns = list(unique_experiments))
    return matrix, unique_sequences, unique_experiments


