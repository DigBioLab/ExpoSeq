import pandas as pd
from .bool_sequences_matrix import find_seq_matches
import numpy as np

#sequencing_report = pd.read_table(r"C:\Users\nilsh\PycharmProjects\ExpoSeq\try_outs\test_data\sequencing_report.txt", sep = ",")
#sequencing_report = sequencing_report.drop('Unnamed: 0', axis=1)
#protein = True
def cleaning_jaccard(sequencing_report, protein, specific_experiments):
    heatmap_axis, unique_sequences, unique_experiments = find_seq_matches(sequencing_report, protein, specific_experiments)
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



  #  if rename_from_dic == True:
   #     matrix.rename(columns = Library_2_to_panning,
    #                  inplace = True)
     #   matrix.rename(index = Library_2_to_panning,
      #                inplace = True)