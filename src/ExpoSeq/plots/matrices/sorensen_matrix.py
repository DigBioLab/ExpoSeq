import numpy as np
import pandas as pd

class PrepareData:
    @staticmethod
    def find_seq_matches(sequencing_report,region_of_interest, protein, specific_experiments):
        if specific_experiments != False:
            sequencing_report = sequencing_report[sequencing_report['Experiment'].isin(specific_experiments)]
        else:
            pass
        unique_experiments = sequencing_report["Experiment"].unique()
        if protein == True:
            strand_column = region_of_interest
        else:
            strand_column =  region_of_interest
        unique_sequences = pd.DataFrame(sequencing_report[strand_column].unique())
        unique_sequences.rename(columns={0: strand_column},
                                inplace=True)

        for i, experiment in enumerate(unique_experiments):
            local_data = sequencing_report[sequencing_report["Experiment"] == experiment][strand_column].unique()
            local_data = pd.DataFrame(local_data, columns=[strand_column])
            unique_sequences = pd.merge(unique_sequences, local_data,
                                        on=strand_column,
                                        how='left',
                                        indicator=True)
            unique_sequences[experiment] = unique_sequences["_merge"].eq("both")
            unique_sequences = unique_sequences.drop("_merge", axis=1)
        unique_sequences.drop(strand_column,
                            inplace=True,
                            axis=1)
        unique_sequences = unique_sequences.astype("int")
        heatmap_axis = len(unique_experiments)
        return heatmap_axis, unique_sequences, unique_experiments
    
    
    def sorensen_matrix(self, sequencing_report, protein,region_of_interest, specific_experiments):
        heatmap_axis, unique_sequences, unique_experiments = self.find_seq_matches(sequencing_report,
                                                                            region_of_interest,
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
        return matrix, unique_experiments
    
    
