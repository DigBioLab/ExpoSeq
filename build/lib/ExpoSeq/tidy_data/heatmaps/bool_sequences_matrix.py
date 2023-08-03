import pandas as pd
def find_seq_matches(sequencing_report, protein, specific_experiments):
    if specific_experiments != False:
        sequencing_report = sequencing_report[sequencing_report['Experiment'].isin(specific_experiments)]
    else:
        pass
    unique_experiments = sequencing_report["Experiment"].unique()
    if protein == True:
        strand_column = "aaSeqCDR3"
    else:
        strand_column = "nSeqCDR3"
    unique_sequences = pd.DataFrame(sequencing_report[strand_column].unique())
    unique_sequences.rename(columns={0: strand_column},
                            inplace=True)

    for i, experiment in enumerate(unique_experiments):
        local_data = sequencing_report[sequencing_report["Experiment"] == experiment][strand_column].unique()
        local_data = pd.DataFrame(local_data, columns=["aaSeqCDR3"])
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