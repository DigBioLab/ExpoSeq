import numpy as np
import pandas as pd


def cleaning_data(sequencing_report, protein = True, specific_experiments = False):
 #   new_fraction = 'new_fraction'
    if specific_experiments != False:
        sequencing_report = sequencing_report[sequencing_report['Experiment'].isin(specific_experiments)]
    else:
        pass
   # new_column = sequencing_report['readCount'] / sequencing_report.groupby('Experiment')['readCount'].transform('sum')
    #sequencing_report[new_fraction] = np.array(new_column)
    if protein == True:
   #     sequencing_report = mapFunc(sequencing_report = sequencing_report,
    #                                column = 'nSeqCDR3',
     #                               func = genetic_dogma,
       #                             column_name = 'peptide_seq')
        strand_column = 'aaSeqCDR3'
    else:
        strand_column = 'nSeqCDR3'

    group_columns = ["Experiment", strand_column]
   # column = sequencing_report.groupby(group_columns)[new_fraction].transform('sum')
    #sequencing_report[new_fraction] = column
    sequencing_report = sequencing_report.drop_duplicates(group_columns,
                                                            keep='last')
    unique_experiments = sequencing_report["Experiment"].unique()
    unique_sequences = pd.DataFrame(sequencing_report[strand_column].unique()) #nseqcdr3 has to be changed since it is a chosen name for column
    unique_sequences.rename(columns={0: strand_column},
                            inplace=True)
    unique_experiments = list(unique_experiments)
    for experiment in unique_experiments:
        local_data = sequencing_report[sequencing_report["Experiment"] == experiment][[strand_column,
                                                                                       "clonesFraction"]]
        unique_sequences = unique_sequences.merge(local_data,
                                                  how='left',
                                                  on = strand_column)
        unique_sequences = unique_sequences.fillna(0)
        unique_sequences = unique_sequences.rename(columns={"clonesFraction": "clonesFraction" + experiment})
    unique_sequences.drop(strand_column,
                          inplace=True,
                          axis=1)
    unique_sequences = unique_sequences.fillna(0)
    return unique_sequences, unique_experiments



def morosita_horn_matrix(unique_sequences, unique_experiments):

    columns = len(unique_experiments)
    morosita_horn_matrix = np.zeros((columns,
                                     columns))
    for index_sample_one in range(columns):
        d1 = unique_sequences.iloc[:, index_sample_one]
        start = 0
        for index_sample_two in range(columns):
            d2 = unique_sequences.iloc[:, index_sample_two]
            product_d1_d2 = d1 * d2
            d1_square = d1**2
            d2_square = d2**2
       #     cal_cols=pd.concat([d1,d2,product_d1_d2,d1_square,d2_square], axis=1).fillna(0)
        #    cal_cols.columns=['d1',
         #                     'd2',
          #                    'product_d1_d2',
           #                   'd1_square',
            #                  'd2_square']
            #sub = cal_cols['product_d1_d2'] > 0
            #subset = cal_cols[sub]
            MH = (2*np.sum(product_d1_d2)) /(np.sum(d1_square) + np.sum(d2_square))
            morosita_horn_matrix[index_sample_one, index_sample_two] = MH
            morosita_horn_matrix[index_sample_two, index_sample_one] = MH
          #  MH = 2*subset['product_d1_d2'].sum()/(subset['d1_square'].sum()+subset['d2_square'].sum())
           # morosita_horn_matrix[index_sample_one, index_sample_two] = MH
            #morosita_horn_matrix[index_sample_two, index_sample_one] = MH
            start += 1
    matrix = pd.DataFrame(morosita_horn_matrix,
                          index = list(unique_experiments),
                          columns = list(unique_experiments))
    return matrix, unique_sequences, unique_experiments


def morosita_horn_matrix_two(unique_sequences, unique_experiments):
    columns = len(unique_experiments)
    morosita_horn_matrix = np.zeros((columns, columns))
    column_names = list(unique_experiments)
    unique_sequences = np.array(unique_sequences)
    for index_sample_one in range(columns):
        d1 = unique_sequences[index_sample_one, :]
        start = 0
        for index_sample_two in range(columns):
            d2 = unique_sequences[index_sample_two, :]
            product_d1_d2 = d1 * d2
            d1_square = d1**2
            d2_square = d2**2
            MH = (2*np.sum(product_d1_d2)) /(np.sum(d1_square) + np.sum(d2_square))
            morosita_horn_matrix[index_sample_one, index_sample_two] = MH
            morosita_horn_matrix[index_sample_two, index_sample_one] = MH
            start += 1
    matrix = pd.DataFrame(morosita_horn_matrix,
                          index = column_names,
                          columns = column_names)

    return morosita_horn_matrix, matrix
