from python_scripts.tidy_data.read_extract_data import read_intermediate_reports
import numpy as np
from python_scripts.genetic_dogma import genetic_dogma
from python_scripts.tidy_data.interpret_data import mapFunc
import pandas as pd
from python_scripts.test_data.rename_labels import Library_2_to_panning


def cleaning_data(protein, sequencing_report_input, tidy_data, specific_experiments = False,
                  experiments = False, new_fraction = 'new_fraction'):
    if specific_experiments != False:
        sequencing_report_input = sequencing_report_input[sequencing_report_input['Experiment'].isin(experiments)]
    else:
        pass

    report = tidy_data.sequencing_report
    new_column = report['cloneCount'] / report.groupby('Experiment')['cloneCount'].transform('sum')
    report[new_fraction] = np.array(new_column)
    tidy_data.sequencing_report = report
    if protein == True:
        tidy_data.sequencing_report = mapFunc(sequencing_report = tidy_data.sequencing_report,
                                                column = 'nSeqCDR3',
                                                func = genetic_dogma,
                                                column_name = 'peptide_seq')
        strang_column = 'peptide_seq'
    else:
        strang_column = 'nSeqCDR3'
    tidy_data.summarize_duplicates(column_to_sum = new_fraction,
                                   duplicate_column = strang_column,
                                    group = ['Experiment']) # group has to be list
    sequencing_report = tidy_data.sequencing_report
    unique_experiments = tidy_data.get_individuals("Experiment")
    unique_sequences = pd.DataFrame(sequencing_report.nSeqCDR3.unique()) #nseqcdr3 has to be changed since it is a chosen name for column
    unique_sequences.rename(columns={0: strang_column},
                            inplace=True)
    unique_experiments = list(unique_experiments)
    for experiment in unique_experiments:
        local_data = sequencing_report[sequencing_report["Experiment"] == experiment][[strang_column, new_fraction]]
        unique_sequences = unique_sequences.merge(local_data,
                                                  how='left',
                                                  on = strang_column)
        unique_sequences = unique_sequences.rename(columns={new_fraction: new_fraction + experiment})
    unique_sequences.drop(strang_column,
                          inplace=True,
                          axis=1)
    unique_sequences = unique_sequences.fillna(0)
    return unique_sequences, unique_experiments



def morosita_horn_matrix(unique_sequences, unique_experiments):
    columns = len(unique_experiments)
    morosita_horn_matrix = np.zeros((columns, columns))

    for index_sample_one in range(columns):
        d1 = unique_sequences.iloc[:, index_sample_one]
        start = 0
        for index_sample_two in range(columns):
            d2 = unique_sequences.iloc[:, index_sample_two]
            product_d1_d2 = d1 * d2
            d1_square = d1**2
            d2_square = d2**2
            cal_cols=pd.concat([d1,d2,product_d1_d2,d1_square,d2_square], axis=1).fillna(0)
            cal_cols.columns=['d1','d2','product_d1_d2','d1_square','d2_square']
            sub = cal_cols['product_d1_d2']>0
            subset = cal_cols[sub]
            #MH = (2*product_d1_d2.sum()) /(d1_square.sum() + d2_square.sum())
            MH = 2*subset['product_d1_d2'].sum()/(subset['d1_square'].sum()+subset['d2_square'].sum())
            morosita_horn_matrix[index_sample_one, index_sample_two] = MH
            morosita_horn_matrix[index_sample_two, index_sample_one] = MH
            start += 1
    matrix = pd.DataFrame(morosita_horn_matrix,
                          index = list(unique_experiments),
                          columns = list(unique_experiments))
    matrix.rename(columns = Library_2_to_panning,
                  inplace = True)
    matrix.rename(index = Library_2_to_panning,
                  inplace = True)
    return matrix


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
    matrix.rename(columns = Library_2_to_panning,
                  inplace = True)
    matrix.rename(index = Library_2_to_panning,
                  inplace = True)
    return morosita_horn_matrix, matrix
