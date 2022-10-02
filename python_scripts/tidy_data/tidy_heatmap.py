from python_scripts.tidy_data.read_extract_data import read_intermediate_reports
import numpy as np
from python_scripts.genetic_dogma import genetic_dogma
from python_scripts.tidy_data.interpret_data import mapFunc

class heatmap_creator(read_intermediate_reports):
    def __init__(self, filename, input_pattern, divisible_by=3, min_count=1):
        super().__init__(filename, input_pattern)
        self.filter_rows_not_divisible(divisor=divisible_by,
                                        column='lengthOfCDR3')
        self.filter_rows_on_min(column = 'cloneCount',
                                min_count = min_count)
        self.summarize_duplicates(column_to_sum = 'cloneFraction',
                                  duplicate_column = ['nSeqCDR3', 'Experiment'])

    def tidy_for_protein(self,):
        self.sequencing_report = mapFunc(self.sequencing_report,
                                         'nSeqCDR3',
                                         genetic_dogma,
                                         'peptide_seq')
        self.nest_data(nest_by="Experiment",
                        column_to_list=['peptide_seq', 'cloneFraction'],
                        unique=False)
        sequencing_report_all = self.sequencing_report
        unique_experiments = self.get_individuals("Experiment")
        max_experiments = unique_experiments.shape[0]
        return sequencing_report_all, max_experiments

    def tidy_for_dna(self):

        self.nest_data(nest_by="Experiment",
                        column_to_list=['nSeqCDR3', 'cloneFraction'],
                        unique=False)
        sequencing_report_all = self.sequencing_report
        unique_experiments = self.get_individuals("Experiment")
        max_experiments = unique_experiments.shape[0]
        return sequencing_report_all, max_experiments

def cleaning_dna(local_pattern_more_digits, sequencing_report_all, min_count = 1, divisible_by = 3):

    tidy_data = read_intermediate_reports(filename = sequencing_report_all,
                                            input_pattern = local_pattern_more_digits)
    tidy_data.filter_rows_not_divisible(divisor = divisible_by,
                                        column = 'lengthOfCDR3')
# if app insert a test and error report if someone fails to insert correct column name
    tidy_data.filter_rows_on_min(column = 'cloneCount',
                                min_count = min_count)
    tidy_data.summarize_duplicates(column_to_sum = 'cloneFraction',
                                   duplicate_column = 'nSeqCDR3')
    tidy_data.nest_data(nest_by = "Experiment",
                        column_to_list = ['nSeqCDR3', 'cloneFraction'],
                        unique = False)
    sequencing_report_all = tidy_data.sequencing_report
    unique_experiments = tidy_data.get_individuals("Experiment")
    max_experiments = unique_experiments.shape[0]

    return sequencing_report_all, max_experiments

def extract_data_heatmap(sequencing_report_all):
    nested_seq_pandas = sequencing_report_all.nSeqCDR3
    nested_seq = np.array(nested_seq_pandas)
    nested_frac_pandas = sequencing_report_all.cloneFraction
    nested_clone_fraction = np.array(nested_frac_pandas)
    return nested_seq, nested_clone_fraction