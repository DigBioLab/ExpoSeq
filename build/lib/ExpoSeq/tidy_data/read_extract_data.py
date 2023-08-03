import re
import numpy as np
import pandas as pd

# output == number of aligned reads
# input == number of total sequencing reads

# create class depending on which data should be analyzed

class read_intermediate_reports:
    def __init__(self, filename):
        self.sequencing_report = filename
        #questions have to be otuside of class, otherwise they'll be asked each time a file is read

    def filter_rows_on_min(self, column, min_count):
        self.sequencing_report = self.sequencing_report[(self.sequencing_report[column] > min_count)]
        return self.sequencing_report

    def filter_rows_not_divisible(self, divisor, column):
        self.sequencing_report = self.sequencing_report[(self.sequencing_report[column] % divisor) == 0]
        return self.sequencing_report

    def get_file_length(self):
        length = self.sequencing_report.shape[0]
        return length

    def get_individuals(self, column):
        unique_values = self.sequencing_report[column].unique()
        return unique_values
    def sort_table(self, sort_params):
        self.sequencing_report = self.sequencing_report.sort_values(by = sort_params)
        return self.sequencing_report

    def nest_data(self, nest_by, column_to_list, unique):
        if unique == True:
            apply_fun = lambda x: pd.unique([z for z in x]).tolist()
        else:
            apply_fun = lambda x: list(x)
        self.sequencing_report = self.sequencing_report.groupby(nest_by)[column_to_list].agg(apply_fun).reset_index()
        return self.sequencing_report
    def sum_nest_data(self, nested_column_name, axis):
        column = self.sequencing_report[nested_column_name]
        clone_counts = np.array(column)
        summed_counts = []
        for nest in clone_counts:
            sum = np.sum(nest, axis = axis)
            summed_counts.append(sum)
        summed_counts = np.array(summed_counts)
        return summed_counts

    def max_row(self, column):
        length_nested_seq = []
        for nested in self.sequencing_report[column]:
            length_nested_seq.append(len(nested))
            length_nested_seq = np.array(length_nested_seq)
        max = max(length_nested_seq)
        return max

    def extract_substring_rows(self, lib_name, col_name = "Experiment"):
        sub_table = self.sequencing_report[self.sequencing_report[col_name].str.contains(lib_name)]
        return sub_table

    def summarize_duplicates(self, column_to_sum, duplicate_column, group):
        if group != '' or False:
            group.append(duplicate_column)
            group_columns = group
        else:
            group_columns = [duplicate_column]
        column = self.sequencing_report.groupby(group_columns)[column_to_sum].transform('sum')
        self.sequencing_report[column_to_sum] = column
        self.sequencing_report = self.sequencing_report.drop_duplicates(group_columns,
                                                                        keep = 'last')
        return self.sequencing_report

    def map_to_column(self, object_to_map, column_to_map, new_column_name):
        self.sequencing_report[new_column_name] = self.sequencing_report[column_to_map].map(object_to_map)
        return self.sequencing_report






class read_alignment_report:
    def __init__(self, alignment_report):
        self.report = alignment_report
        self.digit = "\d+"

    def find_word_align_rep(self, word, regex): #"Total sequencing reads"  # auomate this or ask user how it is called in his / her file
        for line in self.report:
            if line.find(word) != -1:
                preferred_line = line
                total_sequenced_reads = re.search(regex,
                                                preferred_line).group()
        return total_sequenced_reads



