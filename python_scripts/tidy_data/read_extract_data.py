import pandas as pd
import re

# output == number of aligned reads
# input == number of total sequencing reads

# create class depending on which data should be analyzed

class read_intermediate_reports:
    def __init__(self, filename, input_pattern):
        self.sequencing_report = filename
        self.pattern = input_pattern
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



