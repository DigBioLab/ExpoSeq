import pandas as pd
import re

# output == number of aligned reads
# input == number of total sequencing reads

# create class depending on which data should be analyzed
with open(file) as f:
    txtfile = f.readlines()
class read_intermediate_reports:
    def __init__(self, filename, alignment_report, input_pattern):
        self.sequencing_report = pd.read_table(filename)
        self.alignment_report = alignment_report
        self.digit = "\d+"
        self.pattern = input_pattern
        #questions have to be otuside of class, otherwise they'll be asked each time a file is read


    def extract_evalues_and_counts(self, ignore_zeros):
        if ignore_zeros == "Y" or "y":
            delete_zeros = sequencing_report['lengthOfCDR3'] % 3 == 0
            sequencing_report = sequencing_report[delete_zeros]
            evalues = sequencing_report.cloneCount.sum()

        if ignore_one_counts == "Y" or "y":
            clone_counts_tidied = sequencing_report[(sequencing_report['cloneCount'] > min_count)]
        number_clones = sequencing_report.shape[1]
        return evalues, clone_counts_tidied, number_clones

    def find_word_align_rep(self, word, regex): #"Total sequencing reads"  # auomate this or ask user how it is called in his / her file
        for line in alignment_report:
            if line.find(word) != -1:
                preferred_line = line
                total_sequenced_reads = re.search(regex,
                                                preferred_line).group()
                total_sequenced_reads = int(total_sequenced_reads)
        return total_sequenced_reads






