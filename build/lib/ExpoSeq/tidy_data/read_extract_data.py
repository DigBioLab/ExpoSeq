import re
import numpy as np
import pandas as pd



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



