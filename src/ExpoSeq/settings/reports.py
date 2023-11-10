import pickle
from Bio.Seq import Seq
import os
import pandas as pd
import numpy as np
from ..augment_data.binding_data import collect_binding_data
from Bio.Seq import Seq



class SequencingReport:
    def __init__(self,sequencing_report):
        sequencing_report = sequencing_report.dropna(subset = ["Experiment"])
        self.origin_seq_report = sequencing_report.copy()
        self.sequencing_report = sequencing_report

    def get_exp_names(self, experiment_path):
        try:
            with open(experiment_path, "rb") as f:
                unique_experiments = pickle.load(f)
        except:
            unique_experiments = self.sequencing_report["Experiment"].unique().tolist()
            unique_experiments = dict(zip(unique_experiments, unique_experiments))
        return unique_experiments
            
    def map_exp_names(self, unique_experiments):
        self.sequencing_report["Experiment"] = self.sequencing_report["Experiment"].map(unique_experiments)


    def is_divisible_by_three(self, seq):
        return len(seq) % 3 == 0

    @staticmethod
    def translate_seq(nucleotide):
        return Seq(nucleotide).translate()

    def get_fragment(self):
        threshold_count = 0.1 * len(self.origin_seq_report) # to get sequences which have values
        return [col for col in self.origin_seq_report.columns if col.startswith('nSeq') and self.origin_seq_report[col].count() > threshold_count]

    
    def filter_region(self, region_string):
        added_columns = ["nSeq" + region_string, "aaSeq" + region_string] #"minQual" + region_string,
        fixed_cols = ["Experiment", "cloneId", "readCount", "readFraction"]
        cols_of_interest = fixed_cols + added_columns
        self.sequencing_report = self.origin_seq_report[cols_of_interest]
        self.sequencing_report = self.sequencing_report[self.sequencing_report['nSeq' + region_string].apply(self.is_divisible_by_three)]
       # self.sequencing_report["aaSeq" + region_string] = self.sequencing_report["nSeq" + region_string].apply(self.translate_sequence)
    
        

    def trim_data(self, region_string, divisible_by = 3,length_threshold = 9, min_read_count = 0, new_fraction = "cloneFraction"):
        length_threshold = length_threshold - 1
        if min_read_count > 0:
            min_read_count = min_read_count - 1
        sequencing_report = self.sequencing_report.groupby("Experiment", group_keys = True).apply(lambda group: group.drop_duplicates(subset=["nSeq" + region_string], keep="first")).reset_index(drop=True)

        indexes_to_drop = sequencing_report[sequencing_report["aaSeq" + region_string] == 'region_not_covered'].index
        sequencing_report = sequencing_report.drop(indexes_to_drop)

        sequencing_report["lengthOfCDR3"] = sequencing_report["nSeq" + region_string].str.len() ## assumes all that pipeline starts with cdr3 region

        sequencing_report = sequencing_report.reset_index()

        sequencing_report = sequencing_report[(sequencing_report["lengthOfCDR3"] % divisible_by) == 0]
        sequencing_report = sequencing_report[(sequencing_report["aaSeq" + region_string].str.len()) > length_threshold ]
        sequencing_report = sequencing_report[(sequencing_report["readCount"] > (min_read_count + 1))]

        new_column = sequencing_report['readCount'] / sequencing_report.groupby('Experiment')['readCount'].transform('sum')
        self.sequencing_report = sequencing_report.copy()
        self.sequencing_report[new_fraction] = np.array(new_column)
        self.sequencing_report.drop(columns = ["readFraction"], inplace = True)
        self.sequencing_report.drop(columns = ["index"], inplace = True)

    def remove_not_covered(self):
        self.sequencing_report = self.sequencing_report[~self.sequencing_report.applymap(lambda x: x == "region_not_covered").any(axis=1)]

    def prepare_seq_report(self, region_string, divisible_by, length_threshold, min_read_count):

        self.filter_region(region_string)
        self.trim_data(region_string, divisible_by, length_threshold, min_read_count,)
        self.remove_not_covered()
        
    def find_longest_sequence(self):
        columns = self.origin_seq_report.columns.to_list()
        regions_order = ['nSeqFR1', 'nSeqCDR1', 'nSeqFR2', 'nSeqCDR2', 'nSeqFR3', 'nSeqCDR3', 'nSeqFR4']
        # Map regions to their indices in the order
        indices = [regions_order.index(col) for col in columns if col in regions_order]
        # If no regions are present, return an empty list
        if not indices:
            return []
        # Find the longest consecutive subsequence
        longest_seq = []
        current_seq = [indices[0]]
        for i in range(1, len(indices)):
            if indices[i] - current_seq[-1] == 1:  # Check for consecutiveness
                current_seq.append(indices[i])
            else:
                # If not consecutive, update longest sequence if needed and reset current sequence
                if len(current_seq) > len(longest_seq):
                    longest_seq = current_seq
                current_seq = [indices[i]]
        # Handle the case where the longest sequence is at the end
        if len(current_seq) > len(longest_seq):
            longest_seq = current_seq
        # Convert indices back to region names
        return [regions_order[idx] for idx in longest_seq]

        
    def filter_longest_sequence(self, divisible_by= 3, length_threshold= 9, min_read_count = 0):
        cols_longest = self.find_longest_sequence()
        region_string = cols_longest[0]
        for i in cols_longest:
            i = i.replace("nSeq", "")
            self.trim_data(i, divisible_by= divisible_by, length_threshold= length_threshold, min_read_count = min_read_count)
        cols_of_interest = ["Experiment", "cloneId", "readCount", "cloneFraction"] + cols_longest
        found_frags = self.sequencing_report[cols_of_interest]
        found_frags = found_frags[~found_frags.applymap(lambda x: x == "region_not_covered").any(axis=1)]
        found_frags['merged_with_spaces'] = found_frags[cols_longest].apply(lambda row: ' '.join(row), axis=1)
        found_frags["nSeqall"] = found_frags['merged_with_spaces'].str.replace('\\', '').str.replace(' ', '')
        found_frags['aaSeqall'] = found_frags['nSeqall'].apply(self.translate_seq)
        found_frags.reset_index(inplace = True)
        region_string = "all"
        self.sequencing_report = found_frags
        return region_string

        
    def check_sample_name(self, module_dir, experiment_name):
        sample_names = self.origin_seq_report["Experiment"].unique().tolist()
        for sample in sample_names:
            if len(sample) > 22:
                new_sample_name = input(f"Your sample name for {sample} is to long. Please enter a new name with maximum 22 characters.")
                replacement_mapping = {sample: new_sample_name}
                self.origin_seq_report['Experiment'] = self.origin_seq_report['Experiment'].replace(replacement_mapping)
                self.sequencing_report["Experiment"] = self.sequencing_report["Experiment"].replace(replacement_mapping)
        path_file = os.path.join(module_dir, "my_experiments", experiment_name, "sequencing_report.csv")
        if os.path.isdir(os.path.dirname(path_file)):
            self.origin_seq_report.to_csv(path_file, index = False)
        
        
class Test_SequencingReport:
    def __init__():
        divisble_by = 3
        length_threshold = 3
        min_read_count = 2
        example_report = os.path.join("sequencing_report.csv")
        Report = SequencingReport(example_report)
        Report.prepare_seq_report(region_string = "CDR3", divisible_by=divisble_by, length_treshold = length_threshold, min_read_count=min_read_count)
        test = Report.sequencing_report
        assert test["aaSeqCDR3"].str.len() > length_threshold, "Sequence length filter does not work"
        assert test["readCount"] > min_read_count, "Read count filter does not work"
        assert test["nSeqCDR3"].str.len() % 3 == 0, "Filter sequences by divisor does not work"
        unitest_dir = os.path.join(os.getcwd(), "my_experiments", "unitest")
        if not os.path.isdir(unitest_dir):
            os.mkdir(os.path.join(os.getcwd(), "my_experiments", "unitest"))
        Report.check_sample_name(os.getcwd(), "not_existent")
        assert not os.path.isfile(os.path.join(os.getcwd(), "my_experiments", "not_existend", "sequencing_report.csv"))
        
  

class BindingReport:
    def __init__(self, module_dir, experiment):
        self.binding_data_dir = os.path.join(module_dir,
                                "my_experiments",
                                experiment,
                                "binding_data.csv")
        
    def ask_binding_data(self):
        if not os.path.isfile(self.binding_data_dir):
            add_binding = input("Do you have binding Data? Y/n")
            if add_binding.lower() in ["Y", "y"]:
                binding_data = collect_binding_data()

                binding_data.to_csv("binding_data.csv")
            else:
                binding_data = None
        else:
            binding_data = pd.read_csv(self.binding_data_dir)
        if binding_data is not None:
            binding_data.drop_duplicates(subset = "aaSeqCDR3", inplace = True)
        return binding_data





