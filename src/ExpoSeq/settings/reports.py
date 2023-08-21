import pickle
from Bio.Seq import Seq
import os
import pandas as pd
import numpy as np
from ExpoSeq.augment_data.binding_data import collect_binding_data
from Bio.Seq import Seq



class SequencingReport:
    def __init__(self,sequencing_report):
        self.origin_seq_report = sequencing_report.copy()
        self.sequencing_report = sequencing_report.copy()

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
        fixed_cols = ["Experiment", "cloneId", "readCount", "cloneFraction"]
        cols_of_interest = fixed_cols + added_columns
        self.sequencing_report = self.sequencing_report[cols_of_interest]
        self.sequencing_report = self.sequencing_report[self.sequencing_report['nSeq' + region_string].apply(self.is_divisible_by_three)]
       # self.sequencing_report["aaSeq" + region_string] = self.sequencing_report["nSeq" + region_string].apply(self.translate_sequence)
    
        

    def trim_data(self, region_string, divisible_by = 3,length_threshold = 9, min_read_count = 0, new_fraction = "cloneFraction"):
        #new_fractions = clones_sample.groupby("nSeqCDR3")["readFraction"].sum().reset_index()
        self.sequencing_report = self.origin_seq_report.drop_duplicates(subset=["nSeq" + region_string], keep="first")
        indexes_to_drop = self.sequencing_report[self.sequencing_report["aaSeq" + region_string] == 'region_not_covered'].index
        self.sequencing_report = self.sequencing_report.drop(indexes_to_drop)
        # clones_sample = new_fractions.merge(clones_sample,
        #                                    how="left",
        #                                   on="nSeqCDR3")
        self.sequencing_report["lengthOfCDR3"] = self.sequencing_report["nSeqCDR3"].str.len() ## assumes all that pipeline starts with cdr3 region
        self.sequencing_report = self.sequencing_report.sort_values(by="cloneId")
        self.sequencing_report = self.sequencing_report.reset_index()
        #self.sequencing_report = self.sequencing_report.drop(columns=["readFraction_y", "index"])
        #self.sequencing_report = self.sequencing_report.rename(columns={"readFraction_x": "cloneFraction"})
        self.sequencing_report = self.sequencing_report[(self.sequencing_report["lengthOfCDR3"] % divisible_by) == 0]
        self.sequencing_report = self.sequencing_report[(self.sequencing_report["aaSeq" + region_string].str.len()) > length_threshold ]
        self.sequencing_report = self.sequencing_report[(self.sequencing_report["readCount"] > min_read_count)]
        new_column = self.sequencing_report['readCount'] / self.sequencing_report.groupby('Experiment')['readCount'].transform('sum')
        self.sequencing_report[new_fraction] = np.array(new_column)
        
    def remove_not_covered(self):
        self.sequencing_report = self.sequencing_report[~self.sequencing_report.applymap(lambda x: x == "region_not_covered").any(axis=1)]
    
    def prepare_seq_report(self, region_string, divisible_by, length_threshold, min_read_count):
        self.trim_data(region_string, divisible_by, length_threshold, min_read_count,)
        self.filter_region(region_string)
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
        return binding_data





