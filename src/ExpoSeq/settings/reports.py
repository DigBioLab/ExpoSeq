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
        avail_cols = [col for col in self.origin_seq_report.columns if col.startswith('nSeq') and self.origin_seq_report[col].count() > threshold_count]
        if "targetSequences" in self.origin_seq_report.columns.tolist():
            avail_cols += ["targetSequences"]
        else:
            pass
        return avail_cols

    
    def filter_region(self, region_string, remove_gaps = True):
        fixed_cols = ["Experiment", "cloneId", "readCount", "readFraction"]
        if region_string == "targetSequences":
            aa_string = "aaSeqtargetSequences"
            nseq_string = "nSeqtargetSequences"
            self.sequencing_report = self.origin_seq_report.copy()
            self.sequencing_report.rename(columns={region_string: nseq_string}, inplace=True) 
            if remove_gaps:
                self.sequencing_report = self.sequencing_report[self.sequencing_report['nSeq' + region_string].apply(self.is_divisible_by_three)]
            self.sequencing_report[aa_string] = self.sequencing_report[nseq_string].apply(self.translate_nucleotide_to_amino_acid)
            added_columns = ["nSeq" + region_string, "aaSeq" + region_string] #"minQual" + region_string,
            cols_of_interest = fixed_cols + added_columns
            self.sequencing_report = self.sequencing_report[cols_of_interest]
        else:
            added_columns = ["nSeq" + region_string, "aaSeq" + region_string] #"minQual" + region_string,
            cols_of_interest = fixed_cols + added_columns
            self.sequencing_report = self.origin_seq_report[cols_of_interest]
            if remove_gaps:
                self.sequencing_report = self.sequencing_report[self.sequencing_report['nSeq' + region_string].apply(self.is_divisible_by_three)] # removes gaps indirectly
        self.sequencing_report = self.remove_seq_errors(self.sequencing_report, region_string)

       # self.sequencing_report["aaSeq" + region_string] = self.sequencing_report["nSeq" + region_string].apply(self.translate_sequence)
    
    @staticmethod
    def translate_nucleotide_to_amino_acid(nucleotide_sequence):
        seq = Seq(nucleotide_sequence)
        return str(seq.translate())    

    def trim_data(self, region_string, length_threshold = 9, min_read_count = 0, new_fraction = "cloneFraction"):
        aa_string = "aaSeq" + region_string
        nseq_string = "nSeq" + region_string
        assert aa_string in self.sequencing_report.columns.tolist(), f"{self.sequencing_report.columns.tolist()}"
        assert nseq_string in self.sequencing_report.columns.tolist(), f"{self.sequencing_report.columns.tolist()}"
        length_threshold = length_threshold - 1
        if min_read_count > 0:
            min_read_count = min_read_count - 1
        
        sequencing_report = self.sequencing_report.groupby("Experiment", group_keys = True).apply(lambda group: group.drop_duplicates(subset=[aa_string], keep="first")).reset_index(drop=True)

        indexes_to_drop = sequencing_report[sequencing_report[aa_string] == 'region_not_covered'].index
        sequencing_report = sequencing_report.drop(indexes_to_drop)

        sequencing_report["lengthOfCDR3"] = sequencing_report[nseq_string].str.len() ## assumes all that pipeline starts with cdr3 region

        sequencing_report = sequencing_report.reset_index()

        sequencing_report = sequencing_report[(sequencing_report["aaSeq" + region_string].str.len()) > length_threshold ]
        sequencing_report = sequencing_report[(sequencing_report["readCount"] > (min_read_count + 1))]

        new_column = sequencing_report['readCount'] / sequencing_report.groupby('Experiment')['readCount'].transform('sum')
        self.sequencing_report = sequencing_report.copy()
        self.sequencing_report[new_fraction] = np.array(new_column)
        self.sequencing_report.drop(columns = ["readFraction"], inplace = True)
        self.sequencing_report.drop(columns = ["index"], inplace = True)
        
        
    def remove_not_covered(self):
        self.sequencing_report = self.sequencing_report[~self.sequencing_report.applymap(lambda x: x == "region_not_covered").any(axis=1)]

    
    @staticmethod
    def remove_seq_errors(sequencing_report, region_string):
        sequencing_report = sequencing_report.loc[~sequencing_report["aaSeq" + region_string].str.contains("[*]")]
        return sequencing_report
    
    def prepare_seq_report(self, region_string, length_threshold, min_read_count, remove_gaps = True):
        self.filter_region(region_string, remove_gaps) # sequencing errors are removed here and sequences with gaps are indirectly removed with: divisible_by = 3
        self.trim_data(region_string, length_threshold, min_read_count,)
        self.remove_not_covered()
       # self.remove_seq_errors()

    
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
        
    def renew_exp_names_origin(self, replacement_mapping, path_report):
        self.origin_seq_report["Experiment"] = self.origin_seq_report["Experiment"].replace(replacement_mapping)
        user_input = input("Do you want to save and change the new names in the sequencing report? - If so, you will not have to change the names again if you restart the program. Y/n")
        while True:
            if user_input.lower() in ["Y", "y", "n", "N"]:
                if user_input.lower() in ["Y", "y"]:
                    self.origin_seq_report.to_csv(os.path.join(path_report, "sequencing_report.csv"), index = False)
                else:
                    pass
                break
            else:
                print("Please enter Y or n")
            
        


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





