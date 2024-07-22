import pickle
from Bio.Seq import Seq
import os
import pandas as pd
import numpy as np
try:
    from tkinter import filedialog
except:
    pass
import warnings
import glob          
from .full_sequence_finder import FullSequence


class ManageImportFiles:
    def __init__(self):
        """This class can be used to create the sequencing report from a directory containing tsv files. It was primarily designed for the implementation in Platforma.bio. 
            Use the method merge_tsvs and input the directory
        
        """
        self.columns_not_wanted = ['refPoints', 'allVHitsWithScore', 'allDHitsWithScore',
                        'allJHitsWithScore', 'allCHitsWithScore', 'allVAlignments',
                        'allDAlignments', 'allJAlignments', 'allCAlignments']
        
    @staticmethod    
    def check_cols(single_table: pd.DataFrame, filename:str):
        cols = single_table.columns.to_list()
        nuc_regions = ["nSeqCDR1", "nSeqFR2", "nSeqCDR2", "nSeqFR3", "nSeqCDR3", "nSeqFR4"]
        aa_regions = ["aaSeqFR1", "aaSeqCDR1", "aaSeqFR2", "aaSeqCDR2","aaSeqFR3", "aaSeqCDR3", "aaSeqFR4"]
        mandatory_cols = ["readFraction", "readCount", "cloneId", "targetSequences"]
        for col in mandatory_cols:
            assert col in cols, f"{col} does not exist in {filename}"
        

    @staticmethod
    def get_filename(tsv_filename:str):
        whole_basename:str = os.path.basename(tsv_filename)
        experiment_name = whole_basename.split(".tsv")[0]
        return experiment_name
    
    def merge_tsvs(self, dir_tsvs:str):
        assert os.path.isdir(dir_tsvs), f"{dir_tsvs} does not exist"
        tsv_files = glob.glob(os.path.join(dir_tsvs, "*.tsv*"))
        sequencing_report = pd.DataFrame([])
        for tsv_table in tsv_files:
            experiment_name = self.get_filename(tsv_table)
            sample = pd.read_table(tsv_table)
            self.check_cols(sample, tsv_table)
            sample.drop(columns=self.columns_not_wanted, inplace=True)
            sample["Experiment"] = experiment_name #you dont need to worry about the length of the name because in Platforma they will decide how the headers look like
            sample["cloneFraction"] = sample["readFraction"] 
            sequencing_report = pd.concat([sequencing_report, sample])
        
        return sequencing_report

#sequencing_report = ManageImportFiles().merge_tsvs(r"C:\Users\nilsh\OneDrive\Desktop\results_thesis\data\ngs_data")
#sequencing_report.to_csv("sequencing_report.csv", index = False)


class SequencingReport:
    def __init__(self,sequencing_report:pd.DataFrame):
        """This class can be used to tidy and prepare the sequencing report for the pipeline. Alternatively you can input a string as a path to a directory with tsv files to generate a merged object out of these table which can be prepared as sequencing report then.

        Args:
            sequencing_report: _description_

        Returns:
            pd.DataFrame: tidied sequencing report
        """
        if isinstance(sequencing_report, pd.DataFrame) == True: # for pipeline purposes: 
            sequencing_report = sequencing_report.dropna(subset = ["Experiment"])
        elif type(sequencing_report) == str: # fur Platforma purposes to input tsv files
            if os.path.isdir(sequencing_report):
                TSVManager = ManageImportFiles()
                sequencing_report = TSVManager.merge_tsvs(sequencing_report) # contains test for path checking
            elif os.path.isfile(sequencing_report) and sequencing_report.endswith(".csv"):
                sequencing_report = pd.read_csv(sequencing_report)
            else:
                raise ValueError("Please enter a valid path to a directory with tsv files or a csv file")
        else:
            return ValueError
        
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


    def get_longest_region(self, aa_string:str, nseq_string:str):
        """This code modifies the sequencing_report attribute of the class. It will find the longest connecting regions and merges them together to safe them in one column 

        Args:
            aa_string (str): _description_
            nseq_string (str): _description_
        """
        self.sequencing_report = self.origin_seq_report.copy()
        columns_to_check = ["nSeqFR1", "nSeqCDR1", "nSeqFR2", "nSeqCDR2", "nSeqFR3", "nSeqCDR3", "nSeqFR4"]
        present_columns = [col for col in columns_to_check if col in self.sequencing_report]
        present_columns = [s.replace("nSeq", '') for s in present_columns]
        GetFullSeq = FullSequence(present_columns)
        longest_seq_regions = GetFullSeq.find_connecting_seq()
        aaSeq_cols = ["aaSeq" + region for region in longest_seq_regions]
        nSeq_cols = ["nSeq" + region for region in longest_seq_regions]
        mask = self.sequencing_report[nSeq_cols].apply(lambda row: 'region_not_covered' not in row.values, axis=1)
        self.sequencing_report = self.sequencing_report[mask]
        self.sequencing_report[aa_string] = self.sequencing_report[aaSeq_cols].apply(lambda row: ''.join(row.values.astype(str)), axis=1)
        self.sequencing_report[nseq_string] = self.sequencing_report[nSeq_cols].apply(lambda row: ''.join(row.values.astype(str)), axis=1)
        

    
    def filter_region(self, region_string, remove_gaps = True, remove_errors = True):
        fixed_cols = ["Experiment", "cloneId", "readCount", "readFraction"]
        if region_string == "targetSequences":
            aa_string = "aaSeqtargetSequences"
            nseq_string = "nSeqtargetSequences"
            self.get_longest_region(aa_string, nseq_string)
           # self.sequencing_report.rename(columns={region_string: nseq_string}, inplace=True) 
            if remove_gaps:
                self.sequencing_report = self.remove_gaps_in_peptide(self.sequencing_report, 'aaSeq' + region_string) 
            else:
                pass
         #   self.sequencing_report[aa_string] = self.sequencing_report[nseq_string].apply(self.translate_nucleotide_to_amino_acid) 
            
            added_columns = ["nSeq" + region_string, "aaSeq" + region_string] #"minQual" + region_string,
            cols_of_interest = fixed_cols + added_columns
            self.sequencing_report = self.sequencing_report[cols_of_interest]
            self.remove_not_covered()
        else:
            added_columns = ["nSeq" + region_string, "aaSeq" + region_string] #"minQual" + region_string,
            cols_of_interest = fixed_cols + added_columns
            
            self.sequencing_report = self.origin_seq_report[cols_of_interest]
            self.drop_nas()
            self.remove_not_covered()
            if remove_gaps: # removes all sequences with gaps in the peptide sequence for the corresponding region 
                self.sequencing_report = self.remove_gaps_in_peptide(self.sequencing_report, 'aaSeq' + region_string) 
        
        if remove_errors:
            self.sequencing_report = self.remove_seq_errors(self.sequencing_report, region_string)

       # self.sequencing_report["aaSeq" + region_string] = self.sequencing_report["nSeq" + region_string].apply(self.translate_sequence)
    
    @staticmethod
    def translate_nucleotide_to_amino_acid(nucleotide_sequence):
        seq = Seq(nucleotide_sequence)
        return str(seq.translate())    

    def trim_data(self, region_string, length_threshold = 9, min_read_count = 0, new_fraction = "cloneFraction"):
        aa_string = "aaSeq" + region_string
        nseq_string = "nSeq" + region_string

        length_threshold = length_threshold - 1
        if min_read_count > 0:
            min_read_count = min_read_count - 1
        
        sequencing_report = self.sequencing_report.groupby("Experiment", group_keys = True).apply(lambda group: group.drop_duplicates(subset=[aa_string], keep="first")).reset_index(drop=True)

        indexes_to_drop = sequencing_report[sequencing_report[aa_string] == 'region_not_covered'].index
        sequencing_report = sequencing_report.drop(indexes_to_drop)

        
        sequencing_report[f"lengthOf{region_string}"] = sequencing_report[nseq_string].str.len() ## assumes all that pipeline starts with cdr3 region

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
    def remove_gaps_in_peptide(sequencing_report, region_string):
        filter = sequencing_report[ region_string].str.contains("[_]", na=False)
        sequencing_report = sequencing_report.loc[~filter.values]
        return sequencing_report
    
    def drop_nas(self, ):
        self.sequencing_report = self.sequencing_report.dropna()
        
        
    @staticmethod
    def remove_seq_errors(sequencing_report, region_string):
        sequencing_report = sequencing_report.loc[~sequencing_report["aaSeq" + region_string].str.contains("[*]", na = False) ]
        return sequencing_report
    
    def _get_origin_order(self):
        """I could not find the code where the order of the sequencing report object is rearranged, so I created this function which sets the old order again, which is important for some plots.
        """
        order_mapping = {name: i for i, name in enumerate (self.origin_seq_report["Experiment"].unique())}
        self.sequencing_report["order"] = self.sequencing_report["Experiment"].map(order_mapping)
        self.sequencing_report = self.sequencing_report.sort_values("order").drop(columns = "order")
    
    def prepare_seq_report(self, region_string, length_threshold, min_read_count, remove_gaps = True, remove_errors = True):
        self.filter_region(region_string, remove_gaps, remove_errors) # sequencing errors are removed here and sequences with gaps are indirectly removed with: divisible_by = 3
        self.trim_data(region_string, length_threshold, min_read_count,)
        self._get_origin_order()
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
        
    @staticmethod
    def collect_binding_data(binding_data = None):
        if binding_data is None:
            binding_data = pd.DataFrame([])
        else:
            pass
        second_prompt = False
        while True:
            # prompt the user to add a file
            if not second_prompt:
                print("You can either add an excel sheet or a csv file which contains the binding data.\nNote: the first column must contain the CDR3 sequences and its column name has to be aaSeqCDR3.")
            else:
                print("Please modify the file and choose it again.")
            try:
                binding_file = filedialog.askopenfilename()
            except:
                while True:
                    binding_file = input("copy and paste the path to your binding report")
                    if os.path.isfile(os.path.abspath(binding_file)):
                        break
                    else:
                        print("Please enter a valid filepath. ")
            
            if binding_file.endswith(".xlsx") or binding_file.endswith(".csv") or binding_file.endswith(".tsv"):
                if binding_file.endswith(".xlsx"):
                    binding_new = pd.read_excel(binding_file)
                elif binding_file.endswith(".csv"):
                    binding_new = pd.read_csv(binding_file)
                elif binding_file.endswith(".tsv"):
                    binding_new = pd.read_table(binding_file)
                if (binding_new.index == binding_new.iloc[:, 0]).all():
                    binding_new = binding_new.drop(binding_new.columns[0], axis=1) # checks if first column is basically the index column
                else: pass

                if binding_new.columns.to_list()[0] == "aaSeqCDR3":
                    second_prompt = False
                    pass
                elif binding_new.columns.to_list()[0].capitalize() == "Sequences":
                    second_prompt = False
                    pass
                else:
                    second_prompt = True
                    print("Please chechk if your file follows the following requirements:\n1. The first column must contain the CDR3 sequences and its column name has to be aaSeqCDR3.\n2. Check that you separate your data by comma or tab.\n3. Your file needs to end with .xlsx, .csv or .tsv")
            else:
                print("Please enter a valid filepath to a csv or xlsx file")
                second_prompt = True

            if not second_prompt:
                binding_data = pd.concat([binding_data, binding_new])
                response = input("Do you want to continue adding files? (Y/n) ")
                if response.lower() == "n":
                    break
                print("The first five rows of your binding data look like this:")
                print(binding_data.head(5))
        return binding_data
    
    def modify_binding_data(self, binding_data):
        for column in binding_data.columns[1:]:
        # Convert non-numeric values to NaN
            binding_data[column] = pd.to_numeric(binding_data[column], errors='coerce')
        # Drop rows with NaN values
        binding_data = binding_data.dropna(axis=1, how='all')
        return binding_data
    
        
    def ask_binding_data(self, ):
        if not os.path.isfile(self.binding_data_dir):
            add_binding = input("Do you have binding Data? Y/n")
            if add_binding.lower() in ["Y", "y"]:
                binding_data = self.collect_binding_data()
                binding_data = self.modify_binding_data(binding_data)
                binding_data = binding_data.rename(columns={binding_data.columns[0]: "Sequences"})
                
                binding_data.to_csv("binding_data.csv")
            else:
                binding_data = None
        else:
            binding_data = pd.read_csv(self.binding_data_dir)
        if binding_data is not None:
            binding_data.drop_duplicates(subset = "Sequences", inplace = True)
        return binding_data



