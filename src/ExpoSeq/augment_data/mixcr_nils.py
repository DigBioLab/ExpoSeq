
import subprocess
from ast import literal_eval
import pickle
import pandas as pd
import os
from .create_mixcr_commands import get_basename 
import sys
#psutil
try:
    import tkinter as tk
    from tkinter import filedialog
except:
    pass
from ..augment_data.create_mixcr_commands import CreateCommand
from ..augment_data.collect_fastqs import CollectFastq
from ..settings.check_dirs import check_dirs

def check_mixcr(path_to_mixcr, data, settings_dir, testing = False):
    if path_to_mixcr == "":
        while True:
            print("choose the mixcr java file with the file chooser")
            try:
                path_to_mixcr = filedialog.askopenfilename()
            except:
                path_to_mixcr = input("copy and paste the path to the mixcr jar file here.")
                path_to_mixcr = os.path.normpath(path_to_mixcr)
            if testing == True:
                break

            try:
                print("Test Mixcr")
                subprocess.run(["java",
                                "-jar",
                                path_to_mixcr,
                                "--version"
                                ])
                confirmation = True
            except:
                print("Apparently you have not choosen the right path to the mixcr .jar file. Please try again")
                confirmation = False
            if confirmation == True:
                path_to_mixcr = os.path.abspath(path_to_mixcr)
                data["mixcr_path"] = path_to_mixcr
                with open(settings_dir, "w") as f:
                    f.write(str(data))
                break
    else:
        pass
    return path_to_mixcr

def process_mixcr(experiment, method, testing, paired_end_sequencing, add_args):
    module_dir = os.path.abspath("")
    check_dirs(module_dir, experiment)
    if not testing:
        Collect = CollectFastq(paired_end_sequencing)
        Collect.get_files()
        files = Collect.paired
        if paired_end_sequencing:
            while True:
                if len(files) > 0:
                    Collect.call_pairs()
                    correct_matches = input("Are the matches of forward and backward reads correct? Y/n")
                else:
                    correct_matches = "n"
                    print("You will need to restart the collection of the files if you want to proceed with paired end sequencing data.\nMake sure that you have the forward and backward files in two seperate folder and that you have the same number of files.")
                if correct_matches in ["Y", "y", "N", "n"]:
                    if correct_matches in ["Y", "y"]:
                        break
                    else:
                        while True:
                            choose_strat = input("You can restart the collection of the files or manually add them. \nPress 1 for restarting the collection by choosing the directories.\nPress 2 for restarting the collection and manually choose the forward and backward reads.")
                            if choose_strat in ["1", "2"]:
                                Collect = CollectFastq(paired_end_sequencing)
                                if choose_strat == "1":
                                    Collect.get_files()
                                    files = Collect.paired
                                    break
                                else:
                                    Collect.manually_match_pairs()
                                    files = Collect.paired
                                    break
                            else:
                                print("Please enter a valid value")
                else:
                    print("Please enter a valid value")
        else:
            pass
    else:
        return
    settings_dir = os.path.join(module_dir,
                                "settings",
                                "global_vars.txt")
    with open(settings_dir) as f:
        data = f.read()
    data = literal_eval(data)
    path_to_mixcr = data["mixcr_path"]
    path_to_mixcr = check_mixcr(path_to_mixcr, data, settings_dir)

    sequencing_report = pd.DataFrame([])
    while True:
        java_heap_size = input(
            "Enter the amount of memory you want to allocate for java in Mb. Normally 1000 is sufficient.")
        try:
            java_heap_size = int(java_heap_size)
            break
        except:
            print("Please enter the amount as integer")
    while True:
        mixcr_chain = input("\n\nPress enter if you want to choose the default value with heavy chain (IGH) only. Please, type\nIG: For all immunoglobulin chains\n   IGL: For Immunoglobulin lambda locus\n    IGK: For Immunoglobulin kappa locus \nTCR: For all T-cell receptor chains\n     TRA: Only TCR alpha\n   TRB: Only TCR beta.\n Note: If you have not sequenced the coresponding region the pipeline will fail to create the sequencing report.")
        if mixcr_chain == "" or mixcr_chain in ["IG", "IGH", "IGL", "IGK", "TCR", "TRA", "TRB", ]:
            if mixcr_chain == "":
                mixcr_chain = "IGH"
            break
        else:
            print("Please enter a correct value.")
    for filename in files:
        basename = get_basename(filename) # check if basename is correct and aligns to standards, e.g. 22 char
        if " " in basename:
            print("Spaces in the filename are not allowed. Analysis aborted")
            sys.exit(1)
    
    for filename in files:
        for file in os.listdir(os.path.join(module_dir, "temp")):
            os.remove(os.path.join(module_dir,
                                   "temp",
                                   file))
        Commands = CreateCommand(module_dir, path_to_mixcr, paired_end_sequencing, experiment, filename, java_heap_size)
        subprocess.run(Commands.prepare_align(method, add_args))
        subprocess.run(Commands.prepare_assembly())
        subprocess.run(Commands.prepare_clones(mixcr_chain = mixcr_chain))
        basename = Commands.basename
        try:
            columns_not_wanted = ['refPoints', 'allVHitsWithScore', 'allDHitsWithScore',
                                  'allJHitsWithScore', 'allCHitsWithScore', 'allVAlignments',
                                  'allDAlignments', 'allJAlignments', 'allCAlignments']
            clones_sample = pd.read_table(Commands.table_tsv)
            clones_sample.drop(columns=columns_not_wanted, inplace=True)
            clones_sample["Experiment"] = basename
            sequencing_report = pd.concat([sequencing_report, clones_sample])

            for file in os.listdir(os.path.join(module_dir, "temp")):
                os.remove(os.path.join(module_dir,
                                       "temp",
                                       file))
        except:
            print(f"Processing for {basename} failed")

    if not testing:
        
        if sequencing_report.shape[0] == 0:
            raise Exception(
                "The processing failed for all of your files. Most likely because you chose the wrong mixcr method.")
        report_dir = os.path.join(module_dir,
                                  "my_experiments",
                                  experiment,
                                  "sequencing_report.csv")
        sequencing_report.to_csv(report_dir)
        data["last_experiment"] = experiment
        glob_vars_dir = os.path.join(module_dir,
                                     "settings",
                                     "global_vars.txt")
        with open(glob_vars_dir, "w") as f:
            f.write(str(data))
        unique_experiments = sequencing_report["Experiment"].unique()
        experiment_dic = {item: item for item in list(unique_experiments)}
        exp_names_dir = os.path.join(module_dir,
                                     "my_experiments",
                                     experiment,
                                     "experiment_names.pickle")
        with open(exp_names_dir, "wb") as f:
            pickle.dump(experiment_dic, f)
    print(
        "The output of your processing is a sequencing_report which contains the sequences and their sample names of all fastq files. Further, a folder with the alignment reports for each fastq file was created. You can have a look under ~/my_experiments/YOUR_EXPERIMENT_NAME")

class RunMixcr():
    def __init__(self, experiment, method, paired_end_sequencing, add_args):
        self.module_dir = os.path.abspath("")
        self.experiment = experiment
        check_dirs(self.module_dir, self.experiment)
        path_to_mixcr, data = self.validate_mixcr(self.module_dir)
        self.files = self.prep_files(paired_end_sequencing)
        sequencing_report = self.process(path_to_mixcr, method, paired_end_sequencing, add_args)
        self.tidy(sequencing_report, data)
        
        
    @staticmethod
    def prep_files(paired_end_sequencing):
        Collect = CollectFastq(paired_end_sequencing)
        Collect.get_files()
        files = Collect.paired
        if paired_end_sequencing:
            while True:
                if len(files) > 0:
                    Collect.call_pairs()
                    correct_matches = input("Are the matches of forward and backward reads correct? Y/n")
                else:
                    correct_matches = "n"
                    print("You will need to restart the collection of the files if you want to proceed with paired end sequencing data.\nMake sure that you have the forward and backward files in two seperate folder and that you have the same number of files.")
                if correct_matches in ["Y", "y", "N", "n"]:
                    if correct_matches in ["Y", "y"]:
                        break
                    else:
                        while True:
                            choose_strat = input("You can restart the collection of the files or manually add them. \nPress 1 for restarting the collection by choosing the directories.\nPress 2 for restarting the collection and manually choose the forward and backward reads.")
                            if choose_strat in ["1", "2"]:
                                Collect = CollectFastq(paired_end_sequencing)
                                if choose_strat == "1":
                                    Collect.get_files()
                                    files = Collect.paired
                                    break
                                else:
                                    Collect.manually_match_pairs()
                                    files = Collect.paired
                                    break
                            else:
                                print("Please enter a valid value")
                else:
                    print("Please enter a valid value")
        else:
            pass
        return files

    @staticmethod
    def validate_mixcr(module_dir):
        settings_dir = os.path.join(module_dir,
                                "settings",
                                "global_vars.txt")
        with open(settings_dir) as f:
            data = f.read()
        data = literal_eval(data)
        path_to_mixcr = data["mixcr_path"]
        path_to_mixcr = check_mixcr(path_to_mixcr, data, settings_dir)
        return path_to_mixcr, data
    
    @staticmethod
    def ask_heap_size():
        while True:
            java_heap_size = input(
                "Enter the amount of memory you want to allocate for java in Mb. Normally 1000 is sufficient.")
            try:
                java_heap_size = int(java_heap_size)
                break
            except:
                print("Please enter the amount as integer")
        return java_heap_size
    
    @staticmethod
    def ask_chain():
        while True:
            mixcr_chain = input("\n\nPress enter if you want to choose the default value with heavy chain (IGH) only. Please, type\nIG: For all immunoglobulin chains\n   IGL: For Immunoglobulin lambda locus\n    IGK: For Immunoglobulin kappa locus \nTCR: For all T-cell receptor chains\n     TRA: Only TCR alpha\n   TRB: Only TCR beta.\n Note: If you have not sequenced the coresponding region the pipeline will fail to create the sequencing report.")
            if mixcr_chain == "" or mixcr_chain in ["IG", "IGH", "IGL", "IGK", "TCR", "TRA", "TRB", ]:
                if mixcr_chain == "":
                    mixcr_chain = "IGH"
                break
            else:
                print("Please enter a correct value.")
        return mixcr_chain
    
    def check_basenames(self):    
        for filename in self.files:
            basename = get_basename(filename) # check if basename is correct and aligns to standards, e.g. 22 char
            if " " in basename:
                print("Spaces in the filename are not allowed. Analysis aborted")
                sys.exit(1)
    
    
    def process(self,path_to_mixcr,  method, paired_end_sequencing, add_args):         
        
        java_heap_size = self.ask_heap_size()
        mixcr_chain = self.ask_chain()
        
        sequencing_report = pd.DataFrame([])        
        for filename in self.files:
            for file in os.listdir(os.path.join(self.module_dir, "temp")):
                os.remove(os.path.join(self.module_dir,
                                    "temp",
                                    file))
            Commands = CreateCommand(self.module_dir, path_to_mixcr, paired_end_sequencing, self.experiment, filename, java_heap_size)
            subprocess.run(Commands.prepare_align(method, add_args))
            subprocess.run(Commands.prepare_assembly())
            subprocess.run(Commands.prepare_clones(mixcr_chain = mixcr_chain))
            basename = Commands.basename
            try:
                columns_not_wanted = ['refPoints', 'allVHitsWithScore', 'allDHitsWithScore',
                                    'allJHitsWithScore', 'allCHitsWithScore', 'allVAlignments',
                                    'allDAlignments', 'allJAlignments', 'allCAlignments']
                clones_sample = pd.read_table(Commands.table_tsv)
                clones_sample.drop(columns=columns_not_wanted, inplace=True)
                clones_sample["Experiment"] = basename
                sequencing_report = pd.concat([sequencing_report, clones_sample])

                for file in os.listdir(os.path.join(self.module_dir, "temp")):
                    os.remove(os.path.join(self.module_dir,
                                        "temp",
                                        file))
            except:
                print(f"Processing for {basename} failed")
                
        return sequencing_report
    
    def tidy(self, sequencing_report, data):
        if sequencing_report.shape[0] == 0:
            raise Exception(
                "The processing failed for all of your files. Most likely because you chose the wrong mixcr method.")
        report_dir = os.path.join(self.module_dir,
                                  "my_experiments",
                                  self.experiment,
                                  "sequencing_report.csv")
        sequencing_report.to_csv(report_dir)
        data["last_experiment"] = self.experiment
        glob_vars_dir = os.path.join(self.module_dir,
                                     "settings",
                                     "global_vars.txt")
        with open(glob_vars_dir, "w") as f:
            f.write(str(data))
        unique_experiments = sequencing_report["Experiment"].unique()
        experiment_dic = {item: item for item in list(unique_experiments)}
        exp_names_dir = os.path.join(self.module_dir,
                                     "my_experiments",
                                     self.experiment,
                                     "experiment_names.pickle")
        with open(exp_names_dir, "wb") as f:
            pickle.dump(experiment_dic, f)
        print(
            "The output of your processing is a sequencing_report which contains the sequences and their sample names of all fastq files. Further, a folder with the alignment reports for each fastq file was created. You can have a look under ~/my_experiments/YOUR_EXPERIMENT_NAME")
