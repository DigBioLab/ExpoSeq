from .mixcr_nils import process_mixcr
from ast import literal_eval
import os
from .loop_collect_reports import load_mixed_files, load_alignment_reports
import pandas as pd
from glob import glob
from ..augment_data.check_reports import check_completeness
import sys
import pickle
from ..augment_data.trimming import trimming
try:
    import tkinter as tk
    from tkinter import filedialog
except:
    pass
import random
import string
import shutil
from ..settings.available_presets import data_list

def get_random_string(length):
    letters = string.ascii_letters  # this will get all the uppercase and lowercase letters
    result_str = ''.join(random.choice(letters) for i in range(length))
    return result_str

def pull_seq_align_repo(testing, testing_value):
    if not testing:
        try:
            print("choose the directory where you store your alignment reports and your sequencing report")
            filenames_dir = filedialog.askdirectory()
            assert os.path.isdir(filenames_dir), "Please enter a valid directory."
            filenames = glob(filenames_dir + "//*")
            
        except:
            filenames_dir = input(
                "Enter the directory where you store the alignment reports and your sequencing report manually.\nYou should have one file for the sequencing report and multiple files for the alignment reports.")
            while True:
                if os.path.isdir(filenames_dir):
                    break
                else:
                    print("Please enter a valid directory.")
                    filenames_dir = input(
                        "Enter the directory where you store the alignment reports and your sequencing report manually.\nYou should have one file for the sequencing report and multiple files for the alignment reports.")
            filenames_dir = os.path.abspath(filenames_dir + "//*")
            filenames = glob(filenames_dir)
    else:
        module_dir = os.path.abspath("")
        filenames_dir = os.path.join(module_dir,
                                    "tests",
                                    "test_files")
        filenames = glob(filenames_dir + "//*")
    for i in filenames:
        if os.path.isdir(i):
            if "alignment_reports" in i:
                alignment_repos = glob(i)
                filenames.extend(alignment_repos)
        sequencing_report, all_alignment_reports = load_mixed_files(filenames, testing, testing_value)
    return sequencing_report, all_alignment_reports

def check_experiment(module_dir, testing):
    if not testing:
        while True:
            experiment = input("How do you want to call your new experiment?")
            if os.path.isdir(os.path.join(module_dir, "my_experiments", experiment)) == True:
                replace = input("The given directory already exists. Do you want to replace it? (Y/n)")
                if replace in ["Y", "y", "n", "N"]:
                    break
                else:
                    print("Please enter another name.")
            else:
                replace = "n"
                break
        if replace in ["Y", "y"]:
            shutil.rmtree(os.path.join(module_dir, "my_experiments", experiment))

    else: 
        experiment = get_random_string(10)
    
    os.mkdir(os.path.join(module_dir,
                    "my_experiments",
                    experiment))
    return experiment

def create_alignment_report(module_dir, experiment):
    try:
        align_repo_path = os.path.join(module_dir,
                                       "my_experiemnts",
                                       experiment,
                                       "alignment_reports",
                                       "*")
        alignment_reports = glob(align_repo_path)
        all_alignment_reports = load_alignment_reports(alignment_reports)
    except:
        all_alignment_reports = pd.DataFrame([])
        print(
            "No alignment reports could be found in " + experiment + ". You will continue without being able to analyze the Alignment Quality.")
    return all_alignment_reports

def method_one(experiment, repo_path, module_dir, testing, paired_end_test = "n"):
    if not testing:
        for i in data_list:
            print(i)
        use_method = input(
            "Per default you will align your data with the following method: milab-human-tcr-dna-multiplex-cdr3.\nPress enter if you want to continue.\nOtherwise type in the method of your choice.\nIt has to be the exact same string which is given on the Mixcr documentation.The available presets are:")
        if use_method == "":
            use_method = "milab-human-tcr-dna-multiplex-cdr3"
        else:
            while True:
                if use_method in data_list:
                    break
                use_method = input("The method you entered is not in the list. Please enter a valid method.")
                if use_method in data_list:
                    break
                else:
                    print("Please enter a valid method.")
    else:
        use_method = ""
    if use_method == "":
        method = "milab-human-tcr-dna-multiplex-cdr3"
    else:
        method = use_method
    if not testing:
        paired_end = input("Do you want to analyze paired_end_sequencing data? (Y/n)")
        
    else: 
        paired_end = 'n'
    while True:
        if paired_end in ["Y", "y", "n", "N"]:
            break
        else:
            print("Please enter a correct value.")
            paired_end = input("Do you want to analyze paired_end_sequencing data? (Y/n)")
    if paired_end.lower() in ["Y", "y"]:
        paired_end = True
    else:
        paired_end = False
    if repo_path == "":
        repo_path = os.path.join(module_dir,
                                "my_experiments",
                                experiment,
                                "sequencing_report.csv")
    process_mixcr(experiment,
                  method=method,
                  testing = testing,
                  paired_end_sequencing=paired_end)
    repo_path = os.path.join(module_dir,
                             "my_experiments",
                             experiment,
                             "sequencing_report.csv")
    with open(repo_path, "rb") as f:
        sequencing_report = pd.read_table(f, sep=",")
    all_alignment_reports = create_alignment_report(module_dir, experiment)
    return sequencing_report, all_alignment_reports

def method_two(module_dir, experiment, testing = False, testing_value = None):

    print("Choose the folder which contains the sequencing_report and the alignment reports.")

    sequencing_report, all_alignment_reports = pull_seq_align_repo(testing, testing_value)

    while True:
        if sequencing_report.shape[0] == 0:
            print(
                "Something went wrong. Unfortunately you could not collect the data in the sequencing report.")
            sequencing_report, all_alignment_reports = pull_seq_align_repo(testing, testing_value)
        else:
            break
        #  sys.exit()
    if not testing:
        while True:
            trim_data = input("Do you need to trim your data? (Y/n)")
            if trim_data in ["Y", "y", "n", "N"]:
                break
            else:
                print("Please enter a correct value.")
    else: trim_data = "Y"
    if trim_data.lower() in ["Y", "y"]:
        sequencing_report = trimming(sequencing_report,
                                     divisible_by=3,  # might be possible that columns are different
                                     min_count=1,
                                     new_fraction="cloneFraction")
    else:
        pass
    seq_report_dir = os.path.join(module_dir,
                                  "my_experiments",
                                  experiment,
                                  "sequencing_report.csv")
    sequencing_report.to_csv(seq_report_dir)
    unique_experiments = sequencing_report["Experiment"].unique()
    experiment_dic = {item: item for item in list(unique_experiments)}
    exp_names_dir = os.path.join(module_dir,
                                 "my_experiments",
                                 experiment,
                                 "experiment_names.pickle")
    with open(exp_names_dir, "wb") as f:
        pickle.dump(experiment_dic, f)

    if all_alignment_reports.shape[0] == 0:
        print(
            "No Alignment Reports were uploaded. You will continue the analysis without being able to analyze the Alignment Quality.")
    else:
        try:
            all_alignment_reports = check_completeness(all_alignment_reports, sequencing_report)
        except:
            pass
    exp_name_path = os.path.join(module_dir,
                                 "my_experiments",
                                 experiment,
                                 "experiment_names.pickle")

    unique_experiments = sequencing_report["Experiment"].unique()
    experiment_dic = {item: item for item in list(unique_experiments)}

    with open(exp_name_path, "wb") as f:
        pickle.dump(experiment_dic, f)
    return sequencing_report, all_alignment_reports

def method_three(module_dir, testing ):
    
    if not testing:
        try:
            print("choose the directory with your experiment name and which has the sequencing_report.csv file in its root.")
            while True:
                f = filedialog.askdirectory()
                if os.path.isdir(f):
                    if "sequencing_report.csv" in os.listdir(f):
                        break
                    
                else:
                    print("Please enter a valid path to the directory")
        except:
            f = input("give the path to the directory with the sequencing_report.csv file in its root.")
            f = os.path.abspath(f)
            while True:
                if os.path.isdir(f):
                    if "sequencing_report.csv" in os.listdir(f):
                        break
                else:
                    print("Please enter a valid path to the directory")
        experiment = os.path.basename(f)
        shutil.copytree(f, os.path.join(module_dir, "my_experiments", experiment))
        f = os.path.join(f, "sequencing_report.csv")
    else:
        f = os.path.join(module_dir,
                        "my_experiments",
                        "test_directory",
                        "sequencing_report.csv")
        experiment = "test_directory"
    sequencing_report = pd.read_table(f, sep=",")
    if not testing:
        all_alignment_reports = None
    else:
        all_alignment_reports = None
    return sequencing_report,all_alignment_reports, experiment


def find_file_in_subdirectories(directory, filename):
    for dirpath, dirnames, files in os.walk(directory):
        if filename in files:
            if os.path.basename(dirpath) != "test_directory":
                return os.path.join(dirpath, filename)
    return None

def check_last_exp(module_dir, testing):
    if not testing:
        glob_vars = os.path.join(module_dir,
                                "settings",
                                "global_vars.txt")
        with open(glob_vars, "r") as f:
            data = f.read()
        data = literal_eval(data)

        last_experiment = data["last_experiment"]
        if last_experiment == "":
            if not bool(os.listdir(os.path.join(module_dir, "my_experiments"))):
                repo_path = ""
                print("You have no experiments in your folder.")
            else:
                subdirectory = find_file_in_subdirectories(os.path.join(module_dir, "my_experiments"), "sequencing_report.csv")
                if subdirectory is not None:
                    repo_path = subdirectory
                    dirpath = os.path.dirname(subdirectory)
                    last_experiment = os.path.basename(dirpath)
                else:
                    print("The experiments in my experiment do not contain the sequencing report.\nYou need to upload a new experiment.")
                    repo_path = ""
            
        else:
            repo_path = os.path.join(module_dir,
                                   "my_experiments",
                                    last_experiment,
                                    "sequencing_report.csv")
    else:
        repo_path = os.path.join(module_dir,
                                "my_experiments",
                                "test_directory",
                                "sequencing_report.csv")
        last_experiment = "test_directory"
    return repo_path, last_experiment



def get_continue_analysis_input(last_experiment, testing=False, testing_value = None):
    if not testing:
        while True:
            continue_input = input("Do you want to continue to analyze with " + last_experiment + "? Y/n")
            if continue_input in ["Y", "y", "n", "N"]:
                break
            else:
                print("Please enter a correct value.")
    else:
        continue_input = testing_value
    return continue_input

def get_next_step_input(testing=False, testing_value = None):
    if not testing:
        while True:
            next_step = input("If you want to upload a new experiment press 1.\nIf you want to choose another experiment press 2. ")
            if next_step in ["1", "2"]:
                break
            else:
                print("Please enter a correct value.")
    else:
        next_step = testing_value
    return next_step

def get_choose_method_input(testing, testing_value = None):
    if not testing:
        while True:
            processing = input("If you want to process your fastq files with mixcr press 1.\nIf you have processed your files externally and want to upload the directory with the Experiment name and the processing data, press 2.")
            if processing in ["1", "2"]:
                break
            else:
                print("Please enter a correct value.")
    else:
        processing = testing_value
    return processing

def get_experiment_name_input(testing = False):
    if not testing:
        experiment_name = input("Enter the name of the experiment you want to analyze")
    else:
        experiment_name = get_random_string(10)
    return experiment_name  

def upload_new_experiment(module_dir, repo_path, testing, testing_value, paired_end_test, experiment_column):
    while True:
        choose_method = get_choose_method_input(testing, testing_value)
        if choose_method == "1":
            experiment = check_experiment(module_dir, testing)
            
        if choose_method in ["1", "2", "3"]:
            if choose_method == "1":
            
                sequencing_report, all_alignment_reports = method_one(experiment, repo_path,module_dir,testing, paired_end_test )
         #   elif choose_method == "2":
          #      sequencing_report, all_alignment_reports = method_two(module_dir, experiment,testing, experiment_column)
            elif choose_method == "2":
                sequencing_report, all_alignment_reports, experiment = method_three(module_dir, testing)
            break
        else:
            print("Please enter a correct value.")

        

    return sequencing_report, all_alignment_reports, experiment

def choose_existing_experiment(module_dir, testing):
    if not testing:
        while True:
            user_input = get_experiment_name_input()
            spec_exp_name_path = os.path.join(module_dir, "my_experiments", user_input)
            if os.path.isdir(spec_exp_name_path): 
                break
            else:
                print("The experiment name does not exist in my_experiments. Please enter the correct name")
    else:
        user_input = "test_directory"
    return user_input  # experiment



def retrieve_reports(module_dir, experiment, testing):
    if not testing:
        seq_report_path = os.path.join(module_dir, "my_experiments", experiment, "sequencing_report.csv")
        align_repo_path = os.path.join(module_dir,
                            "my_experiments",
                            experiment,
                            "alignment_reports", "*")
    else:
        seq_report_path = os.path.join(module_dir, "my_experiments", "test_directory", "sequencing_report.csv")
        align_repo_path = os.path.join(module_dir, "tests", "test_directory", "all_alignment_reports.pickle")
    if os.path.isfile(seq_report_path):
        with open(seq_report_path, "rb") as f:
            sequencing_report = pd.read_table(f, sep=",")
        try:
            align_repo_path = os.path.join(module_dir,
                                        "my_experiments",
                                        experiment,
                                        "alignment_reports", "*")
            alignment_reports = glob(align_repo_path)
            all_alignment_reports = load_alignment_reports(alignment_reports)
        except:
            all_alignment_reports = pd.DataFrame([])
            print("No alignment reports could be found in " + experiment + ".\nYou will continue without being able to analyze the Alignment Quality.")
    else:
        print("The sequencing report could not be found in " + experiment + ".\nPlease enter an experiment with an existing sequencing_report.csv file or upload new data.")
        sequencing_report = None
        all_alignment_reports = None
        experiment = None

    return sequencing_report, all_alignment_reports, experiment

def process_data_with_last_experiment(module_dir, last_experiment, testing):
    experiment = last_experiment
    sequencing_report, all_alignment_reports, experiment = retrieve_reports(module_dir, experiment, testing)
    return sequencing_report, all_alignment_reports, experiment



def write_last_exp(module_dir, experiment):
    glob_vars = os.path.join(module_dir,
                             "settings",
                             "global_vars.txt")
    with open(glob_vars, "r") as f:
        data = f.read()
    data = literal_eval(data)
    data["last_experiment"] = experiment
    with open(glob_vars, "w") as f:
        f.write(str(data))

def upload(testing=False, continue_analysis = "n", upload_type = "2", choose_exp = '1', paired_end_test = 'n', experiment_column = '1'):
    module_dir = os.getcwd()

    repo_path, last_experiment = check_last_exp( module_dir, testing)
    
    if os.path.isfile(repo_path):
        continue_analysis = get_continue_analysis_input(last_experiment,
                                                        testing,
                                                        continue_analysis)
        if continue_analysis.lower() in ["n", "no"]:
            next_step = get_next_step_input(testing, choose_exp)
            while True:
                if next_step == "1":
                    
                    sequencing_report, all_alignment_reports, experiment = upload_new_experiment(module_dir,
                                                                                                repo_path,
                                                                                                testing,
                                                                                                upload_type,
                                                                                                paired_end_test,
                                                                                                experiment_column)
                elif next_step == "2":
                    experiment = choose_existing_experiment(module_dir,
                                                            testing)
                    sequencing_report, all_alignment_reports, experiment = retrieve_reports(module_dir,
                                                                                            experiment,
                                                                                            testing)
                if isinstance(sequencing_report, pd.DataFrame):
                    break
                else:
                    print("Something went wrong. Please try again.")
                    

        else:
            sequencing_report, all_alignment_reports, experiment = process_data_with_last_experiment(module_dir,
                                                                                                     last_experiment,
                                                                                                     testing)
    else:
        sequencing_report, all_alignment_reports, experiment = upload_new_experiment(module_dir,
                                                                                     repo_path,
                                                                                     testing,
                                                                                     upload_type,
                                                                                     paired_end_test,
                                                                                     experiment_column)
    if testing and experiment != "test_directory":
        shutil.rmtree(os.path.join(module_dir,
                                   "my_experiments",
                                   experiment))
    
    if not testing and experiment != "test_directory":
        write_last_exp(module_dir, experiment)
    return sequencing_report, all_alignment_reports, experiment 