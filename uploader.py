from mixcr_nils import process_mixcr
from ast import literal_eval
import os
from python_scripts.augment_data.load_data import collectFiles
from python_scripts.augment_data.loop_collect_reports import load_mixed_files, load_alignment_reports
from python_scripts.tidy_data.trimming import trimming
import pandas as pd
from glob import glob
from python_scripts.augment_data.check_reports import check_completeness
import sys
import pickle

def upload():
    with open('global_vars.txt', "r") as f:
        data = f.read()
    data = literal_eval(data)

    last_experiment = data["last_experiment"]
    if os.path.isfile("my_experiments/" + last_experiment + "/" + "sequencing_report.txt"):
        continue_analysis = input("Do you want to continue to analyze with " + last_experiment + "? Y/n")
        if continue_analysis.lower() in ["n", "N"]:
            next_step = input("If you want to upload a new experiment press 1. If you want to choose another experiment press 2")
            if next_step == "1":
                experiment = input("How do you want to call your new experiment?")
                choose_method = input("If you want to process your data with mixcr press 1. If you want to upload already processed txt files with unique clones, press 2")
                if choose_method == "1":
                    use_method = input("Per default you will align your data with the following method: milab-human-tcr-dna-multiplex-cdr3 . Press enter if you want to continue. Otherwise type in the method of your choice which. It has to be the exact same string which is given on the Mixcr documentation.")
                    if use_method == "":
                        method = "milab-human-tcr-dna-multiplex-cdr3"
                    else:
                        method = use_method
                    paired_end = input("Do you want to analyze paired_end_sequencing data? (Y/n)")
                    if paired_end.lower() in ["Y", "y"]:
                        paired_end = True
                    else:
                        paired_end = False
                    process_mixcr(experiment,
                                    method = method,
                                    paired_end_sequencing = paired_end)
                    with open("my_experiments/" + experiment + "/sequencing_report.txt", "rb") as f:
                        sequencing_report = pd.read_table(f, sep=",")
                    try:
                        alignment_reports = glob("my_experiments/" + experiment + "/alignment_reports/*")
                        all_alignment_reports = load_alignment_reports(alignment_reports)
                    except:
                        all_alignment_reports = pd.DataFrame([])
                        print("No alignment reports could be found in " + experiment + ". You will continue without being able to analyze the Alignment Quality.")
                if choose_method == "2":
                    experiment = input("How do you want to call your experiment?")
                    print("Choose the folder which contains the txt files")
                    filenames = collectFiles()
                    try:
                        sequencing_report, all_alignment_reports = load_mixed_files(filenames)
                     #   trim_data = input("Do you need to trim your data? (Y/n)")
                    #        if trim_data.lower() in ["Y", "y"]:
                    #           sequencing_report = trimming(sequencing_report,
                    #                                       divisible_by=3,                   # might be possible that columns are different
                    #                                      min_count=1,
                    #                                     new_fraction="cloneFraction")
                    #   else:
                    #      pass
                    except:
                        sys.exit()
                    if all_alignment_reports.shape[0] == 0:
                        print(
                            "No Alignment Reports were uploaded. You will continue the analysis without being able to analyze the Alignment Quality.")
                    else:
                        all_alignment_reports = check_completeness(all_alignment_reports, sequencing_report)

            if next_step == "2":
                while True:
                    user_input = input("Enter the name of the experiment you want to analyze")

                    user_input = user_input  # Try to convert the input to an integer
                    if os.path.isdir("my_experiments/" + user_input):  # Check if the input is in the correct range
                        break  # If the input is valid, break out of the loop
                    else:
                        print("The experiment name does not exist in my_experiments. Please enter the correct name")
                experiment = user_input
                with open("my_experiments/" + experiment + "/sequencing_report.txt", "rb") as f:
                    sequencing_report = pd.read_table(f, sep=",")
                try:
                    alignment_reports = glob("my_experiments/" + experiment + "/alignment_reports/*")
                    all_alignment_reports = load_alignment_reports(alignment_reports)
                except:
                    all_alignment_reports = pd.DataFrame([])
                    print(
                        "No alignment reports could be found in " + experiment + ". You will continue without being able to analyze the Alignment Quality.")

            else:
                pass
        else:
            experiment = last_experiment
            with open("my_experiments/" + experiment + "/sequencing_report.txt", "rb") as f:
                sequencing_report = pd.read_table(f, sep=",")
            try:
                alignment_reports = glob("my_experiments/" + experiment + "/alignment_reports/*")
                all_alignment_reports = load_alignment_reports(alignment_reports)
            except:
                all_alignment_reports = pd.DataFrame([])
                print(
                    "No alignment reports could be found in " + experiment + ". You will continue without being able to analyze the Alignment Quality.")

    else:
        experiment = input("How do you want to call your experiment?")
        choose_method = input("Welcome to ExpoSeq. If you want to start by processing your fastq files with mixcr press 1. If you want to upload already processed files press 2.")

        if choose_method == "1":
            use_method = input(
                "Per default you will align your data with the following method: milab-human-tcr-dna-multiplex-cdr3 . Press enter if you want to continue. Otherwise type in the method of your choice which. It has to be the exact same string which is given on the Mixcr documentation.")
            if use_method == "":
                method = "milab-human-tcr-dna-multiplex-cdr3"
            else:
                method = use_method
            paired_end = input("Do you want to analyze paired_end_sequencing data? (Y/n)")
            if paired_end.lower() in ["Y", "y"]:
                paired_end = True
            else:
                paired_end = False
            process_mixcr(experiment,
                          method=method,
                          paired_end_sequencing=paired_end)
            with open("my_experiments/" + experiment + "/sequencing_report.txt", "rb") as f:
                sequencing_report = pd.read_table(f, sep=",")
            try:
                alignment_reports = glob("my_experiments/" + experiment + "/alignment_reports/*")
                all_alignment_reports = load_alignment_reports(alignment_reports)
            except:
                all_alignment_reports = pd.DataFrame([])
                print(
                    "No alignment reports could be found in " + experiment + ". You will continue without being able to analyze the Alignment Quality.")
        else:
            pass
        if choose_method == "2":
            experiment = input("How do you want to call your experiment")
            print("Choose the folder which contains the txt files with the Alignment Reports and the Sequencing Reports. Please avoid having other text files in this folder.")
            filenames = collectFiles()
            try:
                sequencing_report, all_alignment_reports = load_mixed_files(filenames)
             #   trim_data = input("Do you need to trim your data? (Y/n)")
        #        if trim_data.lower() in ["Y", "y"]:
         #           sequencing_report = trimming(sequencing_report,
          #                                       divisible_by=3,                   # might be possible that columns are different
           #                                      min_count=1,
            #                                     new_fraction="cloneFraction")
             #   else:
              #      pass
            except:
                sys.exit()
            if all_alignment_reports.shape[0] == 0:
                print(
                    "No Alignment Reports were uploaded. You will continue the analysis without being able to analyze the Alignment Quality.")
            else:
                all_alignment_reports = check_completeness(all_alignment_reports, sequencing_report)
            try:
                unique_experiments = sequencing_report["Experiment"].unique()
                experiment_dic = {item: item for item in list(unique_experiments)}
                with open("my_experiments/" + experiment + "/experiment_names.pickle", "wb") as f:
                    pickle.dump(experiment_dic, f)
            except:
                col_name = input("Please type in the exact column name which contains the name of the samples")
                unique_experiments = sequencing_report[col_name].unique()
                experiment_dic = {item: item for item in list(unique_experiments)}
                with open("my_experiments/" + experiment + "/experiment_names.pickle", "wb") as f:
                    pickle.dump(experiment_dic, f)
        else:
            pass


    return sequencing_report, all_alignment_reports, experiment