from mixcr_nils import process_mixcr
from ast import literal_eval
import os
from python_scripts.augment_data.load_data import collectFiles
from python_scripts.augment_data.loop_collect_reports import collect_nocluster_files, load_alignment_reports
from python_scripts.tidy_data.trimming import trimming
import pandas as pd

def upload():
    with open('global_vars.txt', "r") as f:
        data = f.read()
    data = literal_eval(data)
    try:
        last_experiment = data["last_experiment"]
        if os.path.isfile("my_experiments/" + last_experiment + "/" + "sequencing_report.txt"):
            continue_analysis = input("Do you want to continue to analyze with " + last_experiment + "? Y/n")

            if continue_analysis.lower() in ["n", "N"]:
                next_step = input("If you want to upload a new experiment press 1. If you want to choose another experiment press 2")
                if next_step == "1":
                    experiment = input("How do you want to call your new experiment?")
                    choose_method = input("If you want to process your data with mixcr press 1. If you want to upload already processed txt files with unique clones, press 2")
                    if choose_method == "1":
                        process_mixcr(experiment)
                        with open("my_experiments/" + experiment + "/sequencing_report.txt", "rb") as f:
                            sequencing_report = pd.read_table(f, sep=",")
                    if choose_method == "2":
                        print("Choose the folder which contains the txt files")
                        filenames = collectFiles()
                        try:
                            sequencing_report, each_instance = collect_nocluster_files(filenames)
                            sequencing_report = trimming(sequencing_report,
                                                        divisible_by=3,
                                                        min_count=1,
                                                        new_fraction="cloneFraction")

                        except ValueError:
                            print("It seems that in your directory are either txt files which are not dataframes or you have not uploaded any")
                        try:
                            all_alignment_reports = load_alignment_reports(filenames)
                        except ValueError:
                            print("No files could be found")
                if next_step == "2":
                    while True:
                        user_input = input("Enter the name of the experiment you want to analyze")

                        user_input = user_input  # Try to convert the input to an integer
                        if os.path.isdir("my_experiments/" + user_input):  # Check if the input is in the correct range
                            break  # If the input is valid, break out of the loop
                        else:
                            print("The experiment name does not exist in my_experiments. Please enter the correct name")
                else:
                    pass

            else:
                with open("my_experiments/" + last_experiment + "/sequencing_report.txt", "rb") as f:
                    sequencing_report = pd.read_table(f, sep=",")

                with open("my_experiments" + last_experiment + "/sequencing_report.txt", "rb") as f:

        else:
            pass
    except:
        experiment = input("How do you want to call your new experiment?")
        process_mixcr(experiment)
        with open("my_experiments/" + experiment + "/sequencing_report.txt", "rb") as f:
            sequencing_report = pd.read_table(f, sep=",")

    return sequencing_report, data