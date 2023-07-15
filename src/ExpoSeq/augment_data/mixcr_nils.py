import os
import sys
from glob import glob
import subprocess
from ast import literal_eval
import pickle
import pandas as pd
from ExpoSeq.augment_data.trimming import trimming
import pkg_resources
try:
    import tkinter as tk
    from tkinter import filedialog
except:
    pass
def add_fastq_files():
    filenames = []
    while True:
        try:
            while True:
                path_to_files = filedialog.askdirectory()
                if os.path.isdir(path_to_files) == True:
                    break
                else:
                    pass
            
        except:
            while True:
                path_to_files = input("copy and paste the path to your directory")
                if os.path.isdir(os.path.abspath(path_to_files)) == True:
                    break
                else:
                    print("Please enter a valid directory. Dont add a \ to the end of the directory path")
        path_to_files = os.path.abspath(path_to_files)
        filenames.extend(glob(os.path.join(path_to_files, "*.fastq")))
        if len(filenames) == 0:
            user_choice = input("You have not chosen a directory with fastq files. If you want to try again press 1. If you want to cancel the analysis press 2")
            if user_choice == str(1):
                read_more_files = "Y"
            else:
                read_more_files = "N"
        else:
            print("These are the files you chose:")
            for i in filenames:
                print(os.path.basename(i))
            while True:
                read_more_files = input("do you want to add another folder? (Y/n)")
                if read_more_files in ["Y" "y", "n", "N"]:
                    break
                else:
                    pass
        if read_more_files == "Y" or read_more_files == "y":
            continue
        else:
            break
    sorted(filenames)

    return filenames


def process_mixcr(experiment, method, testing, paired_end_sequencing):
    pkg_path = pkg_resources.resource_filename("ExpoSeq", "")
    module_dir = os.path.abspath("")
    
    if not testing:
        print("Choose the directory where you store your fastq files")
        while True:
            filenames = add_fastq_files()
            if len(filenames) == 0:
                print("You have not added any fastq files the preprocessing is cancelled now.")
            else:
                break
        if not os.path.isdir(os.path.join(module_dir, "my_experiments", experiment, "alignment_reports")):
            os.mkdir(os.path.join(module_dir, "my_experiments", experiment, "alignment_reports"))

    else:

        filenames = glob(os.path.join(pkg_path,"test_data","test_files", "*.fastq"))
        print(filenames)
        experiment = "test_directory"
    
    settings_dir = os.path.join(pkg_path,
                                "settings",
                                "global_vars.txt")
    with open(settings_dir) as f:
        data = f.read()
    data = literal_eval(data)
    path_to_mixcr = data["mixcr_path"]
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
    #kAligner2_4.0
    #nebnext-human-bcr-cdr3
    #milab-human-tcr-dna-multiplex-cdr3


    sequencing_report = pd.DataFrame([])
    if paired_end_sequencing == True:
        filenames_unique_reverse = [unique for unique in filenames if "2.fastq" in unique]
        filenames_unique_forward = [unique for unique in filenames if "1.fastq" in unique]
        combined_filenames = []
        for i in filenames_unique_reverse:
            file_one = i
            basename = os.path.basename(os.path.splitext(file_one)[0])
            basename = basename[:len(basename) - 2]
            file_two = ""
            for j in filenames_unique_forward:
                if basename in filenames_unique_forward:
                    file_two = j
                    break
            if file_two == "":
                print("you are missing the the reverse strand for " + basename)
            else:
                fastq_files = file_one + " " + file_two
                combined_filenames.append(fastq_files)
        if len(combined_filenames) == 0:
            print("the system couldnt find any forward reverse matches. Please check your given directory or continue with single-end analysis.")
            return
    else:
        combined_filenames = filenames

    for filename in combined_filenames:
        basename = os.path.basename(os.path.splitext(filename)[0])
        basename = basename[:len(basename) - 2]
        filename_base = basename
        result = os.path.join(module_dir, "temp", filename_base + ".vdjca")
        clones = os.path.join(module_dir, "temp", basename + "clones.clns")
        subprocess.run(["java",
                        "-jar",
                        path_to_mixcr,
                        "align",
                        "-p " + method,
                        os.path.normpath(filename),
                        result,
                        "--report",
                        os.path.normpath(os.path.join(module_dir,
                                           "my_experiments",
                                           experiment,
                                           "alignment_reports",
                                           filename_base + "_AlignmentReport.txt"))
                        ])

        subprocess.run(["java",
                        "-jar",
                        path_to_mixcr,
                        "assemble",
                        "-OseparateByC=true",
                        "-OseparateByV=true",
                        "-OseparateByJ=true",
                        result,
                        clones
                        ])

        subprocess.run(["java",
                        "-jar",
                        path_to_mixcr,
                        "exportClones",
                        "-cloneId",
                        "-readCount",
                        "-readFraction",
                        "-lengthOf CDR3",
                        "-nFeature CDR3",
                        "-aaFeature CDR3",
                        "-avrgFeatureQuality CDR3",
                        clones,
                        os.path.join(module_dir,
                                     "temp",
                                     basename + ".tsv")
                        ])

        clones_sample = pd.read_table(os.path.join(module_dir,
                                                   "temp",
                                                   basename + "_IGH.tsv"))
        clones_sample = clones_sample[["cloneId",
                                        "readCount",
                                        "readFraction",
                                        "nSeqCDR3",
                                        "aaSeqCDR3",
                                        "minQualCDR3",
                                        "lengthOfCDR3",
                                        "meanQualCDR3"]]
        new_fractions = clones_sample.groupby("nSeqCDR3")["readFraction"].sum().reset_index()
        clones_sample = clones_sample.drop_duplicates(subset=["nSeqCDR3"], keep="first")
        clones_sample = new_fractions.merge(clones_sample,
                                            how = "left",
                                            on = "nSeqCDR3")
        clones_sample = clones_sample.sort_values(by="cloneId")
        clones_sample = clones_sample.reset_index()
        clones_sample = clones_sample.drop(columns = ["readFraction_y", "index"])
        clones_sample = clones_sample.rename(columns={"readFraction_x": "cloneFraction"})
        clones_sample["Experiment"] = filename_base
        clones_sample = trimming(clones_sample,
                                 divisible_by = 3,
                                 min_count = 3,
                                 new_fraction = "clonesFraction")
        sequencing_report = pd.concat([sequencing_report, clones_sample])
        files_to_remove = os.listdir(os.path.join(module_dir, "temp"))
        for file in files_to_remove:
            os.remove(os.path.join(module_dir,
                                   "temp",
                                   file))
    if not testing:
        report_dir = os.path.join(module_dir,
                                  "my_experiments",
                                  experiment,
                                  "sequencing_report.csv")
        sequencing_report.to_csv(report_dir)
        data["last_experiment"] = experiment
        glob_vars_dir = os.path.join(pkg_path,
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

