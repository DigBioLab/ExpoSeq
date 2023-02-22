import os
import sys
from glob import glob
from tkinter import filedialog
import subprocess
from os.path import dirname, abspath
from ast import literal_eval
import pickle
import pandas as pd
from .trimming import trimming


def add_fastq_files():
    filenames = []

    while True:
        path_to_files = filedialog.askdirectory()
        filenames.extend(glob(path_to_files + "/*.fastq"))
        print("These are the files you chose:")
        for i in filenames:
            print(os.path.basename(i))
        read_more_files = input("do you want to add another folder? (Y/n)")

        if read_more_files == "Y" or read_more_files == "y":
            continue
        else:
            break
    sorted(filenames)

    return filenames


def process_mixcr(experiment,method, paired_end_sequencing = False):
    print("Choose the directory where you store your fastq files")
    filenames = add_fastq_files()
    with open('settings/global_vars.txt') as f:
        data = f.read()
    data = literal_eval(data)
    path_to_mixcr = data["mixcr_path"]
    if path_to_mixcr == "":
        while True:
            print("choose the mixcr java file with the file chooser")
            path_to_mixcr = filedialog.askopenfilename()
            confirmation = input("Is this the right file? Please make sure that this is the right path." +path_to_mixcr + "Otherwise, you might face issues with conducting the further analysis. Type Y or n")
            if confirmation.lower() in ["Y", "y"]:
                data["mixcr_path"] = path_to_mixcr
                with open("settings/global_vars.txt", "w") as f:
                    f.write(str(data))
                break
    else:
        pass
    #kAligner2_4.0
    #nebnext-human-bcr-cdr3
    #milab-human-tcr-dna-multiplex-cdr3
    while True:
        if not os.path.isdir("my_experiments/" + experiment):
            os.mkdir("my_experiments/" + experiment)
            os.mkdir("my_experiments/" + experiment + "/alignment_reports")  # raise error if already exists
            break
        else:
            experiment = input("The given experiment name already exists. Please type another one.")

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
            sys.exit()
    else:
        combined_filenames = filenames

    for filename in combined_filenames:
        basename = os.path.basename(os.path.splitext(filename)[0])
        basename = basename[:len(basename) - 2]
        filename_base = basename
        result = "temp/" + filename_base + ".vdjca"
        subprocess.run(["java",
                        "-jar",
                        path_to_mixcr,
                        "align",
                        "-p " + method,
                        filename,
                        result,
                        "-r" + "my_experiments/" + experiment + "/alignment_reports/" + filename_base + "_AlignmentReport.txt"])

        clones = "temp/" + basename + "clones.clns"
        subprocess.run(["java",
                        "-jar",
                        path_to_mixcr,
                        "assemble",
                        "-OseparateByC=true",
                        "-OseparateByV=true",
                        "-OseparateByJ=true",
                        result,
                        "temp/" + basename + "clones.clns"])

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
                        "temp/" + basename + ".tsv"
                        ])

        clones_sample = pd.read_table(r"temp/" + basename + "_IGH.tsv")
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
        clones_sample = trimming(clones_sample, divisible_by = 3, min_count = 3, new_fraction = "clonesFraction")
        sequencing_report = pd.concat([sequencing_report, clones_sample])
        files_to_remove = os.listdir("temp")
        for file in files_to_remove:
            os.remove("temp/" + file)

    sequencing_report.to_csv("my_experiments/" + experiment + "/sequencing_report.txt")
    data["last_experiment"] = experiment
    with open("settings/global_vars.txt", "w") as f:
        f.write(str(data))
    unique_experiments = sequencing_report["Experiment"].unique()
    experiment_dic = {item: item for item in list(unique_experiments)}
    with open("my_experiments/" + experiment + "/experiment_names.pickle", "wb") as f:
        pickle.dump(experiment_dic, f)
