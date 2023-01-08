import os
from glob import glob
from tkinter import filedialog
import subprocess
from os.path import dirname, abspath
from ast import literal_eval

import pandas as pd


def add_fastq_files():
    filenames = []

    while True:
        path_to_files = filedialog.askdirectory()
        filenames.extend(glob(path_to_files + "/*.fastq"))
        read_more_files = input("do you want to add another folder? (Y/n)")
        if read_more_files == "Y" or read_more_files == "y":
            continue
        else:
            break
    sorted(filenames)
    return filenames


def process_mixcr(experiment):
    print("Choose the directory where you store your fastq files")
    filenames = add_fastq_files()
    ROOT = dirname(abspath('ExpoSeq'))
    with open(ROOT + 'global_vars.txt') as f:
        data = f.read()
    data = literal_eval(data)
    path_to_mixcr = data["mixcr_path"]
    if path_to_mixcr == "":
        print("choose the mixcr java file with the file chooser")
        path_to_mixcr = filedialog.askopenfilename()
        data["mixcr_path"] = path_to_mixcr
        f = open(ROOT + 'global_vars.txt' 'w')
        f.write(str(data))
    else:
        pass

    #kAligner2_4.0
    #nebnext-human-bcr-cdr3
    #milab-human-tcr-dna-multiplex-cdr3

    paired_end_sequencing = False
    filenames_unique = [unique for unique in filenames if "_1.fastq" in unique]
    sequencing_report = pd.DataFrame([])
    for i in filenames_unique:
        file_one = i
        filenames.remove(file_one)
        basename = os.path.basename(os.path.splitext(file_one)[0])
        basename = basename[:len(basename) - 2]
        if paired_end_sequencing == True:
            file_two = ""
            for filename in filenames:
                if basename in filename:
                    file_two = filename
                    break
            if file_two == "":
                print("you are missing the the reverse strand for " + basename)
            fastq_files = file_one + " " + file_two
        else:
            fastq_files = file_one
        os.mkdir("my_experiments/" + experiment)
        os.mkdir("my_experiments/" + experiment + "/alignment_reports") # raise error if already exists
        filename_base = basename

        result = "temp/" + filename_base + ".vdjca"
        subprocess.run(["java",
                        "-jar",
                        path_to_mixcr,
                        "align",
                        "-p milab-human-tcr-dna-multiplex-cdr3",
                        fastq_files,
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
        clones_sample = clones_sample.rename(columns={"readFraction_x": "clonesFraction"})
        sequencing_report = pd.concat([sequencing_report, clones_sample])
        files_to_remove = os.listdir("temp")
        for file in files_to_remove:
            os.remove("temp/" + file)

    sequencing_report.to_csv("my_experiments/" + experiment + "/sequencing_report.txt")
    data["last_experiment"] = path_to_mixcr
    f = open(ROOT + 'global_vars.txt' 'w')
    f.write(str(experiment))