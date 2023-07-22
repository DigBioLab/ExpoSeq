#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import os
import subprocess
import pickle
from glob import glob
import shutil
from ExpoSeq.augment_data.trimming import trimming
import sys
import argparse
#from .trimming import trimming

#
def remove_directories(directory_path):
    # Iterate over all items in the directory
    for item in os.listdir(directory_path):
        shutil.rmtree(os.path.join(directory_path, item))


def collect_fastq_files(fastq_directory):
    filenames = []
    path_to_files = os.path.abspath(fastq_directory)
    filenames.extend(glob(os.path.join(path_to_files, "*.fastq")))
    if len(filenames) == 0:
       print(
            "You have not chosen a directory with fastq files.")

    else:
        print("These are the files you chose:")
        for i in filenames:
            print(os.path.basename(i))

    return filenames


def mixcr_(fastq_directory, path_to_mixcr, save_dir, paired_end_sequencing, threads, method, trim_div_by, trim_min_count):
    print(path_to_mixcr)
    assert os.path.isdir(fastq_directory), "The given directory for the fastq files does not exist."
    #assert os.path.isfile(path_to_mixcr), "The given path to the mixcr jar file is not correct."
    #assert os.path.isdir(save_dir), "The given directory for the save directory does not exist."
    #assert isinstance(paired_end_sequencing, bool), "The given paired_end_sequencing is not a boolean."
    assert isinstance(threads, int), "The given threads is not an integer."
    assert trim_div_by > 0, "The given trim_div_by is not greater than 0."
    assert trim_min_count > 0, "The given trim_min_count is not greater than 0."


    filenames = collect_fastq_files(fastq_directory)

    if len(filenames) == 0:
        print("You have not added any fastq files.")
        return
    print("Test Mixcr")
    subprocess.run(["java",
                    "-jar",
                    path_to_mixcr,
                    "--version"
                    ])
    confirmation = True
    if confirmation == True:
        if not os.path.isdir("temp"):
            os.mkdir("temp")

        module_dir = os.path.abspath("")
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
            try:
                basename = os.path.basename(os.path.splitext(filename)[0])
                basename = basename[:len(basename) - 2]

                filename_base = basename
                alignment_report_path = os.path.normpath(os.path.join(module_dir, filename_base + "_AlignmentReport.txt"))
                tsv_path = os.path.join(module_dir,
                                "temp",
                                basename + ".tsv")
                result = os.path.join(module_dir, "temp", filename_base + ".vdjca")
                clones = os.path.join(module_dir, "temp", basename + "clones.clns")
                subprocess.run(["java",
                                "-jar",
                                path_to_mixcr,
                                "align",
                                "-p " + method,
                                filename,
                                result,
                                "--report",
                                alignment_report_path,
                                "--threads",
                                str(threads)
                                ])

                subprocess.run(["java",
                                "-jar",
                                path_to_mixcr,
                                "assemble",
                                "-OseparateByC=true",
                                "-OseparateByV=true",
                                "-OseparateByJ=true",
                                result,
                                clones,
                                "--threads",
                                str(threads)
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
                                tsv_path,
                                "--threads",
                                str(threads)
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
                                                    how="left",
                                                    on="nSeqCDR3")
                clones_sample = clones_sample.sort_values(by="cloneId")
                clones_sample = clones_sample.reset_index()
                clones_sample = clones_sample.drop(columns=["readFraction_y", "index"])
                clones_sample = clones_sample.rename(columns={"readFraction_x": "cloneFraction"})
                clones_sample["Experiment"] = filename_base
                clones_sample = trimming(clones_sample,
                                         divisible_by=trim_div_by,
                                         min_count=trim_min_count,
                                         new_fraction="clonesFraction")
                sequencing_report = pd.concat([sequencing_report, clones_sample])
                files_to_remove = os.listdir(os.path.join(module_dir, "temp"))
                for file in files_to_remove:
                    os.remove(os.path.join(module_dir,
                                           "temp",
                                           file))
                report_dir = os.path.join(save_dir,
                                      "sequencing_report.txt")
                sequencing_report.to_csv(report_dir)
            except:
                remove_directories(os.path.join(module_dir, "temp"))

    else:
        print("The path to the mixcr jar file is not correct or you are missing the licence")
        



parser = argparse.ArgumentParser(description='Processing')

# Add arguments to the parser
parser.add_argument('fastq_directory',type = str, help='Path to the directroy containg the fastq files')
parser.add_argument('path_to_mixcr', type=float, help="Path to the mixcr jar file")
parser.add_argument('--save_dir', default='', help='Path to the directory wherethe sequencing report should be saved')
parser.add_argument('--paired_end_sequencing',default=False, default=32, help='sequencing technique you used')
parser.add_argument('--threads', default=1, type=int, help='Number of threads to use')
parser.add_argument("--method", default="milab-human-tcr-dna-multiplex-cdr3", type=str, help="Method to use for the alignment")
parser.add_argument("--trim_div_by", default=3, type=int, help="Trim all sequences that are not divisible by this number")
parser.add_argument("--trim_min_count", default=3, type=int, help="Trim all sequences that have a count lower than this number")
args = parser.parse_args()

fastq_directory = args.fastq_directory
path_to_mixcr = args.path_to_mixcr
save_dir = args.save_dir
paired_end_sequencing = args.paired_end_sequencing
threads = args.threads
method = args.method
trim_div_by = args.trim_div_by
trim_min_count = args.trim_min_count

# Call the function
mixcr_(fastq_directory, path_to_mixcr, save_dir, paired_end_sequencing, threads, method, trim_div_by, trim_min_count)
