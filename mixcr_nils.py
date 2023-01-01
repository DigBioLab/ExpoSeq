import os
from glob import glob
from tkinter import filedialog
import subprocess
from os.path import dirname, abspath
from ast import literal_eval

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


def process_mixcr():
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

    paired_end_sequencing = True
    filenames_unique = [unique for unique in filenames if "_1.fastq" in unique]
    for i in filenames_unique:
        file_one = i
        filenames.remove(file_one)
        if paired_end_sequencing == True:
            file_two = ""
            basename = os.path.basename(os.path.splitext(file_one)[0])
            basename = basename[:len(basename)-2]
            for filename in filenames:
                if basename in filename:
                    file_two = filename
                    break
            if file_two == "":
                print("you are missing the the reverse strand for " + basename)
            fastq_files = file_one + " " + file_two
        else:
            fastq_files = file_one
        filename_base = basename
        reports = "alignment_reports/"
        result = "temp/" + filename_base + ".vdjca"
        subprocess.run(["java",
                        "-jar",
                        path_to_mixcr,
                        "align",
                        "-p milab-human-tcr-dna-multiplex-cdr3",
                        fastq_files,
                        result,
                        "-r alignment_reports/" + filename_base + "_AlignmentReport.txt"])

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
