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


from glob import glob
import subprocess
from ast import literal_eval
import pickle
import pandas as pd
from ExpoSeq.augment_data.trimming import trimming
import pkg_resources
import editdistance
import os
import platform
#psutil
try:
    import tkinter as tk
    from tkinter import filedialog
except:
    pass

class CollectFastq():
    def __init__(self, paired_end_sequencing):
        self.forward = []
        self.backward = []
        self.paired = []
        self.paired_end_sequencing = paired_end_sequencing
    def add_fastq_files(self):
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
                user_choice = input(
                    "You have not chosen a directory with fastq files.\n If you want to try again press 1.\n If you want to cancel the analysis press 2")
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


    def get_filenames(self):
        while True:
            filenames = self.add_fastq_files()
            if len(filenames) == 0:
                print("You have not added any fastq files the preprocessing is cancelled now.")
            else:
                break
        return filenames
    def get_filename(self):
        try:
            while True:
                path_to_file = filedialog.askopenfilename()
                if os.path.isfile(path_to_file) == True:
                    break
                else:
                    pass

        except:
            while True:
                path_to_file = input("copy and paste the path to your directory")
                if os.path.isfile(os.path.abspath(path_to_file)) == True:
                    break
                else:
                    print("Please enter a valid directory. Dont add a \ to the end of the directory path")
        return path_to_file
    def __len__(self, file_list):
        return len(file_list)
    def get_pairs(self):
        if self.__len__(self.forward) != self.__len__(self.backward):
            print("Files matching aborted. Your have not collected the same number of files for the forward and backward reads.")
        else:
            max_similarity = 0
            best_pair = (None, None)
            # Pairwise comparison of filenames
            for i, file1 in enumerate(self.forward):
                for j, file2 in enumerate(self.backward):
                    similarity = 1 - editdistance.distance(file1, file2) / max(len(file1), len(file2))
                    if similarity > max_similarity:
                        max_similarity = similarity
                        best_pair = [file1, file2]
                self.paired.append(best_pair)
                
    def get_files(self, cmd = False, path_to_files = None, path_to_backward = None):
        if cmd == False:
            print("Choose the directory where you store the fastq files with the forward reads or single end sequencing data. \nIf you want to continue with paired end sequencing data make sure that you store your reverse reads in a seperate folder. \nFurther make sure your chosen directory does not contain fastq files from other experiments.")
            self.forward = self.get_filenames()
            if self.paired_end_sequencing:
                print("Now choose the directory where you store the fastq files with the backward reads.")
                self.backward = self.get_filenames()
                self.get_pairs()
        else:
            path_to_forward= os.path.abspath(path_to_files)
            self.forward.extend(glob(os.path.join(path_to_forward, "*.fastq")))
            if self.paired_end_sequencing:
                path_to_backward = os.path.abspath(path_to_backward)
                self.backward.extend(glob(os.path.join(path_to_backward, "*.fastq")))
                self.get_pairs()
            
    def call_pairs(self):
         if self.paired_end_sequencing:
             for i in self.paired:
                 print("You will process " + i[0] + " with " + i[1])

    def manually_match_pairs(self):
        self.paired = []
        self.backward = []
        self.forward = []
        while True:
            print("Choose the file of your fastq file with the forward reads.")
            forward_file = self.get_filename()
            self.forward.append(forward_file)
            print("Choose the file of your fastq file with the backward reads")
            backward_file = self.get_filename()
            self.backward.append(backward_file)
            self.paired.append([forward_file, backward_file])
            read_more_files = input("do you want to add more files? (Y/n)")
            if read_more_files in ["Y" "y", "n", "N"]:
                break
            else:
                pass



class CreateCommand:
    def __init__(self, module_dir, path_to_mixcr, paired_end_sequencing, experiment, files):
        self.module_dir = module_dir
        self.path_to_mixcr = path_to_mixcr
        self.paired_end_sequencing = paired_end_sequencing
        basename = os.path.basename(os.path.splitext(files[0])[0])
        self.basename = basename[:len(basename) - 2]
        self.filename_base = basename
        self.result = os.path.join(self.module_dir, "temp", self.filename_base + ".vdjca")
        self.alignment_path = os.path.normpath(os.path.join(self.module_dir,"my_experiments",experiment,"alignment_reports",self.filename_base + "_AlignmentReport.txt"))
        self.clones = os.path.join(self.module_dir, "temp", self.basename + "clones.clns")
        self.table_tsv = os.path.join(module_dir,
                                     "temp",
                                     basename + ".tsv")
        self.files = [os.path.normpath(i) for i in files]

    def get_free_memory_linux(self):
        with open("/proc/meminfo", "r") as f:
            lines = f.readlines()
        for line in lines:
            if line.startswith("MemAvailable"):
                return int(line.split()[1]) // 1024  # Convert from KB to MB

    def get_free_memory_windows(self):
        try:
            import psutil
            return psutil.virtual_memory().available // (1024 ** 2)  # Convert from Bytes to MB
        except ImportError:
            print("Please install psutil library: pip install psutil")


    def get_free_memory_mac(self):
        cmd = "vm_stat | grep 'Pages free' | awk '{print $3}'"
        pages_free = int(os.popen(cmd).read().replace(".", ""))
        # Typically, a page in macOS is 4096 bytes
        return (pages_free * 4096) // (1024 ** 2)  # Convert from Bytes to MB

    def get_free_memory(self):
        system_platform = platform.system()
        if system_platform == "Linux":
            return self.get_free_memory_linux()
        elif system_platform == "Windows":
            return self.get_free_memory_windows()
        elif system_platform == "Darwin":  # This is for Mac OS
            return self.get_free_memory_mac()
        else:
            print("Unsupported platform")
            exit(1)

    def create_parser(self, java_heap_size = None):
        if java_heap_size is None:
            try:
                free_memory = self.get_free_memory()
                java_heap_size = free_memory // 2  # Using half of the available RAM as an example
                print(f"50 % of currently available memory: {java_heap_size}\n")
            except:
                print("automatic detection of memory failed. Processing continues with 1000MB RAM.")
                java_heap_size = 1000
        else:
            pass
        commands = []
        commands.extend(["java", f"-Xms{java_heap_size}M", "-jar"]) # enable change of para
        commands.extend([self.path_to_mixcr])
        return commands


    def prepare_align(self, method):
        align_commands = self.create_parser()
        align_commands.extend(["align"])
        align_commands.extend(["--preset", method])
        align_commands.extend(self.files)
        align_commands.extend([self.result])
        align_commands.extend(["--report", self.alignment_path])
        return align_commands

    def prepare_assembly(self):
        assembly_commands = self.create_parser()
        assembly_commands.extend(["assemble"])
       # assembly_commands.extend(["-OseparateByC=true", "-OseparateByV=true","-OseparateByJ=true"])
        assembly_commands.extend([self.result])
        assembly_commands.extend([self.clones])

        return assembly_commands


    def prepare_clones(self):
        clones_commands = self.create_parser()
        clones_commands.extend(["exportClones"])
        clones_commands.extend([self.clones])
        clones_commands.extend([self.table_tsv])
        return clones_commands


def mixcr_(fastq_directory, path_to_mixcr, path_to_forward, path_to_backward, paired_end_sequencing, threads, method, java_heap_size = None):
    print(path_to_mixcr)
    assert os.path.isdir(fastq_directory), "The given directory for the fastq files does not exist."
    assert isinstance(threads, int), "The given threads is not an integer."    
    assert subprocess.run(["java", "-jar",path_to_mixcr,"--version"]), "The given path to mixcr is not correct. Or you have not installed mixcr correctly."
    

    
    confirmation = True
    Collect = CollectFastq(paired_end_sequencing)
    module_dir = Collect.module_dir
    Collect.get_files(cmd = True, path_to_files = path_to_forward, path_to_backward = path_to_backward)
    files = Collect.paired
    if paired_end_sequencing:
        print("Found pairs for the forward and reverse files")
        Collect.call_pairs()
    if not os.path.isdir(os.path.join(module_dir, "alignment_reports")):
        os.mkdir(os.path.join(module_dir, "alignment_reports"))
            
    sequencing_report = pd.DataFrame([])
    for filename in files:
        Commands = CreateCommand(module_dir, path_to_mixcr, paired_end_sequencing, experiment, filename)
        subprocess.run(Commands.prepare_align(method, java_heap_size))
        subprocess.run(Commands.prepare_assembly())
        subprocess.run(Commands.prepare_clones())
        basename = Commands.basename
        clones_sample = pd.read_table(os.path.join(module_dir,
                                                   "temp",
                                                   basename + "_IGH.tsv"))
        sequencing_report = pd.concat([sequencing_report, clones_sample])
        files_to_remove = os.listdir(os.path.join(module_dir, "temp"))
        for file in files_to_remove:
            os.remove(os.path.join(module_dir,
                                   "temp",
                                   file))
    sequencing_report.to_csv('sequencing_report.csv')

parser = argparse.ArgumentParser(description='Processing')

# Add arguments to the parser
parser.add_argument('fastq_directory',type = str, help='Path to the directroy containg the fastq files')
parser.add_argument('path_to_mixcr', type=float, help="Path to the mixcr jar file")

parser.add_argument('--paired_end_sequencing',default=False, default=32, help='sequencing technique you used')
parser.add_argument('--threads', default=1, type=int, help='Number of threads to use')
parser.add_argument("--method", default="milab-human-tcr-dna-multiplex-cdr3", type=str, help="Method to use for the alignment")
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
