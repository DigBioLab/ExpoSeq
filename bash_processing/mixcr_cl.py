

import subprocess
import pandas as pd
import pickle
import argparse
from glob import glob
from tkinter import filedialog
#from ExpoSeq.pipeline import PlotManager
import os
from tkinter import filedialog

from glob import glob
import time


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
            best_pair = (None, None)
            # Pairwise comparison of filenames
            for i, file1 in enumerate(self.forward):
                for j, file2 in enumerate(self.backward):
                    with open(file1, "r") as one:
                        one_first = one.readline().strip()
                        one_sub_first = one_first.rfind(":")
                        one_first_sub = one_first[one_sub_first+1:]
                    with open(file2, "r") as two:
                        two_first = two.readline().strip()
                        two_sub_first = two_first.rfind(":")
                        two_first_sub = two_first[two_sub_first+1:]
                    if one_first_sub == two_first_sub:
                        best_pair = [file1, file2]

                if not best_pair:
                    print(f"Could not find match for {file1}")
                else:
                    self.paired.append(best_pair)
    def get_files(self,  cmd = False, path_to_forward = None, path_to_backward = None):
        if cmd == False: 
            print("Choose the directory where you store the fastq files with the forward reads or single end sequencing data. \nIf you want to continue with paired end sequencing data make sure that you store your reverse reads in a seperate folder. \nFurther make sure your chosen directory does not contain fastq files from other experiments.")
            self.forward = self.get_filenames()
            if self.paired_end_sequencing:
                print("Now choose the directory where you store the fastq files with the backward reads.")
                self.backward = self.get_filenames()
                self.get_pairs()
            else:
                self.paired = self.forward
        else:
            path_to_forward= os.path.abspath(path_to_forward)
            self.forward.extend(glob(os.path.join(path_to_forward, "*.fastq*")))
            if self.paired_end_sequencing == True:
                path_to_backward = os.path.abspath(path_to_backward)
                self.backward.extend(glob(os.path.join(path_to_backward, "*.fastq*")))
                self.get_pairs()
            else:
                self.paired = self.forward

                
 

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

def check_dirs(module_dir, experiment):
    if not os.path.isdir(os.path.join(module_dir, "my_experiments")):
        os.mkdir(os.path.join(module_dir, "my_experiments"))
    if not os.path.isdir(os.path.join(module_dir, "my_experiments", experiment)):
        os.mkdir(os.path.join(module_dir, "my_experiments", experiment))      
    if not os.path.isdir(os.path.join(module_dir, "my_experiments", experiment, "alignment_reports")):
        os.mkdir(os.path.join(module_dir, "my_experiments", experiment, "alignment_reports"))
    if not os.path.isdir(os.path.join(module_dir, "temp")):
        os.mkdir(os.path.join(module_dir, "temp"))   
    if not os.path.isdir(os.path.join(module_dir, "my_experiments", experiment, "assembly_reports")):
        os.mkdir(os.path.join(module_dir, "my_experiments", experiment, "assembly_reports"))
    if not os.path.isdir(os.path.join(module_dir, "my_experiments", experiment, "clones_result")):
        os.mkdir(os.path.join(module_dir, "my_experiments", experiment, "clones_result"))
    if not os.path.isdir(os.path.join(module_dir,
                                     "my_experiments",
                                     experiment,
                                     "tables_mixcr")):
        os.mkdir(os.path.join(module_dir, "my_experiments", experiment, "tables_mixcr"))



class CreateCommand:
    def __init__(self, module_dir, path_to_mixcr, paired_end_sequencing, experiment, files, java_heap_size, threads = None):
        self.java_heap_size = java_heap_size
        self.module_dir = module_dir
        self.path_to_mixcr = path_to_mixcr
        if threads != None:
            if type(threads) == str:
                self.threads = threads
            elif type(threads) == int:
                self.threads = str(threads)
        self.paired_end_sequencing = paired_end_sequencing
        if paired_end_sequencing:
            self.files = [os.path.normpath(j) for i in files for j in i]
            basename = os.path.basename(files[0][0]).split(".")[0]
        else:
            self.files = files
            basename = os.path.basename(files[0]).split(".")[0]
        self.basename = basename
        self.result = os.path.join(self.module_dir,
                                   "temp",
                                   self.basename + ".vdjca")
        self.alignment_path = os.path.normpath(os.path.join(self.module_dir,
                                                            "my_experiments",
                                                            experiment,
                                                            "alignment_reports",
                                                            self.basename + "_AlignmentReport.txt"))
        self.assembly_path = os.path.normpath(os.path.join(self.module_dir,
                                                           "my_experiments",
                                                           experiment,
                                                           "assembly_reports",
                                                           self.basename + "_AssemblyReport.txt"))
        self.clones = os.path.join(self.module_dir,
                                   "my_experiments",
                                   experiment,
                                   "clones_result",
                                   self.basename + "_clones.clns")
        self.clones_base = os.path.dirname(self.clones)
        self.table_tsv = os.path.join(module_dir,
                                     "my_experiments",
                                     experiment,
                                     "tables_mixcr",
                                     self.basename + ".tsv")
        self.mixcr_plots_path = os.path.normpath(os.path.join(self.module_dir,
                                                              "my_experiments",
                                                              experiment,
                                                              "mixcr_plots"))

    def read_aligned_reads(self):
        with open(self.alignment_path, "r") as f:
            alignment = f.readlines()
        for i in alignment:
            if "Successfully aligned reads" in i:
                quality = i.replace("Successfully aligned reads: ", "")
        align_quality = quality.split(" ")[0]
        percent_quality = quality.split(" ")[1]
        align_quality = int(align_quality)
        return align_quality, percent_quality

    def read_clones(self):
        with open(self.assembly_path, "r") as f:
            alignment = f.readlines()
        for i in alignment:
            if "Reads used in clonotypes, percent of total" in i:
                quality = i.replace("Reads used in clonotypes, percent of total: ", "")
        align_quality = quality.split(" ")[0]
        percent_quality = quality.split(" ")[1]
        assembly_quality = int(align_quality)
        return assembly_quality, percent_quality

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



    def create_parser(self):
     #   try:
      #      free_memory = self.get_free_memory()
      #      java_heap_size = int(free_memory / 4)  # Using half of the available RAM as an example
       #     print(f"25 % of currently available memory: {java_heap_size}\n")
        #except:
         #   print("automatic detection of memory failed. Processing continues with 1000MB RAM.")
          #  java_heap_size = 1000
        commands = []
        commands.extend(["java", f"-Xms{self.java_heap_size}M", "-jar"]) # enable change of para
        commands.extend([self.path_to_mixcr])
        return commands


    def prepare_align(self, method):
        align_commands = self.create_parser()
        align_commands.extend(["align"])
        print(method)
        align_commands.extend(["--preset", method])
        align_commands.extend(self.files)
        align_commands.extend([self.result])
        align_commands.extend(["--no-warnings"])
        align_commands.extend(["--force-overwrite"])
        if self.threads != None:
            align_commands.extend(["--threads", self.threads])
        align_commands.extend(["--report", self.alignment_path])
        return align_commands

    def prepare_assembly(self):

        assembly_commands = self.create_parser()
        assembly_commands.extend(["assemble"])
       # assembly_commands.extend(["-OseparateByC=true", "-OseparateByV=true","-OseparateByJ=true"])
        assembly_commands.extend([self.result])
        assembly_commands.extend([self.clones])
        assembly_commands.extend(["--force-overwrite"])

        assembly_commands.extend(["--report", self.assembly_path])

        return assembly_commands


    def prepare_clones(self):
        clones_commands = self.create_parser()
        clones_commands.extend(["exportClones"])
        clones_commands.extend(["-c IGH"])
        clones_commands.extend([self.clones])
        clones_commands.extend(["--force-overwrite"])
        clones_commands.extend([self.table_tsv])
        return clones_commands


 
def mixcr_(path_to_mixcr, experiment_name,  path_to_forward, path_to_backward = None, method = "ampliseq-tcrb-plus-full-length", threads = 1, java_heap_size = None):
    print(path_to_backward)
    print(path_to_forward)
    print(experiment_name)
    print(method)
    print(threads)
    print(java_heap_size    )
    module_dir = os.path.abspath("")
    assert os.path.isfile(path_to_mixcr), "File path to mixcr.jar file is incorrect."
    assert path_to_forward != path_to_backward, "Directory to forward and backward files needs to be different"
    assert os.path.isdir(path_to_forward), "The given directory for the forward files does not exist."
    assert subprocess.run(["java", "-jar",path_to_mixcr,"--version"]), "The given path to mixcr is not correct. Or you have not installed mixcr correctly."
    assert type(experiment_name) == str, "Please enter a string as experiment name"
  #  assert os.path.isfile(module_dir, "my_experiments", experiment_name) == False, "The experiment name already exists, please enter another one."
    check_dirs(module_dir, experiment_name)
    if java_heap_size == None:
        java_heap_size = 1000
    if path_to_backward == None:
        paired_end_sequencing = False
    else:
        paired_end_sequencing = True

    Collect = CollectFastq(paired_end_sequencing)
    
    Collect.module_dir = module_dir

    Collect.get_files(cmd = True,
                      path_to_forward = path_to_forward,
                      path_to_backward = path_to_backward)
    files = Collect.paired
    if paired_end_sequencing:
        print("Found pairs for the forward and reverse files")
        Collect.call_pairs()

            
    sequencing_report = pd.DataFrame([])
    for filename in files:
        for file in os.listdir(os.path.join(module_dir, "temp")):
            os.remove(os.path.join(module_dir,
                                        "temp",
                                        file))
        filename = [filename]

        Commands = CreateCommand(module_dir, path_to_mixcr, paired_end_sequencing, experiment_name, filename, java_heap_size, threads = threads)

        subprocess.run(Commands.prepare_align(method))
        skip_sample = False
        try:
            align_quality, percent_quality = Commands.read_aligned_reads()
            if align_quality < 1000 and align_quality >10:
                print(f"Only {align_quality} reads could be aligned ({percent_quality}).")
            if align_quality < 10:
                print("Too less reads were aligned. Please choose the correct mixcr method.")
                skip_sample = True
            if skip_sample == False:
                subprocess.run(Commands.prepare_assembly())
                assembly_qual, realtive_assembly_quality = Commands.read_clones()
                if assembly_qual < 1000 and realtive_assembly_quality > 10:
                    print(f"Only {assembly_qual} reads could be assembled ({realtive_assembly_quality}).")
                if assembly_qual < 10:
                    print("Too less reads were assembled. Please choose the correct ")
                    skip_sample = True
                if skip_sample == False:
                    subprocess.run(Commands.prepare_clones())
                    basename = Commands.basename
                    columns_not_wanted = ['refPoints','allVHitsWithScore', 'allDHitsWithScore',
                                            'allJHitsWithScore', 'allCHitsWithScore', 'allVAlignments',
                                            'allDAlignments', 'allJAlignments', 'allCAlignments']
                    clones_sample = pd.read_table(Commands.table_tsv)
                    clones_sample.drop(columns = columns_not_wanted, inplace = True)
                    clones_sample["Experiment"] = basename
                    sequencing_report = pd.concat([sequencing_report, clones_sample])

                    for file in os.listdir(os.path.join(module_dir, "temp")):
                        os.remove(os.path.join(module_dir,
                                            "temp",
                                            file))
        except:
            print(f"File {filename} could not be processed. Please check if you have chosen the correct mixcr method or verify that the sequencing was successful.")

        
     #   except:
      #      print(f"File {filename} could not be processed. Please check if you have chosen the correct mixcr method or verify that the sequencing was successful.")

    sequencing_report.to_csv(os.path.join(module_dir,"my_experiments", experiment_name, 'sequencing_report.csv'))
    if sequencing_report.shape[0] == 0:
        raise Exception(
            "The processing failed for all of your files. Most likely because you chose the wrong mixcr method.")

    unique_experiments = sequencing_report["Experiment"].unique()
    experiment_dic = {item: item for item in list(unique_experiments)}
    exp_names_dir = os.path.join(module_dir,
                                    "my_experiments",
                                    experiment_name,
                                    "experiment_names.pickle")
    with open(exp_names_dir, "wb") as f:
        pickle.dump(experiment_dic, f)
    print(
    "The output of your processing is a sequencing_report which contains the sequences and their sample names of all fastq files. Further, a folder with the alignment reports for each fastq file was created. You can have a look under ~/my_experiments/YOUR_EXPERIMENT_NAME")



parser = argparse.ArgumentParser(description='Processing')
parser.add_argument('path_to_mixcr', type=str, help="Path to the mixcr jar file")
parser.add_argument('experiment_name', type=str, help="Name of the experiment")
parser.add_argument('path_to_forward',type = str, help='Directory to the fastq files with the forward reads')
parser.add_argument('--path_to_backward',type = str,default=None, help='Directory to the fastq files with the backward reads')
parser.add_argument('--threads', default=1, type=int, help='Number of threads to use')
parser.add_argument("--method", default="milab-human-tcr-dna-multiplex-cdr3", type=str, help="Method to use for the alignment")
parser.add_argument("--java_heap_size", default="1000", type=int, help="Memory allocation for the processing in MB")

args = parser.parse_args()
path_to_mixcr = args.path_to_mixcr
experiment_name = args.experiment_name
path_to_forward = args.path_to_forward
path_to_backward = args.path_to_backward
threads = args.threads
method = args.method
java_heap_size = args.java_heap_size

start_time = time.time()
# Call the function
mixcr_(path_to_mixcr, experiment_name, path_to_forward, path_to_backward,  method,threads,  java_heap_size)
#PlotManager(experiment = experiment_name)
end_time = time.time()
execution_time = end_time - start_time
print(f"Execution Time: {execution_time} seconds")
print(f"Launch ExpoSeq with from ExpoSeq.pipeline import PlotManager\n plot = PlotManager()\n and press 2 when it asks you for your new experiment.\n Choose the folder in my_experiments/{experiment_name} which will then be copy and pasted in ExpoSeq's working directory.")
  