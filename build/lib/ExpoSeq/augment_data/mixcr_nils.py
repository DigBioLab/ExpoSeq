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
import time
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
    def get_files(self):
        print("Choose the directory where you store the fastq files with the forward reads or single end sequencing data. \nIf you want to continue with paired end sequencing data make sure that you store your reverse reads in a seperate folder. \nFurther make sure your chosen directory does not contain fastq files from other experiments.")
        self.forward = self.get_filenames()
        if self.paired_end_sequencing:
            print("Now choose the directory where you store the fastq files with the backward reads.")
            self.backward = self.get_filenames()
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



class CreateCommand:
    def __init__(self, module_dir, path_to_mixcr, paired_end_sequencing, experiment, files, java_heap_size):
        self.java_heap_size = java_heap_size
        self.module_dir = module_dir
        self.path_to_mixcr = path_to_mixcr
        self.paired_end_sequencing = paired_end_sequencing


        if paired_end_sequencing:
            basename = os.path.basename(os.path.splitext(files[0][0])[0])
        else:
            print(files[0])
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
        clones_commands.extend(["-c IGH"])
        clones_commands.extend([self.clones])
        clones_commands.extend([self.table_tsv])
        return clones_commands

def process_mixcr(experiment, method, testing, paired_end_sequencing):
    pkg_path = pkg_resources.resource_filename("ExpoSeq", "")
    module_dir = os.path.abspath("")

    if not testing:
        Collect = CollectFastq(paired_end_sequencing)
        Collect.get_files()
        files = Collect.paired
        files = [[i] for i in files]
        if paired_end_sequencing:
            while True:
                if len(files) > 0:
                    Collect.call_pairs()
                    correct_matches = input("Are the matches of forward and backward reads correct? Y/n")
                else:
                    correct_matches = "n"
                    print("You will need to restart the collection of the files if you want to proceed with paired end sequencing data.\nMake sure that you have the forward and backward files in two seperate folder and that you have the same number of files.")
                if correct_matches in ["Y", "y", "N", "n"]:
                    if correct_matches in ["Y", "y"]:
                        break
                    else:
                        while True:
                            choose_strat = input("You can restart the collection of the files or manually add them. \nPress 1 for restarting the collection by choosing the directories.\nPress 2 for restarting the collection and manually choose the forward and backward reads.")
                            if choose_strat in ["1", "2"]:
                                Collect = CollectFastq(paired_end_sequencing)
                                if choose_strat == "1":
                                    Collect.get_files()
                                    files = Collect.paired
                                    break
                                else:
                                    Collect.manually_match_pairs()
                                    files = Collect.paired
                                    break
                            else:
                                print("Please enter a valid value")
                else:
                    print("Please enter a valid value")
        else:
            pass
        if not os.path.isdir(os.path.join(module_dir, "my_experiments", experiment, "alignment_reports")):
            os.mkdir(os.path.join(module_dir, "my_experiments", experiment, "alignment_reports"))
    else:
        files = glob(os.path.join(pkg_path, "test_data", "test_files", "*.fastq"))
        print(files)
        files = [[i] for i in files]
        experiment = "test_directory"
    settings_dir = os.path.join(module_dir,
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

    sequencing_report = pd.DataFrame([])
    while True:
        java_heap_size = input(
            "Enter the amount of memory you want to allocate for java in Mb. Normally 1000 is sufficient.")
        try:
            java_heap_size = int(java_heap_size)
            break
        except:
            print("Please enter the amount as integer")
    files_to_remove = os.listdir(os.path.join(module_dir, "temp"))
    for filename in files:
        for file in files_to_remove:
            os.remove(os.path.join(module_dir,
                                   "temp",
                                   file))
        Commands = CreateCommand(module_dir, path_to_mixcr, paired_end_sequencing, experiment, filename, java_heap_size)
        subprocess.run(Commands.prepare_align(method))
        subprocess.run(Commands.prepare_assembly())
        subprocess.run(Commands.prepare_clones())
        basename = Commands.basename

        for filename in os.listdir(os.path.join(module_dir, "temp")):
            if filename.endswith(".tsv"):
                clones_sample = pd.read_table(os.path.join(module_dir, "temp", filename))
        clones_sample["Experiment"] = basename
        sequencing_report = pd.concat([sequencing_report, clones_sample])

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
        glob_vars_dir = os.path.join(module_dir,
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
    print(
        "The output of your processing is a sequencing_report which contains the sequences and their sample names of all fastq files. Further, a folder with the alignment reports for each fastq file was created. You can have a look under ~/my_experiments/YOUR_EXPERIMENT_NAME")

