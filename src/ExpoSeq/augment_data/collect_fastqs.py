from tkinter import filedialog
import os
from glob import glob
import gzip


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
            filenames.extend(glob(os.path.join(path_to_files, "*.fastq*")))
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
                    if file1.endswith(".gz"):
                        with gzip.open(file1, 'rb') as f:
                            one_first = f.readline().strip().decode('utf-8')
                            one_sub_first = one_first.rfind(":")
                            one_first_sub = one_first[one_sub_first + 1:]
                    else:
                        with open(file1, "r") as one:
                            one_first = one.readline().strip()
                            one_sub_first = one_first.rfind(":")
                            one_first_sub = one_first[one_sub_first + 1:]
                    if file2.endswith(".gz"):
                        with gzip.open(file2, 'rb') as f:
                            two_first = f.readline().strip().decode('utf-8')
                            two_sub_first = two_first.rfind(":")
                            two_first_sub = two_first[two_sub_first + 1:]
                    else:
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