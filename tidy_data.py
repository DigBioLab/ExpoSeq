import pandas as pd
from tkinter import filedialog
from glob import glob
from itertools import groupby
from os import path


def collectFiles():
    filenames = []
    basenames = []
    ask_input = 'Y'
    while ask_input != 'N' or 'n' or 'no' or 'No':
        data_folder = filedialog.askdirectory()
        data_folder = path.abspath(data_folder) + '\\*'
        filenames_raw_data = glob(data_folder ) # might not be the same name for everyone, has to be changed if package should be used by everyone

        filenames.append(filenames_raw_data)
        ask_input = input("Do you want to add more folders? Press: (Y/n)")

    filenames = filenames.sort()
    return filenames


def group_files(filenames):
    grouped_filenames = [list(group) for key, group in groupby(
                                        filenames,
                                        lambda a: a.split('.')[0])
                         ]
    return grouped_filenames



#next steps: ideas:
# look for txt file in inner list
# extract sample name and save it in variable
# read txt file in pandas data table
# for the future: identify key variables for plotting and data analysis and tell whether they are present
# read other AlignmentReport and extract key variables and attach to table for each row
#



summed_evalues = []
for file in filenames_raw_data:
    Export = pd.read_tabel(file)
    Export2 = Export[(Export['lengthOfCDR3'] % 3) == 0] #This removes rows with a nt length not divisible by 3 (think aa translation)
    evalues = (str(Export2.cloneCount.sum())) #\n = go to the next line
    summed_evalues.append([file[0:file.index('_UniqueCDR3_Exp_UniqueCDR3.txt')], evalues])

