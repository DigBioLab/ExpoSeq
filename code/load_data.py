from os import path
from tkinter import filedialog
from glob import glob


def collectFiles():
    filenames = []

    ask_input = 'Y'
    while ask_input != 'n':
        data_folder = filedialog.askdirectory()
        data_folder = path.abspath(data_folder) + '\\*'
        filetypes = ('.txt', '.clns', '.fastq', '.fastq_AlignmentReport') # maybe create user input who can choose which filetypes he/she wants to include
        for files in filetypes:
            filenames.extend(glob(data_folder+files)) # might not be the same name for everyone, has to be changed if package should be used by everyone
        ask_input = input("Do you want to add more folders? Press: (Y/n)")
    filenames.sort(key = len)
    return filenames
