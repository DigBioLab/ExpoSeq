from os import path
import tkinter as tk
from tkinter import filedialog
from glob import glob



def collectFiles(filetypes = ('.txt', '.fastq_AlignmentReport')):
    filenames = []
    data_folder = filedialog.askdirectory()
    data_folder = path.abspath(data_folder) + '\\*'
        # maybe create user input who can choose which filetypes he/she wants to include
    for files in filetypes:
        filenames.extend(glob(data_folder+files)) # might not be the same name for everyone, has to be changed if package should be used by everyone
    filenames.sort(key = len)
    return filenames





