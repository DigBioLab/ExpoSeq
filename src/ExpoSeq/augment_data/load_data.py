import os.path
from os import path

from glob import glob
try:
    import tkinter as tk
    from tkinter import filedialog
except:
    pass
def collectFiles(filetypes = ('.txt', ".tsv", '.fastq_AlignmentReport')):
    filenames = []
    try:
        data_folder = filedialog.askdirectory()
    except:
        while True:
            data_folder = input("copy and paste the path to the directory without any / as ending.")
            if path.isdir(os.path.abspath(data_folder)) == True:
                break
    data_folder = path.abspath(data_folder) + '\\*'
        # maybe create user input who can choose which filetypes he/she wants to include
    for files in filetypes:
        filenames.extend(glob(data_folder+files)) # might not be the same name for everyone, has to be changed if package should be used by everyone
    filenames.sort(key = len)
    return filenames





