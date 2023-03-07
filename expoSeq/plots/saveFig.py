import os.path
from matplotlib.pyplot import savefig

from ast import literal_eval
import tkinter as tk
from tkinter import filedialog
import time


def saveFig():
    module_dir = os.path.abspath("expoSeq")
    save_settings_file = os.path.join(module_dir, "settings", "save_settings.txt")
    with open(save_settings_file, "r") as f:
        save_settings = f.read()
    save_settings = literal_eval(save_settings)
    filename = input("how do you want to call the plot?")
    if save_settings["fname"] == "":
    #    print("Please choose a directory where you want to store your files")
        save_dir = filedialog.askdirectory()
        save_dir = os.path.abspath(save_dir)
        save_settings["fname"] = save_dir
        with open(save_settings_file, "w") as f:
            f.write(str(save_settings))
    else:
        pass
    filepath = save_settings['fname']
    form = save_settings['format']
    path = os.path.join(filepath, filename + "." + form)
    del save_settings['fname']
    savefig(path,
            **save_settings)