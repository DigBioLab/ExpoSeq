import os.path
from matplotlib.pyplot import savefig
from ast import literal_eval
import pkg_resources
try:
    import tkinter as tk
    from tkinter import filedialog
except:
    pass

def saveFig(name):
    pkg_path = pkg_resources.resource_filename("ExpoSeq", "")
    save_settings_file = os.path.join(pkg_path,
                                      "settings",
                                      "save_settings.txt")
    with open(save_settings_file, "r") as f:
        save_settings = f.read()
    save_settings = literal_eval(save_settings)
    if name == False:
        filename = input("how do you want to call the plot?")
    else:
        filename = name
    if save_settings["fname"] == "":
    #    print("Please choose a directory where you want to store your files")
        try:
            print("Choose the directory where you want to save the file with the filechooser")
            save_dir = filedialog.askdirectory()
        except:
            while True:
                save_dir = input("copy and paste the path to the directory where you want to save the files here.")
                if os.path.isdir(os.path.abspath(save_dir)):
                    break
                else:
                    print("Try again and input a valid directory")
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