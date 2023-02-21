import os.path

from matplotlib.pyplot import savefig
from tkinter.filedialog import askdirectory
from ast import literal_eval

def saveFig():

    with open(os.path.abspath("settings/save_settings.txt"), "r") as f:
        save_settings = f.read()
    save_settings = literal_eval(save_settings)
    filename = input("how do you want to call the plot?")
    if save_settings["fname"] == "":
        print("Please choose a directory where you want to store your files")
        directory = askdirectory()
        save_settings["fname"] = directory
        with open(os.path.abspath("settings/save_settings.txt"), "w") as f:
            f.write(str(save_settings))
    else: pass
    filepath = save_settings['fname']
    format = save_settings['format']
    path = os.path.join(filepath, filename + "." + format)
    del save_settings['fname']
    savefig(path, **save_settings)