import os
from ast import literal_eval

import pkg_resources
try:
    import tkinter as tk
    from tkinter import filedialog
except:
    pass
class Change_save_settings():
    def __init__(self):
        self.pkg_path = pkg_resources.resource_filename("ExpoSeq", "")

        self.save_settings_path = os.path.join(self.pkg_path,
                                          "settings",
                                        "save_settings.txt")
        with open(self.save_settings_path, "r") as f:
            save_settings = f.read()
        self.save_settings = literal_eval(save_settings)
    def change_dpi(self, dpi):
        self.save_settings["dpi"] = dpi
        with open(self.save_settings_path, "w") as f:
            f.write(str(self.save_settings))
    def change_save_path(self, path = None):
        if path == None:
            try:
                self.save_settings["fname"] = filedialog.askopenfilename()
            except:
                while True:
                    save_dir = input("copy and paste the path to the directory where you want to save your plots")
                    if os.path.isdir(os.path.abspath(save_dir)):
                        self.save_settings["fname"] = save_dir
                        break
                    else:
                        print("Try again and input a valid directory")
        else:
            self.save_settings["fname"] = path
            print("The new path you safe your results in is: " + path)
        with open(self.save_settings_path, "w") as f:
            f.write(str(self.save_settings))

    def change_format(self, format):
        self.save_settings["format"] = format
        with open(self.save_settings_path, "w") as f:
            f.write(str(self.save_settings))

    def print_save_path(self):
        print(self.save_settings["fname"])

