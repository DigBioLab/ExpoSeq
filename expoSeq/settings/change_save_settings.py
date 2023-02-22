import os
from ast import literal_eval
from tkinter.filedialog import askopenfilename
class Change_save_settings():
    def __init__(self):
        self.module_dir = os.path.abspath("expoSeq")
        self.save_settings_path = os.path.join(self.module_dir,
                                          "settings",
                                          "save_settings.txt")
        with open(self.save_settings_path, "r") as f:
            save_settings = f.read()
        self.save_settings = literal_eval(save_settings)
    def change_dpi(self, dpi):
        self.save_settings["dpi"] = dpi
        with open(self.save_settings_path, "w") as f:
            f.write(str(self.save_settings))
    def change_save_path(self, path = None, filechooser = True):
        if filechooser == True:
            self.save_settings["fname"] = askopenfilename()
        else:
            self.save_settings["fname"] = path
        with open(self.save_settings_path, "w") as f:
            f.write(str(self.save_settings))
    def change_format(self, format):
        self.save_settings["format"] = format
        with open(self.save_settings_path, "w") as f:
            f.write(str(self.save_settings))