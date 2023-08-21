import subprocess
import os
import shutil
import os
from ast import literal_eval
import pkg_resources
    


class Settings:
    def __init__(self):
        self.module_dir = os.path.abspath("")
        self.pkg_path = pkg_resources.resource_filename("ExpoSeq", "")
        self.common_vars =  os.path.join(self.pkg_path,"settings", "global_vars.txt")
        with open(self.common_vars, "r") as f:
            global_params = f.read()
        self.global_params = literal_eval(global_params)
        
    def change_ram(self,ram):
        self.global_params["RAM"] = ram
        with open(self.common_vars, "w") as f:
            f.write(str(self.global_params))
    
    def change_global_vars(self, param, input):
        self.global_params[param] = input
        with open(self.common_vars, "w") as f:
            f.write(str(self.global_params))
            
    def change_mixcr_path(self, path):
        self.global_params["mixcr_path"] = path
        with open(self.common_vars, "w") as f:
            f.write(str(self.global_params))
            
    def check_dirs(self):
        if not os.path.isdir(os.path.join(self.module_dir, "my_experiments")):
            print("create my_experiments directory")
            os.mkdir(os.path.join(self.module_dir, "my_experiments"))
        if not os.path.isdir(os.path.join(self.module_dir, "temp")):
            print("create temp directory")
            os.mkdir(os.path.join(self.module_dir, "temp"))
        
    def save_settings(self, path, params):
            with open(path, "w") as f:
                f.write(str(params))
    def remove_directories(directory_path):
    # Iterate over all items in the directory
        for item in os.listdir(directory_path):
            shutil.rmtree(os.path.join(directory_path, item))
            
    def read_font_settings(self):
        font_settings_path = os.path.join(self.pkg_path, "settings", "font_settings.txt")
        assert os.path.isfile(font_settings_path), "The font settings file does not exist in the given filepath"
        with open(font_settings_path, "r") as f:
            font_settings = f.read()
        font_settings = literal_eval(font_settings)
        return font_settings
    
    def read_legend_settings(self):
        legend_settings_path = os.path.join(self.pkg_path,
                                        "settings",
                                        "legend_settings.txt")
        assert os.path.isfile(legend_settings_path), "The legend settings file does not exist in the given filepath"
        with open(legend_settings_path, "r") as f:
            legend_settings = f.read()
        legend_settings = literal_eval(legend_settings)
        return legend_settings
    
    def read_colorbar_settings(self):
        colorbar_path = os.path.join(self.pkg_path,
                                    "settings",
                                    "colorbar.txt")
        assert os.path.isfile(colorbar_path), "The colorbar file does not exist in the given filepath"
        with open(colorbar_path, "r") as f:
            colorbar_settings = f.read()
        colorbar_settings = literal_eval(colorbar_settings)
        return colorbar_settings
    
    def get_experiment_path(self, experiment):
        experiment_path = os.path.join(self.module_dir,
                                "my_experiments",
                                experiment,
                                "experiment_names.pickle")
        return experiment_path


    def reset_pipeline(self):
        global_params = {'mixcr_path': '', 'last_experiment': '', 'api_gpt3': '', 'region_of_interest': '', 'RAM': 42}
        self.save_settings(self.global_params, global_params)
        self.save_settings(os.path.join(self.pkg_path,"settings","save_settings.txt"), {'fname': '', 'format': 'png', 'dpi': 300})
        legend_params = {"loc": "upper right", #'upper left', 'upper right', 'lower left', 'lower right'
                    "bbox_to_anchor": (1, 1),
                    "ncols": 1,
                    "fontsize": "medium", #'xx-small', 'x-small', 'small', 'medium', 'large', 'x-large', 'xx-large'
                    "frameon": True,
                    "framealpha": 1,
                    "facecolor": "white", #background color
                    "mode": None,
                    "title_fontsize": 'small',
                    }
        self.save_settings(os.path.join(self.pkg_path,"settings","legend_settings.txt"), legend_params)
        self.save_settings(os.path.join(self.pkg_path,"settings","font_settings.txt"), {"fontfamily": "serif","fontsize": "14","fontstyle": "normal","fontweight": "bold"})
        self.save_settings(os.path.join(self.pkg_path,"settings","colorbar.txt"), {"cmap": "inferno","orientation": "vertical","spacing": "proportional","extend": "neither"})
        

    def reinstall_pipeline(self,):
        subprocess.run(["pip","-f", "uninstall", "ExpoSeq"])
        subprocess.run(["pip","-f", "install", "ExpoSeq"])
        delete = input("You will now delete all your analyzed experiments. Press enter to cancel this operation or type 'delete' to continue.")
        if delete == "delete":
            self.remove_directories(os.path.join(self.module_dir, "my_experiments"))
            self.remove_directories(os.path.join(self.module_dir, "temp"))
            
            