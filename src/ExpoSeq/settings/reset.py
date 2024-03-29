import shutil
import os
from ast import literal_eval

def remove_directories(directory_path):
    # Iterate over all items in the directory
    for item in os.listdir(directory_path):
        shutil.rmtree(os.path.join(directory_path, item))

def original_settings():
    module_dir = os.getcwd()
    experiments_dir = os.path.join(module_dir,
                                   "my_experiments")
    remove_directories(experiments_dir)
    common_vars_path = os.path.join(module_dir,
                                    "settings",
                                    "global_vars.txt")
    with open(common_vars_path, "r") as f:
        global_params = f.read()
    global_params = {'mixcr_path': '', 'last_experiment': '', 'api_gpt3': '', 'region_of_interest': '', 'RAM': 42}
    with open(common_vars_path, "w") as f:
        f.write(str(global_params))
    save_path = os.path.join(module_dir,
                             "settings",
                             "save_settings.txt")
    save_params = {'fname': '', 'format': 'png', 'dpi': 300}
    with open(save_path, "w") as f:
        f.write(str(save_params))
    legend_settings_path = os.path.join(module_dir,
                                        "settings",
                                        "legend_settings.txt")
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
    with open(legend_settings_path, "w") as f:
        f.write(str(legend_params))
    font_settings_path = os.path.join(module_dir,
                                      "settings",
                                      "font_settings.txt")
    font_settings_params = {
                            "fontfamily": "serif",
                            "fontsize": "14",
                            "fontstyle": "normal",
                            "fontweight": "bold",
                            #rotation
                            }
    with open(font_settings_path, "w") as f:
        f.write(str(font_settings_params))
    colorbar_path = os.path.join(module_dir,
                                 "settings",
                                 "colorbar.txt")
    colorbar_params = {
                        "cmap": "inferno",
                        "orientation": "vertical",
                        "spacing": "proportional",
                        "extend": "neither"
                        }
    with open(colorbar_path, "w") as f:
        f.write(str(colorbar_params))




