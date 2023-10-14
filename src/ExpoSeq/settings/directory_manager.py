
import os
from ast import literal_eval


class Directories:
    def __init__(self):
        self.module_dir = os.path.abspath("")
        self.common_vars = os.path.join(self.module_dir, "settings", "global_vars.txt")
        self.font_settings_path = os.path.join(self.module_dir, "settings", "font_settings.txt")
        self.legend_settings_path = os.path.join(self.module_dir,
                                                 "settings",
                                                 "legend_settings.txt")
        self.colorbar_path = os.path.join(self.module_dir,
                                          "settings",
                                          "colorbar.txt")

    def check_dirs(self):
        if not os.path.isdir(os.path.join(self.module_dir, "my_experiments")):
            print("create my_experiments directory")
            os.mkdir(os.path.join(self.module_dir, "my_experiments"))
        if not os.path.isdir(os.path.join(self.module_dir, "temp")):
            print("create temp directory")
            os.mkdir(os.path.join(self.module_dir, "temp"))
        if not os.path.isdir(os.path.join(self.module_dir, "settings")):
            print("create settings directory")
            os.mkdir(os.path.join(self.module_dir, "settings"))

    def create_global_params(self):
        common_vars = {'mixcr_path': '', 'last_experiment': '', 'api_gpt3': '', 'region_of_interest': '', 'RAM': '',
                       'clustalw_path': ''}
        with open(self.common_vars, "w") as f:
            f.write(str(common_vars))

    def create_font_settings(self):
        font_settings = {'fontfamily': 'serif', 'fontsize': '18', 'fontstyle': 'normal', 'fontweight': 'bold'}
        with open(os.path.join(self.module_dir, "settings", "font_settings.txt"), "w") as f:
            f.write(str(font_settings))

    def create_legend_settings(self):
        legend_settings = {'loc': 'upper right', 'bbox_to_anchor': (1, 1), 'ncols': 1, 'fontsize': 16, 'frameon': True,
                           'framealpha': 1, 'facecolor': 'white', 'mode': None, 'title_fontsize': 'small'}
        with open(os.path.join(self.module_dir, "settings", "legend_settings.txt"), "w") as f:
            f.write(str(legend_settings))

    def create_colorbar_settings(self):
        colorbar = {'cmap': 'inferno', 'orientation': 'vertical', 'spacing': 'proportional', 'extend': 'neither'}
        with open(self.colorbar_path, "w") as f:
            f.write(str(colorbar))

    def read_global_params(self):
        if not os.path.isfile(self.common_vars):
            self.create_global_params()
        with open(self.common_vars, "r") as f:
            global_params = f.read()
        global_params = literal_eval(global_params)
        return global_params

    def read_font_settings(self):
        if not os.path.isfile(self.font_settings_path):
            self.create_font_settings()
        with open(self.font_settings_path, "r") as f:
            font_settings = f.read()
        font_settings = literal_eval(font_settings)
        return font_settings

    def read_legend_settings(self):
        if not os.path.isfile(self.legend_settings_path):
            self.create_legend_settings()
        with open(self.legend_settings_path, "r") as f:
            legend_settings = f.read()
        legend_settings = literal_eval(legend_settings)
        return legend_settings

    def read_colorbar_settings(self):
        if not os.path.isfile(self.colorbar_path):
            self.create_colorbar_settings()
        with open(self.colorbar_path, "r") as f:
            colorbar_settings = f.read()
        colorbar_settings = literal_eval(colorbar_settings)
        return colorbar_settings

    def get_experiment_path(self, experiment):
        experiment_path = os.path.join(self.module_dir,
                                       "my_experiments",
                                       experiment,
                                       "experiment_names.pickle")
        return experiment_path

