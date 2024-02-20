import numpy as np 
import matplotlib.pyplot as plt
from textwrap import wrap
import pandas as pd

class PrepareData:
    @staticmethod
    def check_input(sequencing_report, region_of_interest, method):
        assert isinstance(sequencing_report, pd.DataFrame)
        assert type(region_of_interest) == str
        assert region_of_interest in sequencing_report.columns.to_list()
        assert type(method) == str
        assert method in ["Shannon", "InverseSimpson"]
        
    @staticmethod
    def calc_simpson_index(clones):
        return 1/ (np.sum(clones**2))
    
    @staticmethod
    def calc_shannon_index(clones):
        return -np.sum(clones * np.log(clones))
    
    def cleaning(self, sequencing_report, region_of_interest, method):
        self.check_input(sequencing_report, region_of_interest, method)
        values = []
        unique_experiments = sequencing_report["Experiment"].unique().tolist()
        for experiment in unique_experiments:
            exp_spec = sequencing_report[sequencing_report["Experiment"] == experiment]
            aa_seqs = exp_spec[region_of_interest]
            clones = exp_spec["cloneFraction"]
            if method == "InverseSimpson":
                values.append(self.calc_simpson_index(clones))
            if method == "Shannon":
                values.append(self.calc_shannon_index(clones))
        return values, unique_experiments
    
    
    

class DiversityPlot:
    def __init__(self, sequencing_report,  region_of_interest,ax = None, font_settings = {},  method = "InverseSimpson",):
        self.ax = ax
        self.method = method
        self.font_settings = font_settings
        values, unique_experiments = PrepareData().cleaning(sequencing_report, region_of_interest, method )
        if ax != None:
            self.ax = ax
            self.create_base_plot(values, unique_experiments)
            self.add_plot_addons(unique_experiments)
            self.set_title()

    
    def create_base_plot(self, values, unique_experiments, alpha = 1, color = "lightskyblue"):     
        self.ax.bar(x = unique_experiments,
                    height = values,
                    label = "Diversity",
                    color = color,
                    alpha = alpha)
        
    def add_plot_addons(self, unique_experiments):
        if self.method == "InverseSimpson":
            self.ax.set_ylabel("Inverse Simpson Index", **self.font_settings)
            plt.yscale("log")
        if self.method == "Shannon":
            self.ax.set_ylabel("Shannon Index", **self.font_settings)
        self.ax.set_xlabel("Sample", **self.font_settings)
        self.ax.set_xticklabels(
                labels = unique_experiments,
                rotation=45,
                ha = 'right',
                size=12)
        
    def set_title(self):
        if len(self.font_settings) != 0:
            original_fontsize = self.font_settings["fontsize"]
            self.font_settings["fontsize"] = 22
            if self.method == "InverseSimpson":
                title = "Diversity based on Inverse Simpson Index"
            if self.method == "Shannon":
                title = "Diversity based on Shannon Index"
            title = "\n".join(wrap(title, 40))
            self.ax.set_title(title,pad = 12, **self.font_settings)
            self.font_settings["fontsize"] = original_fontsize


