from python_scripts.plots.usq_plot import plot_USQ
from python_scripts.plots.plt_heatmap import plot_heatmap
from python_scripts.plots.logo_plot import plot_logo_multi, plot_logo_single
from python_scripts.plots.length_distribution import length_distribution_multi, length_distribution_single
from python_scripts.augment_data.binding_data import collect_binding_data
from python_scripts.plots.levenshtein_clustering import clusterSeq, cluster_single_AG, cluster_antigens
from python_scripts.plots.cluster_embedding import show_difference
from python_scripts.plots.saveFig import saveFig
import matplotlib.pyplot as plt
from uploader import upload
from tkinter import filedialog
from ast import literal_eval
from python_scripts.plots.barplot import barplot
from plot_styler import PlotStyle
import pandas as pd
import pickle

class PlotManager:
    def __init__(self, test_version = False):
        if test_version == False:
            self.sequencing_report, self.alignment_report, self.experiment = upload()
            with open("my_experiments/" + self.experiment + "/experiment_names.pickle", "rb") as f:
                self.unique_experiments = pickle.load(f)

        else:
            with open("test_data/sequencing_report.txt", "rb") as f:
                self.sequencing_report = pd.read_table(f, sep = ",")
            with open("global_vars.txt", "r") as f:
                self.global_params = f.read()
            self.global_params = literal_eval(self.global_params)
            with open("test_data/experiment_names.pickle", "rb") as f:
                self.unique_experiments = pickle.load(f)
        with open('font_settings.txt', "r") as f:
            font_settings = f.read()
        self.font_settings = literal_eval(font_settings)
        with open('legend_settings.txt', "r") as f:
            legend_settings = f.read()
        self.legend_settings = literal_eval(legend_settings)
        self.zero = 0
        self.batch_size = 300
        self.add_binding = input("Do you have binding Data? Y/n")
        self.fig = plt.figure(1)
        self.ax = self.fig.gca()
        self.style = PlotStyle(self.ax)
        if self.add_binding.lower() in ["Y", "y"]:
            self.binding_data = collect_binding_data()
        else:
            self.binding_data = self.add_binding
    def print_samples(self):
        print(self.unique_experiments)
    def save(self):
        saveFig()
    def close(self):
        plt.close()

    def change_experiment_names(self, specific = False):
        if specific == False:
            for key in self.unique_experiments:
                print(f"Current value for {key}: {self.unique_experiments[key]}")
                new_value = input("Enter a new value or press any key to skip")
                if len(new_value) > 1:
                    self.unique_experiments[key] = new_value
                else:
                    pass
        else:
            new_value = input("Enter the new name for " + specific)
            self.unique_experiments[specific] = new_value


    def usqPlot(self, library):
        self.fig.clear()
        self.ax = self.fig.gca()
        plot_USQ(fig = self.fig,
                sequencing_report = self.sequencing_report,
                 library = library,
                 font_settings = self.font_settings,
                 legend_settings = self.legend_settings)
        self.plot_type = "single"
        self.ax = self.fig.gca()
        self.style = PlotStyle(self.ax, self.plot_type)

    def logoPlot_single(self,
                        sample,
                        highlight_specific_pos = False,
                        highlight_pos_range = False,
                        chosen_seq_length = 16):
        self.fig.clear()
        plot_logo_single(self.ax,
                         self.sequencing_report,
                         sample,
                         self.font_settings,
                         highlight_specific_pos,
                         highlight_pos_range,
                         chosen_seq_length)
        self.plot_type = "single"
        self.ax = self.fig.gca()
        self.style = PlotStyle(self.ax, self.plot_type)

    def logoPlot_multi(self,
                 samples = "all",
                 chosen_seq_length = 16):
        self.fig.clear()
        plot_logo_multi(self.fig,
                  self.sequencing_report,
                  samples,
                  self.font_settings,
                  chosen_seq_length
                  )

        self.plot_type = "multi"
        self.ax= self.fig.gca()
        self.style = PlotStyle(self.ax, self.plot_type)
    def lengthDistribution_multi(self, samples = "all"):
        self.fig.clear()
        length_distribution(self.fig,
                            self.sequencing_report,
                            samples,
                            font_settings=self.font_settings)


    def basic_cluster(self, sample):
        self.fig.clear()
        self.ax = self.fig.gca()
        clusterSeq(
                   self.ax,
                   self.sequencing_report,
                   sample,
                   self.batch_size)
        self.zero = 1
        self.ax = self.fig.gca()
        self.plot_type = "single"
        self.style = PlotStyle(self.ax, self.plot_type)
    def cluster_one_AG(self, antigen, specific_experiments=False):
        self.fig.clear()
        self.ax = self.fig.gca()
        cluster_single_AG(self.sequencing_report,
                          antigen,
                          self.binding_data,
                          specific_experiments,
                          self.batch_size)
        self.zero = 1
    def cluster_multiple_AG(self, sample, antigens):
        self.fig.clear()
        self.ax = self.fig.gca()
        cluster_antigens(self.sequencing_report,
                         sample,
                         antigens,
                         self.batch_size)
        self.plot_type = "multi"
        self.ax = self.fig.gca()
        self.style = PlotStyle(self.ax, self.plot_type)
        self.zero = 1
    def embedding_tsne(self,
                       list_experiments,
                       strands = True,
                       pca_components = 80,
                       perplexity = 30,
                       iterations_tsne = 2500):
        self.fig.clear()
        self.ax = self.fig.gca()
        show_difference(self.sequencing_report,
                        list_experiments,
                        strands,
                        self.batch_size,
                        pca_components,
                        perplexity,
                        iterations_tsne,
                        self.ax)
        self.plot_type = "single"
        self.ax = self.fig.gca()
        self.style = PlotStyle(self.ax, self.plot_type)
        self.zero = 1
    def morosita_horn(self):
        self.fig.clear()
        self.ax = self.fig.gca()
        plot_heatmap(self.sequencing_report,
                     True,
                     "morosita_horn",
                     self.ax,
                     specific_experiments = False,
                     )
        self.plot_type = "single"
        self.ax = self.fig.gca()
        self.style = PlotStyle(self.ax, self.plot_type)


    def jaccard(self):
        self.fig.clear()
        self.ax = self.fig.gca()
        plot_heatmap(self.sequencing_report,
                     True,
                     "jaccard",
                     self.ax,
                     specific_experiments = False,
                     )
        self.plot_type = "single"
        self.ax = self.fig.gca()
        self.style = PlotStyle(self.ax, self.plot_type)
    def sorensen(self):
        self.fig.clear()
        self.ax = self.fig.gca()
        plot_heatmap(self.sequencing_report,
                     True,
                     "sorensen",
                     self.ax,
                     specific_experiments = False,
                     )
        self.plot_type = "single"
        self.ax = self.fig.gca()
        self.style = PlotStyle(self.ax, self.plot_type)
    def relative(self):
        self.fig.clear()
        self.ax = self.fig.gca()
        plot_heatmap(self.sequencing_report,
                     True,
                     "relative",
                     self.ax,
                     specific_experiments = False,
                     )
        self.plot_type = "single"
        self.ax = self.fig.gca()
        self.style = PlotStyle(self.ax, self.plot_type)










