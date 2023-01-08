import pandas as pd
from python_scripts.plots.usq_plot import plot_USQ
from python_scripts.plots.plt_heatmap import plot_heatmap
from python_scripts.plots.logo_plot import plot_logo
from python_scripts.plots.length_distribution import length_distribution
from python_scripts.augment_data.binding_data import collect_binding_data
from python_scripts.plots.levenshtein_clustering import clusterSeq, cluster_single_AG, cluster_antigens
from python_scripts.plots.saveFig import saveFig
from matplotlib.pyplot import close
from mixcr_nils import process_mixcr
import os
from python_scripts.augment_data.load_data import collectFiles
from python_scripts.augment_data.loop_collect_reports import collect_nocluster_files, load_alignment_reports
from python_scripts.tidy_data.trimming import trimming
from tkinter import filedialog
from ast import literal_eval

def upload():
    with open('global_vars.txt', "r") as f:
        data = f.read()
    data = literal_eval(data)
    try:
        last_experiment = data["last_experiment"]
        if os.path.isfile("my_experiments/"  + last_experiment + "/" + "sequencing_report.txt"):
            continue_analysis = input("Do you want to continue to analyze with " + last_experiment + "? Y/n")

            if continue_analysis.lower() in ["n", "N"]:
                next_step = input("If you want to upload a new experiment press 1. If you want to choose another experiment press 2")
                if next_step == "1":
                    experiment = input("How do you want to call your new experiment?")
                    choose_method = input("If you want to process your data with mixcr press 1. If you want to upload already processed txt files with unique clones, press 2")
                    if choose_method == "1":
                        process_mixcr(experiment)
                        with open("my_experiments/" + experiment + "/sequencing_report.txt", "rb") as f:
                            sequencing_report = pd.read_table(f, sep=",")
                    if choose_method == "2":
                        print("Choose the folder which contains the txt files")
                        filenames = collectFiles()
                        try:
                            sequencing_report, each_instance = collect_nocluster_files(filenames)
                            sequencing_report = trimming(sequencing_report,
                                                        divisible_by=3,
                                                        min_count=1,
                                                        new_fraction="cloneFraction")

                        except ValueError:
                            print("It seems that in your directory are either txt files which are not dataframes or you have not uploaded any")
                        try:
                            all_alignment_reports = load_alignment_reports(filenames)
                        except ValueError:
                            print("No files could be found")
                if next_step == "2":
                    while True:
                        user_input = input("Enter the name of the experiment you want to analyze")

                        user_input = user_input  # Try to convert the input to an integer
                        if os.path.isdir("my_experiments/" + user_input):  # Check if the input is in the correct range
                            break  # If the input is valid, break out of the loop
                        else:
                            print("The experiment name does not exist in my_experiments. Please enter the correct name")
                else:
                    pass

            else:
                with open("my_experiments/" + last_experiment + "/sequencing_report.txt", "rb") as f:
                    sequencing_report = pd.read_table(f, sep=",")
        else:
            pass
    except:
        experiment = input("How do you want to call your new experiment?")
        process_mixcr(experiment)
        with open("my_experiments/" + experiment + "/sequencing_report.txt", "rb") as f:
            sequencing_report = pd.read_table(f, sep=",")

    return sequencing_report




class PlotManager:
    def __init__(self):
        self.sequencing_report = upload()
        self.zero = 0
        self.batch_size = 300
        self.add_binding = input("Do you have binding Data? Y/n")
        if self.add_binding.lower() in ["Y", "y"]:
            self.binding_data = collect_binding_data()
        else:
            self.binding_data = self.add_binding

    def print_samples(self):
        print(self.sequencing_report["Experiment"].unique())

    def save_and_close(self):
        if self.zero != 0:
            _save = input("Do you want to save the former plot? Y/n")
            if _save.lower() in ["Y", "y"]:
                saveFig()
            else:
                pass
            close()
        else:
            pass


    def usqPlot(self, library):
        self.save_and_close()
        plot_USQ(sequencing_report = self.sequencing_report,
                 library = library)
        self.zero = 1
    def logoPlot(self, samples = "all", highlight_specific_pos = False, highlight_pos_range = False, chosen_seq_length = 16):
        self.save_and_close()
        plot_logo(self.sequencing_report,
                  samples,
                  highlight_specific_pos,
                  highlight_pos_range)
        self.zero = 1
    def lengthDistribution(self, samples = "all"):
        self.save_and_close()
        length_distribution(self.sequencing_report,
                            samples)
        self.zero = 1

    def basic_cluster(self, sample):
        self.save_and_close()
        clusterSeq(self.sequencing_report, sample, self.batch_size)
        self.zero = 1

    def cluster_one_AG(self, antigen, specific_experiments=False):
        self.save_and_close()
        cluster_single_AG(self.sequencing_report,
                          antigen,
                          self.binding_data,
                          specific_experiments,
                          self.batch_size)
        self.zero = 1

    def cluster_multiple_AG(self, sample, antigens):
        self.save_and_close()
        cluster_antigens(self.sequencing_report, sample, antigens, self.batch_size)
        self.zero = 1

    def morosita_horn(self):
        self.save_and_close()
        plot_heatmap(self.sequencing_report,
                     True,
                     "morosita_horn",
                     specific_experiments = False,
                     )
        self.zero = 1

    def jaccard(self):
        self.save_and_close()
        plot_heatmap(self.sequencing_report,
                     True,
                     "jaccard",
                     specific_experiments = False,
                     )
        self.zero = 1

    def sorensen(self):
        self.save_and_close()
        plot_heatmap(self.sequencing_report,
                     True,
                     "sorensen",
                     specific_experiments = False,
                     )
        self.zero = 1

    def relative(self):
        self.save_and_close()
        plot_heatmap(self.sequencing_report,
                     True,
                     "relative",
                     specific_experiments = False,
                     )
        self.zero = 1









