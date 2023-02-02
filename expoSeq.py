from python_scripts.plots.usq_plot import plot_USQ
from python_scripts.plots.plt_heatmap import plot_heatmap
from python_scripts.plots.logo_plot import plot_logo_multi, plot_logo_single
from python_scripts.plots.length_distribution import length_distribution_multi
from python_scripts.augment_data.binding_data import collect_binding_data
from python_scripts.plots.levenshtein_clustering import clusterSeq, cluster_single_AG
from python_scripts.plots.embedding_with_binding import cluster_toxins_tsne
from python_scripts.plots.cluster_embedding import show_difference
from python_scripts.plots.saveFig import saveFig
import matplotlib.pyplot as plt
from python_scripts.augment_data.uploader import upload
from ast import literal_eval
from settings.plot_styler import PlotStyle
import pandas as pd
import pickle
from python_scripts import chat_bot
import os
from settings.change_save_settings import Change_save_settings
class PlotManager:
    def __init__(self, test_version = False):
        if os.path.isdir("my_experiments"):
            pass
        else:
            os.mkdir("my_experiments")
        if os.path.isdir("temp"):
            pass
        else:
            os.mkdir("temp")
        if test_version == False:
            self.sequencing_report, self.alignment_report, self.experiment = upload()
            with open("my_experiments/" + self.experiment + "/experiment_names.pickle", "rb") as f:
                self.unique_experiments = pickle.load(f)
            self.sequencing_report["Experiment"] = self.sequencing_report["Experiment"].map(self.unique_experiments)
            self.add_binding = input("Do you have binding Data? Y/n")
            if self.add_binding.lower() in ["Y", "y"]:
                self.binding_data = collect_binding_data()
            else:
                self.binding_data = self.add_binding
        else:
            with open("test_data/sequencing_report.csv", "rb") as f:
                self.sequencing_report = pd.read_csv(f, sep = ",")
            with open("settings/global_vars.txt", "r") as f:
                self.global_params = f.read()
            self.global_params = literal_eval(self.global_params)
            with open("test_data/experiment_names.pickle", "rb") as f:
                self.unique_experiments = pickle.load(f)
            self.sequencing_report["Experiment"] = self.sequencing_report["Experiment"].map(self.unique_experiments)
            self.binding_data = pd.read_table("test_data/binding_data.txt", sep=",")
            self.binding_data.drop(self.binding_data.columns[0], axis=1, inplace=True)
        with open('settings/font_settings.txt', "r") as f:
            font_settings = f.read()
        self.font_settings = literal_eval(font_settings)
        with open('settings/legend_settings.txt', "r") as f:
            legend_settings = f.read()
        self.legend_settings = literal_eval(legend_settings)
        self.zero = 0
        self.batch_size = 300
        self.fig = plt.figure(1)
        self.ax = self.fig.gca()
        self.plot_type = "multi"
        self.style = PlotStyle(self.ax, self.plot_type)
        self.settings_saver = Change_save_settings()

    def askMe(self):
        """
        :return: calls the chatbot which can help you to customize your plots or with other question in life and science.
        """
        chat_bot.askMe(self.global_params)
    def print_antigens(self):
        """
        :return: prints the antigens (columns) of your binding data
        """
        print(self.binding_data.columns.to_list()[1:-1])
    def print_samples(self):
        """
        :return: prints the names of your samples, so you can insert them in lists or similar for the analysis with some plots
        """
        print(self.unique_experiments.keys())
    def save(self):
        """
        :return: can be used to save your plots. It will ask you automatically for the directory where you want to save it
        """
        saveFig()
    def close(self):
        plt.close()

    def change_experiment_names(self, specific = None):
        """
        :param specific: optional parameter. You can use this function to change the names of the samples.
        :return:
        """
        if specific == None:
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
        self.sequencing_report["Experiment"] = self.sequencing_report["Experiment"].map(self.unique_experiments)

    def usqPlot(self, library):
        """
        :param library: you insert a string which is a substring of a sample family
        :return: USQ stands for unique sequences quality and the plot shows you the depth of unique sequences which can be used for evaluating your sequencing quality.
        """
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
        """
        :param sample: insert the sample name
        :param highlight_specific_pos: optional. you can highlight a specific position
        :param highlight_pos_range: optional. you can highlight a position range
        :param chosen_seq_length: 16 per default. You always analyze online one sequence length! You can change it if you would like to.
        :return: A logo Plot which shows you the composition of aminoacids per position
        """
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
        """
        :param samples: You analyze all samples per default. If you want to analyze specific samples it has to be a list with the corresponding sample names
        :param chosen_seq_length: 16 per default. You always analyze online one sequence length! You can change it if you would like
        :return: Gives you in one figure one logoPlot per sample.
        """
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
        """
        :param samples: You analyze all samples per default. If you want to analyze specific samples it has to be a list with the corresponding sample names
        :return: Outputs one figurewith one subplot per sample which shows you the distribution of sequence length
        """
        self.fig.clear()
        length_distribution_multi(self.fig,
                            self.sequencing_report,
                            samples,
                            font_settings=self.font_settings)
        self.plot_type = "multi"
        self.ax = self.fig.gca()
        self.style = PlotStyle(self.ax, self.plot_type)


    def basic_cluster(self, sample):
        """

        :param sample: type in a sample name you want to analyze
        :return:
        """
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
        self.plot_type = "multi"
        self.ax = self.fig.gca()
        self.style = PlotStyle(self.ax, self.plot_type)
   # def cluster_multiple_AG(self, sample, antigens):
    #    self.fig.clear()
     #   self.ax = self.fig.gca()
      #  cluster_antigens(self.sequencing_report,
       #                  sample,
        #                 antigens,
         #                self.batch_size)
        #self.plot_type = "multi"
        #self.ax = self.fig.gca()
        #self.style = PlotStyle(self.ax, self.plot_type)
        #self.zero = 1
    def tsne_cluster_AG(self, sample, toxins, pca_components = 70, perplexity = 25, iterations_tsne = 2500):
        """
        :param sample: the sample you would like to analyze
        :param toxins: the toxins you would like to cluster
        :param pca_components: optional. Default is 70
        :param perplexity: optional. Default 25
        :param iterations_tsne: optional. Default is 2500
        :return: It first embeds the sequences in a vector space and then clusters them with PCA and TSNE. The sequences with the binding data are processed with the input sequences, to enable the plotting of the binding data.
        """
        self.fig.clear()
        self.ax = self.fig.gca()
        cluster_toxins_tsne(self.fig,
                            self.ax,
                            self.sequencing_report,
                            sample,
                            toxins,
                            self.binding_data,
                            self.font_settings,
                            pca_components,
                            perplexity,
                            iterations_tsne)
        self.plot_type = "single"
        self.ax = self.fig.gca()
        self.style = PlotStyle(self.ax, self.plot_type)
    def embedding_tsne(self,
                       samples,
                       strands = True,
                       pca_components = 80,
                       perplexity = 30,
                       iterations_tsne = 2500):
        """
        :param samples: the samples you would like to compare towards their sequences
        :param strands: Default is True. It means that you will plot a batch of the strands in your plot
        :param pca_components: Default is 80. Has to be applied for better accuracy of t-SNE. You can indirectly change the described variance with this.
        :param perplexity: Default is 30. It roughly determines the number of nearest neighbors that are considered in the embedding. A higher perplexity value results in a more global structure in the low-dimensional embedding, while a lower perplexity value emphasizes local structure. The optimal perplexity value for a given dataset depends on the dataset's intrinsic dimensionality, and it is usually determined by trial and err
        :param iterations_tsne: Default is 2500. number of times that the algorithm will repeat the optimization process for reducing the cost function. The optimization process aims to minimize the difference between the high-dimensional and low-dimensional representations of the data. More iterations result in a more optimized low-dimensional representation, but also increases the computational cost.
        :return:
        """
        self.fig.clear()
        self.ax = self.fig.gca()
        show_difference(self.sequencing_report,
                        samples,
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
    def morosita_horn(self, specific_samples = False):
        """
        :param specific_experiments: you can give a list with specific experiments
        :return: Returns a matrix of the identity between your samples based on the Morosita Horn Index
        """
        self.fig.clear()
        self.ax = self.fig.gca()
        plot_heatmap(self.sequencing_report,
                     True,
                     "morosita_horn",
                     self.ax,
                     specific_experiments = specific_samples,
                     )
        self.plot_type = "single"
        self.ax = self.fig.gca()
        self.style = PlotStyle(self.ax, self.plot_type)

    def jaccard(self, specific_samples = False):
        """
        :param specific_experiments: give a list with specific samples you would like to analyze
        :return: Returns a matrix of the identity between your samples based on the Jaccard Index
        """
        self.fig.clear()
        self.ax = self.fig.gca()
        plot_heatmap(self.sequencing_report,
                     True,
                     "jaccard",
                     self.ax,
                     specific_experiments = specific_samples,
                     )
        self.plot_type = "single"
        self.ax = self.fig.gca()
        self.style = PlotStyle(self.ax, self.plot_type)
    def sorensen(self, specific_samples = False):
        """
        :param specific_samples: give a list with specific samples you would like to analyze
        :return: Returns a matrix of the identity between your samples based on the Sorensen Dice Index
        """
        self.fig.clear()
        self.ax = self.fig.gca()
        plot_heatmap(self.sequencing_report,
                     True,
                     "sorensen",
                     self.ax,
                     specific_experiments = specific_samples,
                     )
        self.plot_type = "single"
        self.ax = self.fig.gca()
        self.style = PlotStyle(self.ax, self.plot_type)
    def relative(self, specific_samples = False):
        """
        :param specific_samples: give a list with specific samples you would like to analyze
        :return: Returns a matrix of the identity between your samples based on the proportion of identical sequences
        """
        self.fig.clear()
        self.ax = self.fig.gca()
        plot_heatmap(self.sequencing_report,
                     True,
                     "relative",
                     self.ax,
                     specific_experiments = specific_samples,
                     )
        self.plot_type = "single"
        self.ax = self.fig.gca()
        self.style = PlotStyle(self.ax, self.plot_type)









