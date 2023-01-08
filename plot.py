import pandas as pd
from python_scripts.plots.usq_plot import plot_USQ
from python_scripts.plots.plt_heatmap import plot_heatmap
from python_scripts.plots.logo_plot import plot_logo
from python_scripts.plots.length_distribution import length_distribution
from python_scripts.augment_data.binding_data import collect_binding_data
from python_scripts.plots.levenshtein_clustering import clusterSeq, cluster_single_AG, cluster_antigens
class PlotManager:
    def __init__(self, experiment):
        try:
            with open("my_experiments/" + experiment + "/panning_report", "rb") as f:
                self.sequencing_report = pd.read_table(f,
                                                       sep=",")
        except:
            raise("The Report could not be uploaded. Check the correct experiment name which you can find in the my experiments directory. ")
        self.rename_from_dic = rename_from_dic

    def usqPlot(self, library):
        plot_USQ(sequencing_report = self.sequencing_report,
                 lib_name = library)
    def logoPlot(self, samples = "all", highlight_specific_pos = False, highlight_pos_range = False):
        plot_logo(self.sequencing_report,
                  samples,
                  highlight_specific_pos,
                  highlight_pos_range)

    def lengthDistribution(self, samples = "all"):
        length_distribution(self.sequencing_report,
                            samples)

class Clustering(PlotManager):
        def __init__(self, add_binding, batch_size = 300):
            super().__init__()
            if add_binding == True:
                self.binding_data = collect_binding_data()
            else:
                self.binding_data = add_binding
            self.batch_size = batch_size

        def basic_cluster(self, sample):
            clusterSeq(self.sequencing_report, sample, self.batch_size)

        def cluster_one_AG(self, antigen, specific_experiments = False):
            cluster_single_AG(self.sequencing_report,
                              antigen,
                              self.binding_data,
                              specific_experiments,
                              self.batch_size)
        def cluster_multiple_AG(self, sample, antigens):
            cluster_antigens(self.sequencing_report, sample, antigens, self.batch_size)



class Identity(PlotManager):
        def __init__(self, analyze_protein = True, choose_experiments = "all"):
            super().__init__()
            self.specific_experiments = choose_experiments
            self.protein = analyze_protein

        def morosita_horn(self):
            plot_heatmap(self.sequencing_report,
                         self.protein,
                         "morosita_horn",
                         self.specific_experiments,
                         rename_from_dic = self.rename_from_dic)
        def jaccard(self):
            plot_heatmap(self.sequencing_report,
                         self.protein,
                         "jaccard",
                         self.specific_experiments,
                         rename_from_dic = self.rename_from_dic)

        def sorensen(self):
            plot_heatmap(self.sequencing_report,
                         self.protein,
                         "sorensen",
                         self.specific_experiments,
                         rename_from_dic=self.rename_from_dic)
        def relative(self):
            plot_heatmap(self.sequencing_report,
                         self.protein,
                         "relative",
                         self.specific_experiments,
                         rename_from_dic=self.rename_from_dic)









