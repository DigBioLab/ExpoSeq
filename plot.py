import pandas as pd
from python_scripts.plots.usq_plot import plot_USQ
from python_scripts.plots.plt_heatmap import plot_heatmap
from python_scripts.plots.logo_plot import plot_logo
from python_scripts.plots.length_distribution import length_distribution
class PlotManager():
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
        length_distribution(self.sequencing_report, samples)




    class Identity():
        def __init__(self, analyze_protein = True, choose_experiments = "all"):
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









