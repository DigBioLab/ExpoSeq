import pandas as pd
from python_scripts.plots.usq_plot import plot_USQ
from python_scripts.plots.plt_heatmap import plot_heatmap
class Toolbox(experiment = "chris_pannings", rename_from_dic = True):
    def __init__(self, experiment):
        try:
            with open("my_experiments/" + experiment + "/panning_report", "rb") as f:
                self.sequencing_report = pd.read_table(f,
                                                       sep=",")
        except:
            raise("The Report could not be uploaded. Check the correct experiment name which you can find in the my experiments directory. ")
        self.rename_from_dic = rename_from_dic


    def usq_plot(self, library):
        plot_USQ(sequencing_report = self.sequencing_report,
                 lib_name = library)

    class Identity(analyze_protein = True, choose_experiments = "all"):
        def __init__(self):
            self.specific_experiments = choose_experiments
            self.protein = analyze_protein
        def morosita_horn(self):
            plot_heatmap(self.sequencing_report,
                         self.protein,
                         self.specific_experiments,
                         rename_from_dic = self.rename_from_dic)
        def




