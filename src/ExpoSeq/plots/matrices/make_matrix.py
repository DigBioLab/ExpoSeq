import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from textwrap import wrap
import ExpoSeq.plots.matrices.morosita_horn_matrix as morosita_horn_matrix
import ExpoSeq.plots.matrices.sorensen_matrix as sorensen_matrix
import ExpoSeq.plots.matrices.jaccard_matrix as jaccard_matrix
import ExpoSeq.plots.matrices.relative_matrix as relative_matrix

class IdentityMatrix:
    def __init__(self, sequencing_report, region_string,matrix_type, colorbar_settings, specific_experiments = False,ax = None, protein = True, font_settings = {}, 
                 cmap = 'inferno', fmt = ".2f", annotate_cells = False):
        self.ax = ax

        self.matrix, self.unique_experiments = self.tidy_data(matrix_type, specific_experiments, sequencing_report, region_string, protein)

        self.font_settings = font_settings
        if ax != None:
            self.createPlot(cmap = cmap, fmt = fmt, annotate_cells = annotate_cells, colorbar_settings=colorbar_settings)
            if font_settings != {}:
                self.add_style()
                
    @staticmethod
    def tidy_data(matrix_type, specific_experiments, sequencing_report, region_string, protein):
        if matrix_type == "morosita_horn":
            matrix, unique_experiments = morosita_horn_matrix.PrepareData().cleaningPlot(specific_experiments, sequencing_report, protein, region_string)
        if matrix_type == "sorensen":
            matrix, unique_experiments = sorensen_matrix.PrepareData().sorensen_matrix(sequencing_report, protein, region_string, specific_experiments)        
        if matrix_type == "jaccard":
            matrix, unique_experiments = jaccard_matrix.PrepareData().cleaning_jaccard(sequencing_report, region_string, protein, specific_experiments)
        if matrix_type == "relative":
            matrix, unique_experiments = relative_matrix.PrepareData().heatmap_share(sequencing_report, region_string, protein, specific_experiments)
        
        return matrix, unique_experiments

    def createPlot(self,linewidth = 0,line_color = "white", annotate_cells = False, annot_kws = 6, fmt = ".2f", cmap = 'inferno', colorbar_settings = None):
        annot_kws = {"size": annot_kws}
        matrix = self.matrix.sort_index(axis=1)
        matrix = matrix.sort_index(axis = 0)
        if colorbar_settings != None:
            colorbar_sets = {**colorbar_settings,
                            **{'label': "Degree of Identity"}}
        else:
            colorbar_sets = {}
            
        self.morosita_horn_matrix = sns.heatmap(matrix,
                                                cbar_kws=colorbar_sets,
                                                ax = self.ax,
                                                linewidths=linewidth,
                                                linecolor = line_color,
                                                annot = annotate_cells,
                                                annot_kws=annot_kws,
                                                fmt = fmt,
                                                cmap= cmap)

    def add_style(self):
        plt.xticks(ticks = np.arange(0.5, len(self.unique_experiments) + 0.5 ,1),
                labels = self.unique_experiments,
                rotation = 45,
                ha = 'right',
                size = 12) # create a function which finds the perfect size based on counts of xlabels
        plt.yticks(ticks = np.arange(0.5, len(self.unique_experiments) + 0.5 ,1),
                labels = self.unique_experiments,
                rotation = 0,
                va='center',
                size=12)
        original_fontsize = self.font_settings["fontsize"]
        self.font_settings["fontsize"] = 22
        title = "Identity between samples based on Morosita Horn index"
        title = "\n".join(wrap(title, 40))
        self.ax.set_title(title,pad = 12, **self.font_settings)
        self.font_settings["fontsize"] = original_fontsize
