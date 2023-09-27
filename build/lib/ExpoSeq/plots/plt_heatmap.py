import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
from ExpoSeq.tidy_data.heatmaps import tidy_jaccard, tidy_sorensen, tidy_morosita_horn, tidy_heatmap_share

import numpy as np


class MorositaHornMatrix:
    def __init__(self,ax, sequencing_report, region_string, samples,protein, font_settings):
        self.ax = ax
        self.samples = samples
        self.matrix = self.cleaningPlot( samples, sequencing_report, protein, region_string, )
       # self.createPlot()

        self.font_settings = font_settings
        self.func_type = {
        "Transparency": ["alpha", [0.01, 1]],
        "Line Width": ["linewidth",[0.01, 1]],
        "annotation size": ["annot_kws",[1, 10]],
        "Color Map": ["cmap", ["virdis", "plasma", "inferno", "magma", "cividis", "BuPu", "PuRd", "jet", "turbo", "hot", "magma", ]],
        "Line Color": ["line_color"]
            }
        self.button_type = {
        "Transparency": "CustomScale",
            "Line Width": "CustomScale",
            "Annotation size": "CustomScale",
            "cmap": "optionlist",
            "Line Color": "SimpleTextField"
        }
        
        
                
    def cleaningPlot(self, samples, sequencing_report, protein, region_string, ):
        unique_sequences, unique_experiments = tidy_morosita_horn.cleaning_data(sequencing_report,
                                                                                region_string,
                                                                                protein = protein,
                                                                                specific_experiments = self.samples,
                                                                                )
        matrix, unique_sequences, unique_experiments = tidy_morosita_horn.morosita_horn_matrix(unique_sequences,
                                                                                        unique_experiments)
        unique_experiments = list(unique_experiments)
        matrix = matrix.sort_index(axis=1)
        matrix = matrix.sort_index(axis = 0)
        return matrix
    
    def createPlot(self,linewidth = 0,line_color = "white", annotate_cells = False, annot_kws = 1, fmt = ".2f", cmap = 'inferno', orientation = "vertical", spacing = 'proportional'):
        annot_kws = {"size": annot_kws}
        unique_experiments = list(unique_experiments)
        matrix = matrix.sort_index(axis=1)
        matrix = matrix.sort_index(axis = 0)
        matplotlib.use('Qt5Agg')
        #colorbar_sets = {**colorbar_settings,
             #               **{'label': "Degree of Identity"}}
        self.morosita_horn_matrix = sns.heatmap(matrix,
                                                ax = self.ax,
                                                linewidths=linewidth,
                                                linecolor = line_color,
                                                annot = annotate_cells,
                                                annot_kws=annot_kws,
                                                fmt = fmt,
                                                cmap= cmap,
                                                orientation = orientation,
                                                spacing = spacing)


        

    def add_style(self, highlight_specific_pos):
        plt.xticks(ticks = np.arange(0.5, len(self.samples) + 0.5 ,1),
                labels = self.samples,
                rotation = 45,
                ha = 'right',
                size = 12) # create a function which finds the perfect size based on counts of xlabels
        plt.yticks(ticks = np.arange(0.5, len(self.samples) + 0.5 ,1),
                labels = self.samples,
                rotation = 0,
                va='center',
                size=12)
        original_fontsize = self.font_settings["fontsize"]
        self.font_settings["fontsize"] = 22
        title = "Identity between samples based on Morosita Horn index"
        self.ax.set_title(title,pad = 12, **self.font_settings)
        self.font_settings["fontsize"] = original_fontsize






def plot_heatmap(sequencing_report, protein, heatmap, ax, colorbar_settings, font_settings, annotate_cells,region_of_interest, specific_experiments = False):

    if heatmap == "morosita_horn":
        unique_sequences, unique_experiments = tidy_morosita_horn.cleaning_data(sequencing_report,
                                                                                region_of_interest,
                                                                                protein = protein,
                                                                                specific_experiments = specific_experiments,
                                                                                )
        matrix, unique_sequences, unique_experiments = tidy_morosita_horn.morosita_horn_matrix(unique_sequences,
                                                                                               unique_experiments)
        title = "Identity between samples based on Morosita Horn index"
    if heatmap == "jaccard":
        matrix, unique_sequences, unique_experiments = tidy_jaccard.cleaning_jaccard(sequencing_report,
                                                                                     region_of_interest,
                                                                                     protein=protein,
                                                                                     specific_experiments = specific_experiments,
                                                                                     )
        title = "Identity between samples based on Jaccard index"
    if heatmap == "sorensen":
        matrix, unique_sequences, unique_experiments = tidy_sorensen.heatmap_sorensen(sequencing_report,
                                                                                      protein,
                                                                                      region_of_interest,
                                                                                      
                                                                                      specific_experiments = specific_experiments)
        title = "Identity between samples based on Sorensen index"
    if heatmap == "relative":
        matrix, unique_sequences, unique_experiments = tidy_heatmap_share.heatmap_share(sequencing_report,
                                                                                        region_of_interest,
                                                                                        protein=protein,
                                                                                        specific_experiments = specific_experiments)
        title = "Identity between samples based on relative coefficient"

    unique_experiments = list(unique_experiments)
    matrix = matrix.sort_index(axis=1)
    matrix = matrix.sort_index(axis = 0)
    matplotlib.use('Qt5Agg')
    colorbar_sets = {**colorbar_settings,
                         **{'label': "Degree of Identity"}}
    sns.heatmap(matrix,
                ax = ax,
                cbar_kws = colorbar_sets,
                annot = annotate_cells,
                annot_kws={"size": 6},
                fmt = ".2f")

    plt.xticks(ticks = np.arange(0.5, len(unique_experiments) + 0.5 ,1),
               labels = unique_experiments,
               rotation = 45,
               ha = 'right',
               size = 12) # create a function which finds the perfect size based on counts of xlabels
    plt.yticks(ticks = np.arange(0.5, len(unique_experiments) + 0.5 ,1),
               labels = unique_experiments,
               rotation = 0,
               va='center',
               size=12)
    original_fontsize = font_settings["fontsize"]
    font_settings["fontsize"] = 22
    ax.set_title(title,pad = 12, **font_settings)
    font_settings["fontsize"] = original_fontsize
    return matrix   






