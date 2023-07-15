import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
from ExpoSeq.tidy_data.heatmaps import tidy_jaccard, tidy_sorensen, tidy_morosita_horn, tidy_heatmap_share

import numpy as np

def plot_heatmap(sequencing_report, protein, heatmap, ax, colorbar_settings, font_settings, specific_experiments = False):

    if heatmap == "morosita_horn":
        unique_sequences, unique_experiments = tidy_morosita_horn.cleaning_data(sequencing_report,
                                                                                protein = protein,
                                                                                specific_experiments = specific_experiments)
        matrix, unique_sequences, unique_experiments = tidy_morosita_horn.morosita_horn_matrix(unique_sequences,
                                                                                               unique_experiments)
        title = "Identity between samples based on Morosita Horn index"
    if heatmap == "jaccard":
        matrix, unique_sequences, unique_experiments = tidy_jaccard.cleaning_jaccard(sequencing_report,
                                                                                     protein=protein,
                                                                                     specific_experiments = specific_experiments)
        title = "Identity between samples based on Jaccard index"
    if heatmap == "sorensen":
        matrix, unique_sequences, unique_experiments = tidy_sorensen.heatmap_sorensen(sequencing_report,
                                                                                      protein = protein,
                                                                                      specific_experiments = specific_experiments)
        title = "Identity between samples based on Sorensen index"
    if heatmap == "relative":
        matrix, unique_sequences, unique_experiments = tidy_heatmap_share.heatmap_share(sequencing_report,
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
                cbar_kws = colorbar_sets)

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






