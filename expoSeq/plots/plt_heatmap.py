import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
from ..tidy_data.heatmaps import tidy_jaccard, tidy_sorensen, tidy_morosita_horn, tidy_heatmap_share



def plot_heatmap(sequencing_report, protein, heatmap, ax, colorbar_settings, font_settings, specific_experiments = False):

    if heatmap == "morosita_horn":
        unique_sequences, unique_experiments = tidy_morosita_horn.cleaning_data(sequencing_report,
                                                                                protein = protein,
                                                                                specific_experiments = specific_experiments)
        matrix, unique_sequences, unique_experiments = tidy_morosita_horn.morosita_horn_matrix(unique_sequences,
                                                                                               unique_experiments)
    if heatmap == "jaccard":
        matrix, unique_sequences, unique_experiments = tidy_jaccard.cleaning_jaccard(sequencing_report,
                                                                                     protein=protein,
                                                                                     specific_experiments = specific_experiments)

    if heatmap == "sorensen":
        matrix, unique_sequences, unique_experiments = tidy_sorensen.heatmap_sorensen(sequencing_report,
                                                                                      protein = protein,
                                                                                      specific_experiments = specific_experiments)
    if heatmap == "relative":
        matrix, unique_sequences, unique_experiments = tidy_heatmap_share.heatmap_share(sequencing_report,
                                                                                        protein=protein,
                                                                                        specific_experiments = specific_experiments)

    unique_experiments = list(unique_experiments)
    matrix = matrix.sort_index(axis=1)
    matrix = matrix.sort_index(axis = 0)
    matplotlib.use('Qt5Agg')
    colorbar_settings = {**colorbar_settings,
                         **{'label': "Degree of Identity"}}
    sns.heatmap(matrix,
                ax = ax,
                cbar_kws = colorbar_settings)
    plt.xticks(ticks = range(0, len(unique_experiments), 1),
               labels = unique_experiments,
               rotation = 45,
               ha = 'right',
               size = 5) # create a function which finds the perfect size based on counts of xlabels
    plt.yticks(ticks = range(0, len(unique_experiments), 1),
               labels = unique_experiments,
               va='center',
               size=5)
    original_fontsize = font_settings["fontsize"]
    font_settings["fontsize"] = 22
    ax.set_title("Identity between samples based on Morosita Horn Index", **font_settings)
    font_settings["fontsize"] = original_fontsize






