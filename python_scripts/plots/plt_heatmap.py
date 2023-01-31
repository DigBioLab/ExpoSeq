import seaborn as sns
from os.path import dirname, abspath
from ast import literal_eval
import matplotlib
import matplotlib.pyplot as plt
from python_scripts.tidy_data import heatmaps
from python_scripts.tidy_data.heatmaps import tidy_jaccard, tidy_sorensen, tidy_morosita_horn, tidy_heatmap_share

def plot_heatmap(sequencing_report, protein, heatmap, ax, specific_experiments = False, rename_from_dic = True):

    if heatmap == "morosita_horn":
        unique_sequences, unique_experiments = tidy_morosita_horn.cleaning_data(sequencing_report,
                                                             protein = True,
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

    matrix = matrix.sort_index(axis=1)
    matrix = matrix.sort_index(axis = 0)
    matplotlib.use('Qt5Agg')
    sns.heatmap(matrix,
                ax = ax)
    plt.xticks(rotation = 45, ha = 'right', size = 5) # create a function which finds the perfect size based on counts of xlabels
    plt.yticks(va='top', size=5)
    savefig = ""
   # while savefig != "Y" or "n":
    #    savefig = input("Do you want to save the figure? Y/n" )
     #   if savefig != "Y" or "n":
      #      raise Exception("Sorry, your input was neither Y or n. Try Again, please.")

   # if savefig == "Y":
    #    plt.savefig("morosita_horn_matrix_dna.png", dpi=300)
    #else:
     #   pass
    # close figure





