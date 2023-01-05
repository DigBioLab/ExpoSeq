import seaborn as sns
from os.path import dirname, abspath
from ast import literal_eval
import matplotlib
import matplotlib.pyplot as plt
from python_scripts.tidy_data import heatmaps

def plot_heatmap(sequencing_report, protein, heatmap, specific_experiments = False, rename_from_dic = True):

    if heatmap == "morosita_horn":
        unique_sequences, unique_experiments = heatmaps.tidy_morosita_horn.cleaning_data(sequencing_report,
                                                             protein = True,
                                                             specific_experiments = specific_experiments)
        matrix, unique_sequences, unique_experiments = heatmaps.tidy_heatmap.morosita_horn_matrix(unique_sequences,
                                                                                                  unique_experiments)
    if heatmap == "jaccard":
        matrix, unique_sequences, unique_experiments = heatmaps.tidy_jaccard.cleaning_jaccard(sequencing_report,
                                                                                              protein=protein)

    if heatmap == "sorensen":
        matrix, unique_sequences, unique_experiments = heatmaps.tidy_sorensen(sequencing_report,
                                                                              protein = protein)
    if heatmap == "relative":
        matrix, unique_sequences, unique_experiments = heatmaps.tidy_heatmap_share(sequencing_report,
                                                                                    protein=protein)

    matrix = matrix.sort_index(axis=1)
    matrix = matrix.sort_index(axis = 0)
    matplotlib.use('Qt5AGG')
    ROOT = dirname(abspath('ExpoSeq'))
    with open(ROOT + '\python_scripts\plots\plot_params\heatmap_plot_params.txt') as f:
        data = f.read()
    data = literal_eval(data)
   # plt.figure(
    #           dpi = 300)
    sns.heatmap(matrix,
                **data)
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





