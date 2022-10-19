import seaborn as sns
from os.path import dirname, abspath
from ast import literal_eval
import matplotlib
import matplotlib.pyplot as plt
from python_scripts.tidy_data import tidy_heatmap


def plot_heatmap(unique_sequences, unique_experiments):
    matrix = tidy_heatmap.morosita_horn_matrix(unique_sequences, unique_experiments)
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





