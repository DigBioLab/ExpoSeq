import logomaker
from python_scripts.plots.plot_params.open_txtfiles import openParams
import numpy as np
import matplotlib.pyplot as plt
from python_scripts.tidy_data.tidy_seqlogoPlot import cleaning

def plot_logo(sequencing_report, samples, highlight_specific_pos, highlight_pos_range, chosen_seq_length = 16):

    if samples == "all":
        unique_experiments = sequencing_report["Experiment"].unique()
        unique_experiments = np.sort(unique_experiments)
    else:
        unique_experiments = sequencing_report["Experiment"].unique()
        unique_experiments = np.array([i for i in unique_experiments if i in samples])
        unique_experiments = np.sort(unique_experiments)
    # Subplots are organized in a Rows x Cols Grid
    # Tot and Cols are known
    Tot = unique_experiments.shape[0]
    Cols = int(input("How many Columns for the Plots do you want?"))
    # Compute Rows required
    Rows = Tot // Cols
    # EDIT for correct number of rows:
    # If one additional row is necessary -> add one:
    if Tot % Cols != 0:
        Rows += 1
    # Create a Position index
    Position = range(1, Tot + 1)
    n = 0
    fig = plt.figure(1, constrained_layout=True)
    for i in unique_experiments:
        aa_distribution, sequence_length, length_filtered_seqs = cleaning(i, sequencing_report, chosen_seq_length)
        if length_filtered_seqs != 0:
            if length_filtered_seqs < 100:
                print("only " + str(length_filtered_seqs) + " sequences with the given length were found. The results might be biased")
            ax = fig.add_subplot(Rows, Cols, Position[n])
            logo_plot = logomaker.Logo(aa_distribution,
                                        shade_below=.5,
                                        fade_below=.5,
                                        font_name='Arial Rounded MT Bold',
                                        ax=ax)
            plot_style = openParams('plot_style.txt')
            plt.title(i, fontsize = 8) # check out https://matplotlib.org/3.1.0/api/_as_gen/matplotlib.pyplot.title.html
                # logo_plot.ax.set_ylabel("Amino Acid Frequency", fontsize=5)
                # logo_plot.ax.set_xlabel("Sequence Position", fontsize=5)
                #logo_plot.ax.set_xticks(range(0, sequence_length), range(1, sequence_length + 1))
            if highlight_specific_pos != False:
                logo_plot.highlight_position(p = 5, color = 'gold', alpha = .5)
            if highlight_pos_range != False:
                logo_plot.ax.highlight_position_range(pmin = 3, pmax = 5, color = "lightcyan")
            n = n + 1
        else:
            print("Sample " + i + "was skipped because no sequence was found")

    fig.suptitle("Logo Plots for sequence Length " + str(chosen_seq_length))


