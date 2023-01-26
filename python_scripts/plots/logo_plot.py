import logomaker
from python_scripts.plots.plot_params.open_txtfiles import openParams
import numpy as np
import matplotlib.pyplot as plt
from python_scripts.tidy_data.tidy_seqlogoPlot import cleaning

def plot_logo(fig, sequencing_report, samples,font_settings, highlight_specific_pos, highlight_pos_range, chosen_seq_length = 16):
    if samples == "all":
        unique_experiments = sequencing_report["Experiment"].unique()
        unique_experiments = np.sort(unique_experiments)
    else:
        unique_experiments = sequencing_report["Experiment"].unique()
        unique_experiments = np.array([i for i in unique_experiments if i in samples])
        unique_experiments = np.sort(unique_experiments)
    Tot = unique_experiments.shape[0]
    Cols = int(input("How many Columns for the Plots do you want?"))
    # Compute Rows required
    Rows = Tot // Cols
    if Tot % Cols != 0:
        Rows += 1
    # Create a Position index
    Position = range(1, Tot + 1)
    n = 0
    adapted_fontsize = 10 - Cols + 2
    font_settings["fontsize"] = adapted_fontsize
  #  fig = plt.figure(1, constrained_layout=True)
    for i in unique_experiments:
        aa_distribution, sequence_length, length_filtered_seqs = cleaning(i,
                                                                          sequencing_report,
                                                                          chosen_seq_length)
        if length_filtered_seqs != 0:
            if length_filtered_seqs < 100:
                print("only " + str(length_filtered_seqs) + " sequences with the given length were found. The results might be biased")
            ax = fig.add_subplot(Rows,
                                 Cols,
                                 Position[n],
                                 xticks = (np.arange(0, chosen_seq_length, step = 1)))
            logo_plot = logomaker.Logo(aa_distribution,
                                        shade_below=.5,
                                        fade_below=.5,
                                        font_name='Arial Rounded MT Bold',
                                        ax=ax,
                                        )
            logo_plot.set_xticks(range(aa_distribution.shape[0]))
            logo_plot.style_xticks(anchor=0,
                                   spacing=1,
                                   rotation=0)
            plt.title(i, **font_settings) # check out https://matplotlib.org/3.1.0/api/_as_gen/matplotlib.pyplot.title.html

            if highlight_specific_pos != False:
                logo_plot.highlight_position(p = 5,
                                             color = 'gold',
                                             alpha = .5)
            if highlight_pos_range != False:
                logo_plot.ax.highlight_position_range(pmin = 3,
                                                      pmax = 5,
                                                      color = "lightcyan")
            n = n + 1
        else:
            print("Sample " + i + "was skipped because no sequence was found")

    fig.suptitle("Logo Plots for sequence Length " + str(chosen_seq_length))

