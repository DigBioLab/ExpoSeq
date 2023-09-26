import matplotlib.pyplot as plt
import logomaker
from ExpoSeq.tidy_data.tidy_seqlogoPlot import cleaning
import numpy as np




class LogoPlot:
    def __init__(self,ax, sequencing_report, region_string, sample, chosen_seq_length, highlight_spec_position, font_settings):
        self.ax = ax
        self.aa_distribution = self.cleaningPlot(sample, sequencing_report, chosen_seq_length, region_string)
       # self.createPlot()
        self.chosen_seq_length = chosen_seq_length
        self.highlight_spec_position = highlight_spec_position
        self.font_settings = font_settings
        self.func_type = {
    "alpha": [0.01, 0.01, 0.01, 1],
    "shade": [0.01, 0.01, 0.01, 1],
    "fade": [0.01, 0.01, 0.01, 1],
            }
        self.button_type = {
        "alpha": "CustomScale",
            "shade": "CustomScale",
            "fade": "CustomScale",
        }
        
        
                
    def cleaningPlot(self, sample, sequencing_report, chosen_seq_length, region_string):
        aa_distribution, sequence_length, length_filtered_seqs = cleaning(sample,
                                                                    sequencing_report,
                                                                    chosen_seq_length,
                                                                    region_string)
        return aa_distribution, sequence_length, length_filtered_seqs
    
    def createPlot(self,shade_below = .5, fade_below = .5, font_name = 'Arial Rounded MT Bold', color_scheme = "skylign_protein", show_spines = False,):

        self.logo_plot = logomaker.Logo(
                            self.aa_distribution,
                            shade_below=shade_below,
                            fade_below=fade_below,
                            font_name=font_name,
                            color_scheme=color_scheme,
                            show_spines=show_spines,
                            ax=self.ax,
                            )
        self.logo_plot.style_xticks(anchor=1,
                    spacing=1,
                    rotation=0)
        

    def add_style(self, highlight_specific_pos):
        original_fontsize = self.font_settings["fontsize"]
        self.ax.set_xlabel("Frequency", **self.font_settings)
        self.ax.set_ylabel("Position on sequence", **self.font_settings)
        self.font_settings["fontsize"] = 22
        plt.title("Logo Plot of " + self.sample + " with sequence length " + str(self.chosen_seq_length), **self.font_settings)
        self.font_settings["fontsize"] = original_fontsize
        labels_true = list(range(0, self.chosen_seq_length))
        numbers_true = list(range(1, self.chosen_seq_length + 1))
        plt.xticks(labels_true, numbers_true)
        if highlight_specific_pos != False:
            self.logo_plot.highlight_position(p=5,
                                        color='gold',
                                        alpha=.5)




def plot_logo_single(ax, sequencing_report, sample, font_settings, highlight_specific_pos, region_string, chosen_seq_length = 16):
    aa_distribution, sequence_length, length_filtered_seqs = cleaning(sample,
                                                                      sequencing_report,
                                                                      chosen_seq_length,
                                                                      region_string)
    logo_plot = logomaker.Logo(aa_distribution,
                               shade_below=.5,
                               fade_below=.5,
                               font_name='Arial Rounded MT Bold',
                               color_scheme="skylign_protein",
                               show_spines=False,
                               ax=ax,
                               )
    logo_plot.style_xticks(anchor=1,
                           spacing=1,
                           rotation=0)
    original_fontsize = font_settings["fontsize"]
    ax.set_xlabel("Frequency", **font_settings)
    ax.set_ylabel("Position on sequence", **font_settings)
    font_settings["fontsize"] = 22
   # plt.title("Logo Plot of " + sample + " with sequence length " + str(chosen_seq_length), **font_settings)
    font_settings["fontsize"] = original_fontsize
    labels_true = list(range(0, chosen_seq_length))
    numbers_true = list(range(1, chosen_seq_length + 1))
    plt.xticks(labels_true, numbers_true)
    if highlight_specific_pos != False:
        logo_plot.highlight_position(p=5,
                                     color='gold',
                                     alpha=.5)
 #   if highlight_pos_range != False:
  #      logo_plot.ax.highlight_position_range(pmin=3,
   #                                           pmax=5,
    #                                          color="lightcyan")


def plot_logo_multi(fig, sequencing_report, samples, font_settings,region_string, chosen_seq_length = 16,):
    if samples == "all":
        unique_experiments = sequencing_report["Experiment"].unique()
        unique_experiments = np.sort(unique_experiments)
    else:
        unique_experiments = sequencing_report["Experiment"].unique()
        unique_experiments = np.array([i for i in unique_experiments if i in samples])
        unique_experiments = np.sort(unique_experiments)
    Tot = unique_experiments.shape[0]
    Rows, Cols = best_layout(Tot)
    Position = range(1, Tot + 1)
    n = 0
    adapted_fontsize = 10 - int(Cols) + 2
    original_fontsize = font_settings["fontsize"]
    font_settings["fontsize"] = adapted_fontsize
  #  fig = plt.figure(1, constrained_layout=True)
    for i in unique_experiments:
        aa_distribution, sequence_length, length_filtered_seqs = cleaning(i,
                                                                          sequencing_report,
                                                                          chosen_seq_length,
                                                                          region_string)
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
                                        color_scheme="skylign_protein",
                                        show_spines=False,
                                        ax=ax,
                                        )
            #logo_plot.set_xticks(range(aa_distribution.shape[0]))
            logo_plot.style_xticks(anchor=0,
                                   spacing=1,
                                   rotation=0,)
            plt.title(i, **font_settings) # check out https://matplotlib.org/3.1.0/api/_as_gen/matplotlib.pyplot.title.html


            n = n + 1
        else:
            print("Sample " + i + "was skipped because no sequence was found")
    original_fontsize = font_settings["fontsize"]
    font_settings["fontsize"] = 22
    fig.suptitle("Logo Plots for sequence Length " + str(chosen_seq_length), **font_settings)
    font_settings["fontsize"] = original_fontsize


