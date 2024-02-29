import numpy as np
import matplotlib.pyplot as plt
from ..settings.layout_finder import best_layout
from textwrap import wrap



class LengthDistributionSingle:
    def __init__(self, sequencing_report, sample, region_of_interest, ax = None, font_settings = {}, title_type = "single"):
        unique_length, counts_length = self.tidy(sequencing_report, sample, region_of_interest)
        self.plot(unique_length, counts_length, sample, ax, font_settings, title_type)

        if title_type == "single":
            self.title(sample, font_settings)
        else:
            pass
    @staticmethod
    def tidy(sequencing_report, sample, region_string):
        batch = sequencing_report[sequencing_report["Experiment"] == sample]
        length = batch[region_string].str.len()
        unique_length, counts_length = np.unique(np.array(length)
                                                    , return_counts = True)
        return unique_length, counts_length
    @staticmethod
    def plot(unique_length, counts_length, sample, ax,font_settings, title_type = "single"):
        ax.bar(unique_length, counts_length,  color = "lightskyblue")  # Or whatever you want in the subplot
        #ax.set_xticks(range(0, max_length + 1, 1), range(0, max_length + 1, 1))
        ax.title.set_text(sample)
        if title_type == "single":    
            ax.title.set_size(18)
        else:
            ax.title.set_size(12)
        ax.set_ylabel("Read Count",
                        **font_settings)  # Y label
        ax.set_xlabel('Read Length',
                    **font_settings)  # X label

        original_fontsize = font_settings["fontsize"]
        font_settings["fontsize"] = 20
        title = "\n".join(wrap("Length Distribution of " + sample, 40))
        plt.title(title,
                pad=12,
                **font_settings)
        font_settings["fontsize"] = original_fontsize

        
    @staticmethod
    def title(sample, font_settings):
        original_fontsize = font_settings["fontsize"]
        font_settings["fontsize"] = 20
        title = "\n".join(wrap("Length Distribution of " + sample, 40))
        plt.title(title,
                pad=12,
                **font_settings)
        font_settings["fontsize"] = original_fontsize
        


def length_distribution_multi(fig, sequencing_report, samples, font_settings, region_string, test_version = False,):
    if samples == "all":
        unique_experiments = sequencing_report["Experiment"].unique()
        unique_experiments = np.sort(unique_experiments)
    else:
        sequencing_report = sequencing_report[sequencing_report['Experiment'].isin(samples)]
        unique_experiments = sequencing_report["Experiment"].unique()
        unique_experiments = np.sort(unique_experiments)
    Tot = unique_experiments.shape[0]
    Rows, Cols = best_layout(Tot)
    Position = range(1, Tot + 1)
    n = 0

   # fig = plt.figure(1, constrained_layout=True)
    for experiment in unique_experiments:
            # add every single subplot to the figure with a for loop
        ax = fig.add_subplot(Rows, Cols, Position[n])
        adapted_fontsize = 10 - int(Cols) + 2
        font_settings["fontsize"] = adapted_fontsize
        LengthDistributionSingle(sequencing_report, sample = experiment, region_of_interest = region_string, ax = ax, font_settings = font_settings, title_type="multi")
        n += 1
    title = "\n".join(wrap("Length Distribution of given Samples", 40))
    fig.suptitle(title)







 #   plt.show()
   # plt.bar(unique_length, counts_length)
   # p_values = counts_length/np.sum(counts_length)
    #normalized on one:
   # plt.bar(unique_length, counts_length/np.max(counts_length))
    #p-values
  #  plt.bar(unique_length, p_values)