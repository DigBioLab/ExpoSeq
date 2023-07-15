from ExpoSeq.tidy_data.tidy_USQ_plot_ import cleaning_data
import matplotlib.pyplot as plt
import random
from collections import Counter
import warnings
import matplotlib
matplotlib.use('Qt5Agg')
def plot_USQ(fig, sequencing_report, samples, font_settings, legend_settings):
    sub_table = cleaning_data(sequencing_report = sequencing_report,
                              samples = samples)
    x_axis = []
    y_axis = []
    for sample in range(len(sub_table)):
        sequences = list(sub_table.iloc[sample]["nSeqCDR3"])
        fraction = list(sub_table.iloc[sample]["clonesFraction"])
        counter = 0
        temp_list = []
        unique_list = []
        total_list = []
        mean_count = (sum(sub_table.iloc[sample]["readCount"]) / 100)
        while counter < 100:
            temp = random.choices(sequences, weights=fraction, k = int(mean_count))  # randomly sample 100th of the sequences with the recalculated clonefractions as probablities
            temp_list.extend(temp)  ## extend the temp_list with the samples from each 'pull-out' of the 100 'pull-outs'
            unique = Counter(temp_list).keys()  ## count unique seequences from the increasing temp_list
            unique_list.append(unique)  ## for each run append the number of unique sequences to a new cell in this list - will be accumulative
            total_list.append(len(temp_list))  ## append the total sequences amount for plotting
            counter += 1
        plt.plot(total_list,
                 [len(x) for x in unique_list],
                 figure = fig,
                 label = sub_table.iloc[sample].Experiment,
                 alpha = 1,
                fillstyle = 'full',
                linewidth = 0.8,
        )
        x_axis.append(total_list)
        y_axis.append([len(x) for x in unique_list])
    if len(sub_table) > 7:
        warnings.warn("Risk of Overplotting: Many different colors are used. It may be hard to differ between the data")
  #  for i in range(len(x_axis)):
   #     plt.plot(x_axis[i], y_axis[i], **params_plot)
    plt.legend(sub_table.Experiment, title = "Sample Names", **legend_settings)
    plt.xlabel("Total sampled sequences", **font_settings)
    plt.ylabel("Total Unique Sequences", **font_settings)
    original_fontsize = font_settings["fontsize"]
    font_settings["fontsize"] = 22
    plt.title("Sequencing depth of the given samples ",pad = 12, **font_settings)
    font_settings["fontsize"] = original_fontsize
    fig.tight_layout()

