from python_scripts.tidy_data.tidy_USQ_plot_ import cleaning_data
import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
import random
from collections import Counter
import warnings
from python_scripts.plots.plot_params.open_txtfiles import openParams

def plot_USQ(sequencing_report, local_pattern_more_digits):
    sub_table = cleaning_data(sequencing_report_all = sequencing_report,
                              local_pattern_more_digits = local_pattern_more_digits)
    x_axis = []
    y_axis = []
    params_plot = openParams('USQ_plot.txt')
    plot_style = openParams('plot_style.txt')
    fig = plt.figure()
    for sample in range(len(sub_table)):
        sequences = list(sub_table.iloc[sample].nSeqCDR3)
        fraction = list(sub_table.iloc[sample].clonefrac)
        counter = 0
        temp_list = []
        unique_list = []
        total_list = []
        mean_count = (sum(sub_table.iloc[sample]["cloneCount"]) / 100)
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
                 **params_plot)
        x_axis.append(total_list)
        y_axis.append([len(x) for x in unique_list])
    if len(sub_table) > 7:
        warnings.warn("Risk of Overplotting: Many different colors are used. It may be hard to differ between the data")
  #  for i in range(len(x_axis)):
   #     plt.plot(x_axis[i], y_axis[i], **params_plot)
    params_legend = openParams("USQ_plot_legend_params.txt")

    fig.legend(sub_table.Experiment, title = "Sample Names", **params_legend)
    plt.xlabel("Total sampled sequences", **plot_style)
    plt.ylabel("Found Unique Sequences", **plot_style)
   # plt.show(fig)
