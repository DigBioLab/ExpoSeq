import matplotlib.pyplot as plt
from python_scripts.plots.plot_params.open_txtfiles import openParams
from python_scripts.tidy_data.barplot import cleaning_data
import numpy as np
import matplotlib
def barplot(all_alignment_reports, sequencing_report_all, apply_log = True):
    matplotlib.use("Qt5Agg")
    boxplot_data_frame = cleaning_data(all_alignment_reports,
                                       sequencing_report_all)
    plt.bar(boxplot_data_frame.Experiment, np.array(boxplot_data_frame.tot_sequenced_reads).astype(np.float32), label="Total Sequenced Reads",
            color="orange")
    plt.bar(boxplot_data_frame.Experiment,
            np.array(boxplot_data_frame.Aligned_Reads).astype(np.float32),
            label="Aligned Reads",
            color="royalblue")
    plt.xticks(rotation=45, ha = 'right', size=5)
    params_legend = openParams("USQ_plot_legend_params.txt")
    plt.legend(**params_legend)
    plot_style = openParams('plot_style.txt')
    plt.ylabel("Reads Count", **plot_style)
    plt.xlabel("Sample", **plot_style)
    if apply_log == True:
        plt.yscale("log")
    plt.tight_layout()


   # matplotlib.use("Qt5Agg")
   # boxplot_data_frame = cleaning_data(all_alignment_reports, sequencing_report_all)
   # ax.bar(boxplot_data_frame.Experiment, np.array(boxplot_data_frame.tot_sequenced_reads).astype(np.float32),label="Total Sequenced Reads",color="orange", np.array(boxplot_data_frame.Aligned_Reads).astype(np.float32),
   # label="Aligned Reads",
   # color="royalblue")
    #ax.set_xticks(rotation=45, ha = 'right', size=5)
    #params_legend = openParams("USQ_plot_legend_params.txt")
    #ax.legend(**params_legend)
    #plot_style = openParams('plot_style.txt')
    #ax.set_ylabel("Reads Count", **plot_style)
    #ax.set_xlabel("Sample", **plot_style)
    #if apply_log == True:
     #   ax.set_yscale("log")
    #fig.tight_layout()