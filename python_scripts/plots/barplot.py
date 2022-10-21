import matplotlib.pyplot as plt
from python_scripts.plots.plot_params.open_txtfiles import openParams
from python_scripts.tidy_data.barplot import cleaning_data
def barplot(common_vars, local_pattern_more_digits):
    boxplot_data_frame = cleaning_data(common_vars = common_vars,
                                        local_pattern_more_digits = local_pattern_more_digits)
    plt.bar(boxplot_data_frame.Experiment, boxplot_data_frame.tot_sequenced_reads, label="Total Sequenced Reads",
            color="orange")
    plt.bar(boxplot_data_frame.Experiment, boxplot_data_frame.Aligned_Reads, label="Aligned Reads", color="royalblue")
    plt.xticks(rotation=45, ha = 'right', size=5)
    params_legend = openParams("USQ_plot_legend_params.txt")
    plt.legend(**params_legend)
    plot_style = openParams('plot_style.txt')
    plt.ylabel("Reads Count", **plot_style)
    plt.xlabel("Sample", **plot_style)