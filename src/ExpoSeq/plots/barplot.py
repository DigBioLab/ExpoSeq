import matplotlib.pyplot as plt
import numpy as np
from textwrap import wrap
from pandas import DataFrame


class PrepareData:
    def cleaning_data(all_alignment_reports, sequencing_report):
        all_alignment_reports = all_alignment_reports.sort_values("Input file(s)").reset_index(drop=True)
        unique_experiments = sequencing_report.sort_values("Experiment")["Experiment"].unique()
        unique_experiments = list(unique_experiments)
        tot_sequencing_reads = all_alignment_reports["Total sequencing reads"]
        aligned_reads = all_alignment_reports["Successfully aligned reads"].str.split("(", expand = True).iloc[:, 0]
        data = DataFrame([unique_experiments, aligned_reads, tot_sequencing_reads]).T
        data.columns = ["Experiment", "Aligned_Reads", "tot_sequenced_reads"]
        return data

class AlignmentPlot:
    def __init__(self, sequencing_report, all_alignment_reports, ax = None, font_settings = {}, legend_settings = {}, apply_log = False) -> None:
        self.ax = ax
        self.data = PrepareData().cleaning_data(all_alignment_reports, sequencing_report)
        self.base_plot()
        self.log_x_axis(apply_log)
        self.add_labels(font_settings)
        self.add_legend()
        
    def base_plot(self):
        self.ax.bar(self.data.Experiment,
           np.array(self.data.tot_sequenced_reads).astype(np.float32),
           label="Total Sequenced Reads",
            color="lightsalmon",
            alpha = 1)
        self.ax.bar(self.data.Experiment,
            np.array(self.data.Aligned_Reads).astype(np.float32),
            label="Aligned Reads",
            color="lightskyblue",
            alpha = 1)
        plt.xticks(rotation=45,
              ha = 'right',
              size=12)
        
    @staticmethod
    def log_x_axis(apply_log):
        if apply_log:
            plt.yscale("log")
        
    def add_labels(self, font_settings):     
        if font_settings != {}:   
            plt.ylabel("Reads Count", **font_settings)
            plt.xlabel("Sample", **font_settings)
            original_fontsize = font_settings["fontsize"]
            font_settings["fontsize"] = 22
            title = "\n".join(wrap("Alignment Quality of the analyzed samples", 40))
            self.ax.set_title(title,pad = 12, **font_settings)
            font_settings["fontsize"] = original_fontsize
    
    def add_legend(self, legend_settings):
        if legend_settings != {}:
            self.ax.legend(**legend_settings)

        
        
def barplot(ax, all_alignment_reports, sequencing_report_all,font_settings, legend_settings, apply_log):
    boxplot_data_frame = cleaning_data(all_alignment_reports,
                                       sequencing_report_all)
    ax.bar(boxplot_data_frame.Experiment,
           np.array(boxplot_data_frame.tot_sequenced_reads).astype(np.float32),
           label="Total Sequenced Reads",
            color="lightsalmon",
            alpha = 1)
    ax.bar(boxplot_data_frame.Experiment,
            np.array(boxplot_data_frame.Aligned_Reads).astype(np.float32),
            label="Aligned Reads",
            color="lightskyblue",
            alpha = 1)
    plt.xticks(rotation=45,
              ha = 'right',
              size=12)

    ax.legend(**legend_settings)
    plt.ylabel("Reads Count", **font_settings)
    plt.xlabel("Sample", **font_settings)

    if apply_log == True:
        plt.yscale("log")
    plt.tight_layout()
    original_fontsize = font_settings["fontsize"]
    font_settings["fontsize"] = 22
    title = "\n".join(wrap("Alignment Quality of the analyzed samples", 40))
    ax.set_title(title,pad = 12, **font_settings)
    font_settings["fontsize"] = original_fontsize


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
