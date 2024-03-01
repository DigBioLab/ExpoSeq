import seaborn as sns
import matplotlib.pyplot as plt
from textwrap import wrap


class LengthPlotMultiple:
    def __init__(self, sequencing_report,region_of_interest, ax = None, font_settings = {}, plot_type = "boxplot" ):
        self.ax = ax
        
        sequencing_report = self.tidy(sequencing_report, region_of_interest)
        if ax != None:
            if plot_type == "boxplot":
                self.make_boxplot(sequencing_report)
            else:
                self.make_violin_plot(sequencing_report)
        if font_settings != {}:
            self.add_labels(font_settings)
            self.title(font_settings)
            
    @staticmethod
    def tidy(sequencing_report, region_of_interest):
        region = region_of_interest.replace("aaSeq", "nSeq")
        length_region = sequencing_report[region].str.len()
        sequencing_report["length_region"] = length_region
        return sequencing_report
        
    def make_boxplot(self, sequencing_report):
        sns.boxplot(x = "length_region", y = "Experiment", data = sequencing_report, palette = "Dark2", ax = self.ax)
        sns.despine(top=True, right=True, bottom=True, left=True, ax = self.ax)
        
    def make_violin_plot(self, sequencing_report):
        sns.violinplot(sequencing_report, x='length_region', y='Experiment', inner='box', palette='Dark2', ax = self.ax)
        sns.despine(top=True, right=True, bottom=True, left=True, ax = self.ax)
    
    def add_labels(self, font_settings):
        self.ax.set_ylabel("Read Count",
                        **font_settings)  # Y label
        self.ax.set_xlabel('Read Length',
                    **font_settings)  # X label
        
        
    @staticmethod
    def title(font_settings):
        original_fontsize = font_settings["fontsize"]
        font_settings["fontsize"] = 20
        title = "\n".join(wrap("Length Distribution of all samples ", 40))
        plt.title(title,
                pad=12,
                **font_settings)
        font_settings["fontsize"] = original_fontsize
        
        


