import matplotlib.pyplot as plt
import random
from collections import Counter
from textwrap import wrap
import pandas as pd

class PrepareData:
    @staticmethod
    def cleaning_data(sequencing_report, samples, region_of_interest, fraction_column):
        sub_table = sequencing_report.loc[sequencing_report["Experiment"].isin(samples)]
        apply_fun = lambda x: list(x)
        sub_table = sub_table.groupby("Experiment")[[region_of_interest, fraction_column, 'readCount']].agg(apply_fun).reset_index()
        return sub_table    
    
    def tidy(self, sequencing_report, samples, region_of_interest, fraction_column = "cloneFraction"):
        region_of_interest = region_of_interest.replace("aaSeq", "nSeq")
        sub_table = self.cleaning_data(sequencing_report, samples, region_of_interest, fraction_column)
        x_axis = []
        y_axis = []
        sample_names = []
        for sample in range(len(sub_table)):
            sequences = list(sub_table.iloc[int(sample)][region_of_interest])
            fraction = list(sub_table.iloc[int(sample)][fraction_column])
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
            
            x_axis.append(total_list)
            y_axis.append([len(x) for x in unique_list])
            sample_names.append(sample)
        results_plot = pd.DataFrame({"samples": sample_names,
                      "x_axis": x_axis,
                      "y_axis": y_axis})
        return results_plot
        


class RarefractionCurves:
    def __init__(self, sequencing_report, samples, region_of_interest, ax = None, font_settings = {}, legend_settings = {}, fraction_column = "cloneFraction"):
        self.ax = ax
        self.plot_data = self.prepare_data(sequencing_report, samples, region_of_interest, fraction_column)
        self.plot()
        self.customize_axis(font_settings)
        self.legend(legend_settings, samples)
        self.title(font_settings)
        
    @staticmethod
    def prepare_data(sequencing_report, samples, region_of_interest, fraction_column):
        return PrepareData().tidy(sequencing_report, samples, region_of_interest, fraction_column=fraction_column)
    
    def plot(self):
        if self.ax != None:
            for i, row in self.plot_data.iterrows():
                self.ax.plot(row['x_axis'], row['y_axis'], label=row['samples'], alpha = 1, fillstyle = 'full', linewidth = 0.8)
    
    @staticmethod
    def customize_axis(font_settings):
        if font_settings != {}:
            plt.xlabel("Total sampled sequences", **font_settings)
            plt.ylabel("Total Unique Sequences", **font_settings)
            
    def legend(self, legend_settings, samples):
        if legend_settings != {}:
            plt.legend(samples, title = "Sample Names", **legend_settings)
      
    @staticmethod      
    def title(font_settings):
        if font_settings != {}:
            original_fontsize = font_settings["fontsize"]
            font_settings["fontsize"] = 22
            
            title = "\n".join(wrap("Sequencing depth of the given samples", 40))
            plt.title(title,pad = 12, **font_settings)
            font_settings["fontsize"] = original_fontsize

