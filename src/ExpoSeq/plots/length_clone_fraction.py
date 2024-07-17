import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns
from textwrap import wrap
from .global_font import font_settings_title, font_settings_normal
class PrepareData:
    @staticmethod
    def color_setup(sequencing_report:pd.DataFrame, no_sequences = 10):
        sequencing_report = sequencing_report.copy()
        sequencing_report['color'] = 'lightgray'

        # Use a Seaborn color palette
        colors = sns.color_palette("husl", no_sequences)

        # Convert the colors to hex format
        colors = [matplotlib.colors.rgb2hex(color) for color in colors]
        # Set the color for the top 10 rows
        colored_report = pd.DataFrame()

        # Group the DataFrame by sample and assign colors within each group
        for sample, group in sequencing_report.groupby('Experiment'):
            group = group.copy()
            group['color'] = 'lightgray'
            
            # Set the color for the top sequences in the group
            for i in range(min(no_sequences, len(group))):
                group.iloc[i, group.columns.get_loc('color')] = colors[i]
            
            colored_report = pd.concat([colored_report, group])

        # Reset the index of the DataFrame
        colored_report.reset_index(drop=True, inplace=True)
        return colored_report, colors
        
    def tidy(self, sequencing_report:pd.DataFrame,region, no_sequences ):

        sequencing_report["length_aaSeq"] = sequencing_report[region].str.len()
        
        sequencing_report,colors = self.color_setup(sequencing_report, no_sequences)
        
        return sequencing_report, colors
    
    @staticmethod
    def filter_by_sample(sequencing_report, sample:str, ):
        return sequencing_report[sequencing_report["Experiment"] == sample]
    
    

class LengthSeqFraction:
    def __init__(self, sequencing_report:pd.DataFrame, sample:str, region:str,
                 ax = None, no_sequences = 10, legend_params = {}, font_settings = {}):
        """This class creates a plot which visualizes the clone fraction distribution of one sample for all available sequence lengths. 
        

        Args:
            sequencing_report (pd.DataFrame): _description_
            sample (str): _description_
            region (str): region on your heavy chain
            ax (_type_, optional): ax object from matplotlib. Defaults to None.
            no_sequences (int, optional): Number of sequences you want to visualized in the stacked bar plot. Defaults to 10.
            legend_params (dict, optional): dictionary containing kwargs for matplotlib legend. Defaults to {}.
            font_settings (dict, optional): dictinoary containing kwargs for matplotlib font settings. Defaults to {}.
        """
        PrepData = PrepareData()
        tidied_report,colors = PrepData.tidy(sequencing_report, region, no_sequences)
        self.tidied_report = PrepData.filter_by_sample(tidied_report, sample)
        grouped = self.tidied_report.groupby(['length_aaSeq', 'color'])['cloneFraction'].sum().reset_index()
        self.grouped = grouped.sort_values(by='cloneFraction', ascending=False)
        self.ax = ax
        if self.ax:
            self.plot(colors, region, no_sequences)
            if legend_params != {}:
                self.make_legend(legend_params)
            if font_settings != {}:
                self.arange_labels(font_settings, sample)
            
        
    def make_legend(self, legend_params:dict):
        self.ax.legend(**legend_params)

    def arange_labels(self, font_settings:dict, sample):
        self.ax.set_ylabel("Clone Fraction", **font_settings_normal)  # Y label
        self.ax.set_xlabel("Sequence Length", **font_settings_normal)  # X label

        original_fontsize = font_settings_normal["fontsize"]

        title = "\n".join(wrap("Length Distribution of " + sample, 40))
        plt.title(title, pad=12, **font_settings_title)
        font_settings["fontsize"] = original_fontsize

        
    def plot(self, colors:list, region:str, no_sequences:int):

        # Plot the bar chart
        unique_colors = self.grouped['color'].unique()
        sorted_colors = [color for color in unique_colors if color in colors] + [color for color in unique_colors if color not in colors]
        cumulative_clone_fractions = {}

        for color in sorted_colors:
            data = self.grouped[self.grouped['color'] == color]
            bar_bottoms = [cumulative_clone_fractions.get(length, 0) for length in data['length_aaSeq']]
            self.ax.bar(data['length_aaSeq'], data['cloneFraction'], bottom=bar_bottoms, color=color)
            for index, row in data.iterrows():
                cumulative_clone_fractions[row['length_aaSeq']] = cumulative_clone_fractions.get(row['length_aaSeq'], 0) + row['cloneFraction']
    
        # Create a legend with the sequences
        for i in range(min(no_sequences, len(self.tidied_report))):
            self.ax.plot([], [], color=colors[i], label=self.tidied_report.iloc[i, self.tidied_report.columns.get_loc(region)])
        



# Create a new column 'color'. Set it to 'gray' for all rows

