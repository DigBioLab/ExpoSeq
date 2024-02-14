from textwrap import wrap
import pandas as pd

class StackedAADistribution:
    def __init__(self, sequencing_report, sample, region, region_of_interest, protein = True, font_settings = {}, ax = None):
        aa_distribution = self.prepare_data(sequencing_report, sample,region, region_of_interest)
        
        self.plot(aa_distribution, ax, sample, region, protein, font_settings)
    
    @staticmethod
    def prepare_data(sequencing_report, sample,region, region_string):
        sample = sequencing_report[sequencing_report["Experiment"] == sample]
        local_report = sample[["Experiment", "cloneFraction", region_string]]
        aminoacids = "ACDEFGHIKLMNPQRSTVWY"

        sequences = local_report[local_report[region_string].astype(str).str.len() >= region[1]][region_string]
        max_length = local_report[region_string].str.len().max()
        if not region[1] <= max_length:
            raise ValueError("you upper region limit is above the longest sequence. That is not possible. Please reduce it.")
        compDict = {aa: max_length*[0] for aa in aminoacids}
        for seq in sequences:
            for aa_position in range(len(seq)):
                aminoacid = seq[aa_position]
                if aminoacid == '*':
                    pass
                else:
                    compDict[aminoacid][aa_position] += 1
        aa_distribution = pd.DataFrame.from_dict(compDict)
        aa_distribution = aa_distribution.divide(aa_distribution.sum(axis=1), axis=0)
        aa_distribution = aa_distribution[(aa_distribution.index >= region[0]) & (aa_distribution.index <= region[1])]
        return aa_distribution
    
    @staticmethod
    def get_colors(aa_distribution):
        color_scheme = {
        'F': [.16, .99, .18],
        'Y': [.04, .40, .05],
        'L': [.99, .60, .25],
        'V': [1.0, .80, .27],
        'I': [.80, .60, .24],
        'H': [.40, .02, .20],
        'W': [.42, .79, .42],
        'A': [.99, .60, .42],
        'S': [.04, .14, .98],
        'T': [.17, 1.0, 1.0],
        'M': [.80, .60, .80],
        'N': [.21, .40, .40],
        'Q': [.40, .41, .79],
        'R': [.59, .02, .04],
        'K': [.40, .20, .03],
        'E': [.79, .04, .22],
        'G': [.95, .94, .22],
        'D': [.99, .05, .11],
        'P': [.10, .61, .99],
        'C': [.09, .60, .60],
        }
        color_list = [color_scheme[aa] for aa in aa_distribution.columns]
        return color_list

    def plot(self, aa_distribution, ax, sample, region, protein, font_settings):
        color_list = self.get_colors(aa_distribution)
        aa_distribution.plot(kind='bar', stacked=True,color=color_list, ax = ax)
        ax.set_xlabel('Position on amino acid sequence', **font_settings)
        if protein == True:
            ax.set_ylabel('Relatvie Proportion of Amino Acid', **font_settings)
        else:
            ax.set_ylabel('Relatvie Proportion of Nucleotide', **font_settings)

    #    ax.set_xticks(rotation = 360)
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        original_fontsize = font_settings["fontsize"]
        font_settings["fontsize"] = 22
        if protein == True:
            title = "\n".join(wrap("Stacked Amino Acid Distribution of " + sample, 40))
            ax.set_title(title,pad = 12, **font_settings)
        else:
            title = "\n".join(wrap("Stacked Nucleotide Distribution of " + sample, 40))
            ax.set_title(title,pad = 12, **font_settings)
        ax.set_xticklabels([*range(region[0], region[1]+1)], rotation=0, ha='center')
        font_settings["fontsize"] = original_fontsize
        
        

        
