from textwrap import wrap
import pandas as pd


class StackedAADistribution:
    def __init__(
        self,
        sequencing_report,
        sample,
        region,
        region_of_interest,
        protein=True,
        font_settings={},
        ax=None,
    ):
        aa_distribution = self.prepare_data(
            sequencing_report, sample, region, region_of_interest
        )

        self.plot(aa_distribution, ax, sample, region, protein, font_settings)

    @staticmethod
    def prepare_data(sequencing_report, sample, region, region_string):
        sample = sequencing_report[sequencing_report["Experiment"] == sample]
        local_report = sample[["Experiment", "cloneFraction", region_string]]
        aminoacids = "ACDEFGHIKLMNPQRSTVWY"

        sequences = local_report[
            local_report[region_string].astype(str).str.len() >= region[1]
        ][region_string]
        max_length = local_report[region_string].str.len().max()
        if not region[1] <= max_length:
            raise ValueError(
                "you upper region limit is above the longest sequence. That is not possible. Please reduce it."
            )
        compDict = {aa: max_length * [0] for aa in aminoacids}
        for seq in sequences:
            for aa_position in range(len(seq)):
                aminoacid = seq[aa_position]
                if aminoacid == "*":
                    pass
                else:
                    compDict[aminoacid][aa_position] += 1
        aa_distribution = pd.DataFrame.from_dict(compDict)
        aa_distribution = aa_distribution.divide(aa_distribution.sum(axis=1), axis=0)
        aa_distribution = aa_distribution[
            (aa_distribution.index >= region[0]) & (aa_distribution.index <= region[1])
        ]
        return aa_distribution

    @staticmethod
    def get_colors(aa_distribution):
        color_scheme = {
            "F": [0.16, 0.99, 0.18],
            "Y": [0.04, 0.40, 0.05],
            "L": [0.99, 0.60, 0.25],
            "V": [1.0, 0.80, 0.27],
            "I": [0.80, 0.60, 0.24],
            "H": [0.40, 0.02, 0.20],
            "W": [0.42, 0.79, 0.42],
            "A": [0.99, 0.60, 0.42],
            "S": [0.04, 0.14, 0.98],
            "T": [0.17, 1.0, 1.0],
            "M": [0.80, 0.60, 0.80],
            "N": [0.21, 0.40, 0.40],
            "Q": [0.40, 0.41, 0.79],
            "R": [0.59, 0.02, 0.04],
            "K": [0.40, 0.20, 0.03],
            "E": [0.79, 0.04, 0.22],
            "G": [0.95, 0.94, 0.22],
            "D": [0.99, 0.05, 0.11],
            "P": [0.10, 0.61, 0.99],
            "C": [0.09, 0.60, 0.60],
        }
        color_list = [color_scheme[aa] for aa in aa_distribution.columns]
        return color_list

    def plot(self, aa_distribution, ax, sample, region, protein, font_settings):
        color_list = self.get_colors(aa_distribution)
        aa_distribution.plot(kind="bar", stacked=True, color=color_list, ax=ax)
        ax.set_xlabel("Position on amino acid sequence", **font_settings)
        if protein == True:
            ax.set_ylabel("Relatvie Proportion of Amino Acid", **font_settings)
        else:
            ax.set_ylabel("Relatvie Proportion of Nucleotide", **font_settings)

        #    ax.set_xticks(rotation = 360)
        ax.legend(bbox_to_anchor=(1.05, 1), loc="upper left")
        original_fontsize = font_settings["fontsize"]
        font_settings["fontsize"] = 22
        if protein == True:
            title = "\n".join(wrap("Stacked Amino Acid Distribution of " + sample, 40))
            ax.set_title(title, pad=12, **font_settings)
        else:
            title = "\n".join(wrap("Stacked Nucleotide Distribution of " + sample, 40))
            ax.set_title(title, pad=12, **font_settings)
        ax.set_xticklabels([*range(region[0], region[1] + 1)], rotation=0, ha="center")
        font_settings["fontsize"] = original_fontsize
