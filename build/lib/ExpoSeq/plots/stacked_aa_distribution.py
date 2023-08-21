from ExpoSeq.tidy_data.tidy_stacked_aa_distr import cleaning
import matplotlib.pyplot as plt

def stacked_aa_distr(ax, sequencing_report, sample, region, protein, font_settings,region_string):

    aa_distribution = cleaning(sequencing_report, sample, region, region_string)
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
        ax.set_title("Distribution of Aminoacids per Position for your sequences",pad = 12, **font_settings)
    else:
        ax.set_title("Distribution of Nucleotides per Position for your sequences",pad = 12, **font_settings)
    ax.set_xticklabels([*range(region[0], region[1]+1)], rotation=0, ha='center')
    font_settings["fontsize"] = original_fontsize