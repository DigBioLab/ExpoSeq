from ExpoSeq.tidy_data.tidy_stacked_aa_distr import cleaning
import matplotlib.pyplot as plt

def stacked_aa_distr(ax, sequencing_report, sample, region, protein, font_settings, legend_settings):
    aa_distribution = cleaning(sequencing_report, sample, region, protein)
    aa_distribution.plot(kind='bar', stacked=True, ax = ax)
    ax.set_xlabel('Position on amino acid sequence', **font_settings)
    if protein == True:
        ax.set_ylabel('Relatvie Proportion of Amino Acid', **font_settings)
    else:
        ax.set_ylabel('Relatvie Proportion of Nucleotide', **font_settings)

#    ax.set_xticks(rotation = 360)
    ax.legend(**legend_settings)
    original_fontsize = font_settings["fontsize"]
    font_settings["fontsize"] = 22
    if protein == True:
        ax.set_title("Distribution of Aminoacids per Position for your sequences",pad = 12, **font_settings)
    else:
        ax.set_title("Distribution of Nucleotides per Position for your sequences",pad = 12, **font_settings)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=0, ha='center')
    font_settings["fontsize"] = original_fontsize