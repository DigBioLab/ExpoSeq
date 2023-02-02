from python_scripts.tidy_data.tidy_stacked_aa_distr import cleaning
import matplotlib.pyplot as plt
def stacked_aa_distr(ax, sequencing_report, sample, region, protein, font_settings, legend_settings):
    aa_distribution = cleaning(sequencing_report, sample, region, protein)
    aa_distribution.plot(kind='bar', stacked=True, ax = ax)
    ax.set_xlabel('Position on amino acid sequence')
    if protein == True:
        ax.set_ylabel('Relatvie Proportion of Amino Acid', **font_settings)
    else:
        ax.set_ylabel('Relatvie Proportion of Nucleotide', **font_settings)
    plt.legend(**legend_settings)