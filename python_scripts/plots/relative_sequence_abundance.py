from python_scripts.tidy_data.tidy_rel_sequence_abun import cleaning

def relative_sequence_abundance(ax, sequencing_report, samples,max_levenshtein_distance,length_filter,batch, font_settings):
    all_samples = cleaning(sequencing_report,
                                             max_levenshtein_distance,
                                             samples,
                                             length_filter,
                                             batch)
    x_vals = all_samples.index.to_list()
    bar_width = 0.25
    for i in range(all_samples.shape[1]):
        y = all_samples.iloc[:, i]
        x_pos = [j + i * bar_width for j in range(len(x_vals))]
        ax.bar(x_pos, y, width=0.25, align='center', label='bar_' + str(i+1))
    ax.xticks(range(len(x_vals)), x_vals, rotation = 45, ha = "right")
    ax.set_ylabel("Frequency in samples with Levenshtein Distance " + str(max_levenshtein_distance), **font_settings)
    ax.set_ylabel("Sequences", **font_settings)