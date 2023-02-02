from python_scripts.tidy_data.tidy_rel_sequence_abun import tidy_levenshtein_fractions

def relative_sequence_abundance(ax, sequencing_report, experiments, max_levenshtein_distance = 0,length_filter = 5,batch = 3000):
    all_samples = tidy_levenshtein_fractions(sequencing_report,
                                             max_levenshtein_distance,
                                             experiments,
                                             length_filter,
                                             batch)
    x_vals = all_samples.index.to_list()
    bar_width = 0.25
    for i in range(all_samples.shape[1]):
        y = all_samples.iloc[:, i]
        x_pos = [j + i * bar_width for j in range(len(x_vals))]
        ax.bar(x_pos, y, width=0.25, align='center', label='bar_' + str(i+1))
    ax.xticks(range(len(x_vals)), x_vals, rotation = 45, ha = "right")
    ax.ylabel("Frequency in samples with Levenshtein Distance " + max_levenshtein_distance)