from ExpoSeq.tidy_data.tidy_rel_sequence_abun import cleaning
def relative_sequence_abundance(ax, sequencing_report, samples,max_levenshtein_distance,length_filter,batch, font_settings, legend_settings):
    all_samples = cleaning(sequencing_report,
                            max_levenshtein_distance,
                            samples,
                            length_filter,
                            batch)
    x_vals = all_samples.index.to_list()
    bar_width = 0.25
    all_samples_columns = all_samples.columns.to_list()
    for i in range(all_samples.shape[1]):
        y = all_samples.iloc[:, i]
        x_pos = [j + i * bar_width for j in range(len(x_vals))]
        ax.bar(x_pos,
               y,
               width=0.25,
               align='center',
               label=all_samples_columns[i])
    ax.set_xticks(range(len(x_vals)),
                  x_vals,
                  rotation = 45,
                  ha = "right")
    ax.set_ylabel("Frequency in samples with Levenshtein Distance " + str(max_levenshtein_distance),
                  **font_settings)
    ax.set_xlabel("Sequences",
                  **font_settings)
    ax.legend(**legend_settings)
    original_fontsize = font_settings["fontsize"]
    font_settings["fontsize"] = 22
    ax.set_title("Sequence Abundance of given samples",
                 pad = 12,
                 **font_settings)
    font_settings["fontsize"] = original_fontsize
