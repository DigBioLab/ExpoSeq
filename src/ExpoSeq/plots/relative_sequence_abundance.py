from ExpoSeq.tidy_data.tidy_rel_sequence_abun import cleaning
import numpy as np
from ExpoSeq.plots.layout_finder import best_layout
from textwrap import wrap
    

def relative_sequence_abundance(ax, sequencing_report, samples,max_levenshtein_distance,length_filter,batch,region_string, font_settings, legend_settings):
    all_samples = cleaning(sequencing_report,
                            max_levenshtein_distance,
                            samples,
                            length_filter,
                            region_string,
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
    title = "\n".join(wrap("Sequence Abundance of given samples", 40))
    ax.set_title(title,
                 pad = 12,
                 **font_settings)
    font_settings["fontsize"] = original_fontsize
    
def relative_sequence_abundance_all(fig, sequencing_report, samples, font_settings, region_string):
    if samples == "all":
        unique_experiments = sequencing_report["Experiment"].unique()
        unique_experiments = np.sort(unique_experiments)
    else:
        sequencing_report = sequencing_report[sequencing_report['Experiment'].isin(samples)]
        unique_experiments = sequencing_report["Experiment"].unique()
        unique_experiments = np.sort(unique_experiments)
        
    Tot = unique_experiments.shape[0]
    Rows, Cols = best_layout(Tot)
    Position = range(1, Tot + 1)
    n = 0
    
    for experiment in unique_experiments:
        batch = sequencing_report[sequencing_report["Experiment"] == experiment]
        fraction = batch["cloneFraction"]
        clones = batch["cloneId"]
        ax = fig.add_subplot(Rows, Cols, Position[n])
        
        ax.bar(clones, fraction)  # Or whatever you want in the subplot
       # ax.set_xticks(range(0, max_length + 1, 1), range(0, max_length + 1, 1))
        ax.title.set_text(experiment)
        ax.title.set_size(10)

        adapted_fontsize = 10 - int(Cols) + 2
        font_settings["fontsize"] = adapted_fontsize
        ax.set_ylabel("Clone Fraction",
                      **font_settings)  # Y label
        ax.set_xlabel('Clones',
                      **font_settings)  # X label

        ax.set_ylim(0, np.max(np.array(fraction) + np.max(np.array(fraction) * 0.1)))
        n += 1
    fig.suptitle("Clone Fraction of given Samples")

