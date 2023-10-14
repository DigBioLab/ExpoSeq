from scipy.cluster.hierarchy import dendrogram, linkage
from matplotlib import pyplot as plt
from scipy.spatial.distance import squareform
import pandas as pd
import warnings
from ..tidy_data.tidy_dendro import create_distance_matrix, get_clustered_sequences, label_bind_seqs, sort_binding_seqs




def levenshtein_dend(ax, sequencing_report, sample, batch_size,max_cluster_dist, font_settings, region_string):
    sample_report = sequencing_report[sequencing_report["Experiment"] == sample] ## insert test if sample not found
    report = sample_report.head(batch_size)
    aa = list(report[region_string])
    aa_clustered = get_clustered_sequences(aa, max_cluster_dist)
    # Create the distance matrix using the filtered list of sequences
    levenshtein_distance_matrix = create_distance_matrix(aa_clustered)
    if levenshtein_distance_matrix.shape[0] == 0:
        print("0 matches were found please increase the batch size or the levenshtein distance")
        return
    # Convert to condensed form and create the dendrogram
    condensed_matrix = squareform(levenshtein_distance_matrix, checks=False)
    linked = linkage(condensed_matrix, 'single')
    dendrogram(linked,
                orientation='right',
                distance_sort='descending',
                show_leaf_counts=True,
                labels=aa_clustered,
                ax = ax
                )
    ax.set_xlabel("Levenshtein Distance", **font_settings)
    ax.set_ylabel("Sequences", **font_settings)
    title = "Levenshtein Distance between sequences in " + sample 
    ax.set_title(title,pad = 12, **font_settings)
    plt.tight_layout()


def dendo_binding(fig, sequencing_report,binding_data, sample,antigens, batch_size,max_cluster_dist,font_settings, region_string, ascending ):
    ax = fig.gca()

    sample_report = sequencing_report[sequencing_report["Experiment"] == sample] ## insert test if sample not found
    report = sample_report.head(batch_size)
    aa = report[region_string]
    aa = pd.DataFrame(aa)
    pref_columns = antigens + [region_string]
    b_data = binding_data[pref_columns]
    
    mix = pd.concat([aa, b_data])
    mix = mix.fillna(0)
    mix = mix.reset_index()
    
    aa_all = mix[region_string]
    aa_clustered = get_clustered_sequences(aa_all, max_cluster_dist)
    if len(aa_clustered) > 0:
        warnings.warn("More than 30 sequences with Levenshtein distance < " + str(max_cluster_dist) + " found. The resulting plot could be disordered. To change that please reduce the batch size, max_cluster_dist or you can adjust the string size of the sequences on the y axis with: (1) ax = plot.ax and (2) ax.tick_params(axis='y', labelsize=your_desired_size)")
    
    levenshtein_distance_matrix = create_distance_matrix(aa_clustered)
    condensed_matrix = squareform(levenshtein_distance_matrix, checks=False)
    linked = linkage(condensed_matrix, 'single')
    
    aa_clustered, binding_seqs, seq_val = label_bind_seqs(mix,region_string, aa_clustered, antigens, )
    binding_seqs_sorted, binding_values_sorted = sort_binding_seqs(binding_seqs, seq_val, ascending)


    dendrogram(linked,
            orientation='right',
            distance_sort='descending',
            show_leaf_counts=True,
            labels=aa_clustered,
             ax = ax
            )
    ax = plt.gca()
    ax.set_xlabel("Levenshtein Distance", **font_settings)
    ax.set_ylabel("Sequences", **font_settings)
    title = "Levenshtein Distance between sequences in " + sample
    
    ax.set_title(title,pad = 12, **font_settings)
    fig.show()
    fig2 = plt.figure(2)
    ax2 = fig2.gca()
    bars = ax2.barh(binding_seqs_sorted, binding_values_sorted)
    ax2.set_ylabel('Sequences with binding data', **font_settings)
    ax2.set_xlabel('Binding Value', **font_settings)
    fig2.tight_layout()
    fig.tight_layout()
    return fig2

    
