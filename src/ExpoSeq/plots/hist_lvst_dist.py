import editdistance
from scipy.cluster.hierarchy import dendrogram, linkage
from matplotlib import pyplot as plt
import numpy as np
from scipy.spatial.distance import squareform
from ExpoSeq.tidy_data.tidy_hist_lvst_dist import create_distance_matrix, get_clustered_sequences




def levenshtein_hist(ax, sequencing_report, sample, max_dist, batch_size, font_settings):
    sample_report = sequencing_report[sequencing_report["Experiment"] == sample] ## insert test if sample not found
    report = sample_report.head(batch_size)
    aa = list(report["aaSeqCDR3"])
    aa_clustered = get_clustered_sequences(aa)

    # Create the distance matrix using the filtered list of sequences
    levenshtein_distance_matrix = create_distance_matrix(aa_clustered, max_dist)

    # Convert to condensed form and create the dendrogram
    condensed_matrix = squareform(levenshtein_distance_matrix, checks=False)
    linked = linkage(condensed_matrix, 'single')

    dendro = dendrogram(linked,
                        orientation='right',
                        distance_sort='descending',
                        show_leaf_counts=True,
                        labels=aa_clustered,
                        ax = ax)
    ax.set_xlabel("Levenshtein Distance", **font_settings)
    ax.set_ylabel("Sequences", **font_settings)
    plt.tight_layout()

