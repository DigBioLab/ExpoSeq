import editdistance
from scipy.cluster.hierarchy import dendrogram, linkage
from matplotlib import pyplot as plt
import numpy as np
from scipy.spatial.distance import squareform
#from ExpoSeq.tidy_data.tidy_hist_lvst_dist import create_distance_matrix, get_clustered_sequences
import pandas as pd

from matplotlib import collections
import numpy as np
import editdistance


def create_distance_matrix(aa):
    num_sequences = len(aa)
    distance_matrix = np.zeros((num_sequences, num_sequences))  # initialize with zeros

    for i in range(num_sequences):
        for j in range(i+1, num_sequences):
            lev_distance = editdistance.distance(aa[i], aa[j])
            distance_matrix[i, j] = lev_distance
            distance_matrix[j, i] = lev_distance  # The distance matrix is symmetric

    return distance_matrix

def get_clustered_sequences(aa, max_cluster_dist):
    clustered_sequences = set()

    for i in range(len(aa)):
        for j in range(i+1, len(aa)):
            if editdistance.distance(aa[i], aa[j]) < max_cluster_dist:
                clustered_sequences.add(aa[i])
                clustered_sequences.add(aa[j])

    return list(clustered_sequences)


def levenshtein_dend(ax, sequencing_report, sample, batch_size,max_cluster_dist, font_settings, region_of_interest):
    sample_report = sequencing_report[sequencing_report["Experiment"] == sample] ## insert test if sample not found
    report = sample_report.head(batch_size)
    aa = list(report[region_of_interest])
    aa_clustered = get_clustered_sequences(aa, max_cluster_dist)
    # Create the distance matrix using the filtered list of sequences
    levenshtein_distance_matrix = create_distance_matrix(aa_clustered)

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


def dendo_binding(ax, sequencing_report,binding_data, sample,antigens, batch_size,max_cluster_dist, scale_factor_lines, font_settings, region_of_interest ):
    sample_report = sequencing_report[sequencing_report["Experiment"] == sample] ## insert test if sample not found
    report = sample_report.head(batch_size)
    aa = report[region_of_interest]
    aa = pd.DataFrame(aa)
    pref_columns = antigens + [region_of_interest]
    b_data = binding_data[pref_columns]
    
    mix = pd.concat([aa, b_data])
    mix = mix.fillna(0)
    mix = mix.reset_index()
    
    aa_all = mix[region_of_interest]
    aa_clustered = get_clustered_sequences(aa_all, max_cluster_dist)
    levenshtein_distance_matrix = create_distance_matrix(aa_clustered)
    condensed_matrix = squareform(levenshtein_distance_matrix, checks=False)
    linked = linkage(condensed_matrix, 'single')
    values = mix.loc[mix[region_of_interest].isin(aa_clustered), antigens]
    #values = pd.DataFrame(values, columns = ["binding"])
    filtered_df = values[antigens][values[antigens] > 1].dropna(how='all')
    indices = filtered_df.index.tolist()
    key_sequences = []
    key_values = []
    for i in indices:
        key_sequences.append(mix.iloc[i][region_of_interest])
        max_value = max(mix.iloc[i][column] for column in antigens)
        key_values.append(max_value)
    
    if len(key_sequences)!=0:
        sequence_values = []
        for sequence in aa_clustered:
            value = 0
            for key_sequence, key_value in zip(key_sequences, key_values):
                lev_distance = editdistance.distance(sequence, key_sequence)
                if lev_distance > 0:
                    value += key_value / lev_distance  # add the ratio for each key sequence
                else:
                    value += key_value / 1
            sequence_values.append(value)
    else:
        sequence_values = [1] * len(aa_clustered)
    
    # Step 4: Normalize these values to suitable line widths.
    if sequence_values == []:
        sequence_values.append(1)
    max_value = max(sequence_values)
    normalized_values = [value / max_value for value in sequence_values]  # normalize to [0, 1]
    line_widths = [1 + 4 * value for value in normalized_values]  # scale to [1, 5] for line width
    # Step 5: Create the dendrogram and assign line widths
     # or any other value that works best in your case
    scaled_line_widths = [width * scale_factor_lines for width in line_widths]
    dendrogram(linked,
               orientation='right',
               distance_sort='descending',
               show_leaf_counts=True,
               labels=aa_clustered,
                ax = ax
               )
    ax = plt.gca()
    if len(key_sequences)!=0:
        line_collection = collections.LineCollection(ax.collections[0].get_segments())
        line_collection.set_linewidth(scaled_line_widths)  # set the line widths
    
        ax.collections[0].remove()
        ax.add_collection(line_collection)
    else:
        print("No sequence in your binidng data could be found that forms a cluster with a sequences in the sample given the levenshtein distance.")
    ax.set_xlabel("Levenshtein Distance", **font_settings)
    ax.set_ylabel("Sequences", **font_settings)
    title = "Levenshtein Distance between sequences in " + sample 
    ax.set_title(title,pad = 12, **font_settings)
  #  plt.tight_layout()