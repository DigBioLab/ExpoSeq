from scipy.cluster.hierarchy import dendrogram, linkage
from matplotlib import pyplot as plt
from scipy.spatial.distance import squareform
import pandas as pd
import warnings
from ExpoSeq.tidy_data.tidy_dendro import create_distance_matrix, get_clustered_sequences, label_bind_seqs, sort_binding_seqs
from textwrap import wrap
import editdistance
import numpy as np


class LevenshteinDend:
    def __init__(self, sequencing_report, region_of_interest,sample, font_settings = {}, batch_size = 300, max_cluster_dist = 2, ax = None) -> None:
        linked, aa_clustered = self.tidy(sequencing_report, sample, batch_size, region_of_interest,max_cluster_dist)
        if ax != None:
            self.plot(linked, ax, aa_clustered, font_settings, sample)
        
        
    @staticmethod
    def get_clustered_sequences(aa, max_cluster_dist):
        clustered_sequences = set()

        for i in range(len(aa)):
            for j in range(i+1, len(aa)):
                if editdistance.distance(aa[i], aa[j]) < max_cluster_dist:
                    clustered_sequences.add(aa[i])
                    clustered_sequences.add(aa[j])
        clustered_sequences = set(clustered_sequences)
        return list(clustered_sequences)
    
    
    @staticmethod        
    def create_distance_matrix(aa):
        num_sequences = len(aa)
        distance_matrix = np.zeros((num_sequences, num_sequences))  # initialize with zeros
        for i in range(num_sequences):
            for j in range(i+1, num_sequences):
                lev_distance = editdistance.distance(aa[i], aa[j])
                distance_matrix[i, j] = lev_distance
                distance_matrix[j, i] = lev_distance  # The distance matrix is symmetric
        return distance_matrix

    
    def tidy(self, sequencing_report, sample, batch_size, region_string, max_cluster_dist):
        """
        linked: linkage matrix which contains per sequence 4 values.
        aa_clustered: List of sequences which length is the same as the number of arrays in the linkage matrix. 
        """
        sample_report = sequencing_report[sequencing_report["Experiment"] == sample] ## insert test if sample not found
        report = sample_report.head(batch_size)
        aa = list(report[region_string])
        aa_clustered = self.get_clustered_sequences(aa, max_cluster_dist)
        # Create the distance matrix using the filtered list of sequences
        levenshtein_distance_matrix = self.create_distance_matrix(aa_clustered)
        if levenshtein_distance_matrix.shape[0] == 0: # user input test
            raise ValueError("0 matches were found please increase the batch size or the levenshtein distance")
        else:
            # Convert to condensed form and create the dendrogram
            condensed_matrix = squareform(levenshtein_distance_matrix, checks=False)
            linked = linkage(condensed_matrix, 'single') # nearest point algorithm, n-1 matrix returned
        return linked, aa_clustered
    
    def plot(self, linked, ax, aa_clustered, font_settings, sample):
        dendrogram(linked,
                orientation='right',
                distance_sort='descending',
                show_leaf_counts=True,
                labels=aa_clustered,
                ax = ax
                )
        ax.set_xlabel("Levenshtein Distance", **font_settings)
        ax.set_ylabel("Sequences", **font_settings)
        title = "\n".join(wrap("Levenshtein Distance between sequences in " + sample, 40))
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
    if len(binding_seqs) == 0:
        print("No matches between your binding data and your sequences were found. Please increase the batch size or change the antigen.")
        fig2 = False
        return fig2
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
    title = "\n".join(wrap("Levenshtein Distance between sequences in " + sample, 40))
    
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

    
