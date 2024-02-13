from scipy.cluster.hierarchy import dendrogram, linkage
from matplotlib import pyplot as plt
from scipy.spatial.distance import squareform
import pandas as pd
import warnings
from textwrap import wrap
import editdistance
import numpy as np


class PrepareData:
    @staticmethod
    def create_distance_matrix(aa, max_dist=20):  # max_dist is a large finite value
        num_sequences = len(aa)
        distance_matrix = np.full((num_sequences, num_sequences), max_dist)  # initialize with max_dist

        for i in range(num_sequences):
            for j in range(i+1, num_sequences):
                lev_distance = editdistance.distance(aa[i], aa[j])
                if lev_distance < 2:
                    distance_matrix[i, j] = lev_distance
                    distance_matrix[j, i] = lev_distance  # The distance matrix is symmetric

        return distance_matrix
    @staticmethod
    def get_clustered_sequences(aa, max_cluster_dist):
        clustered_sequences = set()

        for i in range(len(aa)):
            for j in range(i+1, len(aa)):
                if editdistance.distance(aa[i], aa[j]) < max_cluster_dist:
                    clustered_sequences.add(aa[i])
                    clustered_sequences.add(aa[j])

        return list(clustered_sequences)
    
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



class LevenshteinDend:
    def __init__(self, sequencing_report, region_of_interest,sample, font_settings = {}, batch_size = 300, max_cluster_dist = 2, ax = None) -> None:
        linked, aa_clustered = PrepareData().tidy(sequencing_report, sample, batch_size, region_of_interest, max_cluster_dist)
        if ax != None:
            self.plot(linked, ax, aa_clustered, font_settings, sample)
        
    
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



