import pandas as pd
import numpy as np
import editdistance
from scipy.spatial.distance import squareform
import warnings
from scipy.cluster.hierarchy import dendrogram, linkage
from textwrap import wrap
from matplotlib import pyplot as plt


class PrepareData:
    @staticmethod
    def check_input(sample, batch_size, max_cluster_dist, antigens):
        assert (
            type(sample) == str
        ), "You have to give a string as input for the sample"  # only one sample
        assert (
            type(max_cluster_dist) == int
        ), "You have to give an integer as input for the maximum levenshtein distance"
        assert (
            max_cluster_dist != 0
        ), "Plot does not make sense if you do not investigate at least SNPs"
        assert (
            type(batch_size) == int
        ), "You have to give an integer as input for the batch size"
        assert batch_size > 1
        assert antigens != None, "Please provide an antigen name"

    @staticmethod
    def tidy_data(
        sequencing_report,
        sample,
        batch_size,
        region_of_interest,
        antigens,
        binding_data,
    ):
        sample_report = sequencing_report[
            sequencing_report["Experiment"] == sample
        ]  ## insert test if sample not found
        report = sample_report.head(batch_size)
        aa = report[region_of_interest]
        aa = pd.DataFrame(aa)
        pref_columns = antigens + [region_of_interest]
        for i in antigens:
            assert i in binding_data.columns.tolist(), f"{i} not in your binding data"
        assert region_of_interest in binding_data.columns.tolist()
        b_data = binding_data[pref_columns]
        mix = aa.merge(b_data, how="outer", on=region_of_interest)
        mix = mix.fillna(0)
        mix = mix.reset_index()
        aa_all = mix[region_of_interest]
        return mix, aa_all

    @staticmethod
    def create_distance_matrix(aa):
        num_sequences = len(aa)
        distance_matrix = np.zeros(
            (num_sequences, num_sequences)
        )  # initialize with zeros

        for i in range(num_sequences):
            for j in range(i + 1, num_sequences):
                lev_distance = editdistance.distance(aa[i], aa[j])
                distance_matrix[i, j] = lev_distance
                distance_matrix[j, i] = lev_distance  # The distance matrix is symmetric

        return distance_matrix

    @staticmethod
    def get_clustered_sequences(aa, max_cluster_dist):
        clustered_sequences = set()

        for i in range(len(aa)):
            for j in range(i + 1, len(aa)):
                if editdistance.distance(aa[i], aa[j]) < max_cluster_dist:
                    clustered_sequences.add(aa[i])
                    clustered_sequences.add(aa[j])

        return list(clustered_sequences)

    @staticmethod
    def label_bind_seqs(mix, region_string, aa_clustered, antigens):
        values = mix.loc[
            mix[region_string].isin(aa_clustered), antigens
        ]  # get the binding values for the sequences which were clustered
        # values = pd.DataFrame(values, columns = ["binding"])
        filtered_df = values[antigens][values[antigens] > 1].dropna(
            how="all"
        )  # filter out the binidng values which are 0
        indices = filtered_df.index.tolist()
        assert (
            len(indices) != 0
        ), "The cluster does not contain binding values. Please increase the lev. distance or the batch size."
        key_sequences = []
        key_values = []
        assert mix.columns.to_list() == ["index", "aaSeqCDR3"] + antigens
        for i in indices:
            key_sequences.append(mix.iloc[i][region_string])
            max_value = max(mix.iloc[i][column] for column in antigens)
            key_values.append(max_value)
        seq_val = pd.DataFrame([key_sequences, key_values]).T
        seq_val.columns = [
            "binding_seqs",
            "binding_values",
        ]  # seq_val object contains only the sequences from sanger - this is important for labelling
        seq_val = seq_val.sort_values(by="binding_values", ascending=False)
        seq_val.reset_index(inplace=True)
        binding_seqs = seq_val.binding_seqs.to_list()
        for label in range(len(aa_clustered)):
            label_str = aa_clustered[label]
            if label_str in binding_seqs:
                pos_label = binding_seqs.index(label_str)
                label_str = "(" + str(pos_label + 1) + ") " + label_str
                binding_seqs[pos_label] = label_str
            aa_clustered[label] = label_str
        return aa_clustered, binding_seqs, seq_val

    def process(
        self,
        sequencing_report,
        sample,
        batch_size,
        region_of_interest,
        antigens,
        binding_data,
        max_cluster_dist,
        sample_column_name="Experiment",
    ):
        assert sample in list(
            sequencing_report[sample_column_name].unique()
        ), "The provided sample name does not exist"
        self.check_input(sample, batch_size, max_cluster_dist, antigens)
        mix, aa_all = self.tidy_data(
            sequencing_report,
            sample,
            batch_size,
            region_of_interest,
            antigens,
            binding_data,
        )
        aa_clustered = self.get_clustered_sequences(aa_all, max_cluster_dist)
        if len(aa_clustered) > 0:
            warnings.warn(
                "More than 30 sequences with Levenshtein distance < "
                + str(max_cluster_dist)
                + " found. The resulting plot could be disordered. To change that please reduce the batch size, max_cluster_dist or you can adjust the string size of the sequences on the y axis with: (1) ax = plot.ax and (2) ax.tick_params(axis='y', labelsize=your_desired_size)"
            )
        levenshtein_distance_matrix = self.create_distance_matrix(aa_clustered)
        condensed_matrix = squareform(
            levenshtein_distance_matrix, checks=False
        )  # linearize, now only 1 dim array
        linked = linkage(condensed_matrix, "single")
        aa_clustered, binding_seqs, seq_val = self.label_bind_seqs(
            mix, region_of_interest, aa_clustered, antigens
        )
        return linked, aa_clustered, binding_seqs, seq_val


class DendroBind:
    def __init__(
        self,
        sequencing_report,
        sample,
        region_of_interest,
        antigens,
        batch_size,
        max_cluster_dist,
        binding_data,
        ascending=True,
        fig=None,
        font_settings={},
    ) -> None:
        linked, aa_clustered, binding_seqs, seq_val = PrepareData().process(
            sequencing_report,
            sample,
            batch_size,
            region_of_interest,
            antigens,
            binding_data,
            max_cluster_dist,
        )
        if fig != None:
            self.ax = fig.gca()
            self.fig2 = plt.figure(2)
            self.ax2 = self.fig2.gca()
            self.create_dendrogram(linked, aa_clustered)
            self.create_bar_plot(binding_seqs, seq_val, ascending)
        if font_settings != {}:
            self.set_dendrogram_labels(font_settings, sample)
            self.set_barplot_labels(font_settings)

    def create_dendrogram(self, linked, aa_clustered):
        dendrogram(
            linked,
            orientation="right",
            distance_sort="descending",
            show_leaf_counts=True,
            labels=aa_clustered,
            ax=self.ax,
        )

    def set_dendrogram_labels(
        self,
        font_settings,
        sample,
    ):
        self.ax.set_xlabel("Levenshtein Distance", **font_settings)
        self.ax.set_ylabel("Sequences", **font_settings)
        title = "\n".join(
            wrap("Levenshtein Distance between sequences in " + sample, 40)
        )

        self.ax.set_title(title, pad=12, **font_settings)

    @staticmethod
    def sort_binding_seqs(binding_seqs, seq_val, ascending=False):

        if ascending:

            zipped = list(zip(binding_seqs, seq_val.binding_values))
            # Sort by binding_values
            zipped_sorted = sorted(zipped, key=lambda x: x[1])

            # Unzip the sorted list
            binding_seqs_sorted, binding_values_sorted = zip(*zipped_sorted)
        else:
            binding_seqs_sorted = binding_seqs
            binding_values_sorted = seq_val.binding_values.to_list()
        return binding_seqs_sorted, binding_values_sorted

    def create_bar_plot(self, binding_seqs, seq_val, ascending):
        binding_seqs_sorted, binding_values_sorted = self.sort_binding_seqs(
            binding_seqs, seq_val, ascending
        )
        bars = self.ax2.barh(binding_seqs_sorted, binding_values_sorted)

    def set_barplot_labels(self, font_settings):
        self.ax2.set_ylabel("Sequences with binding data", **font_settings)
        self.ax2.set_xlabel("Binding Value", **font_settings)
