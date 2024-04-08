import re
from transformers import AutoTokenizer, AutoModel
import numpy as np
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
import seaborn as sns
from textwrap import wrap

# pip install rjieba
class PrepareData:
    @staticmethod
    def check_input(seq):
        assert type(seq) == list
        for single_seq in seq:
            for aa in single_seq:
                assert aa in ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "V", "W", "Y", "_", "*" ], f"{seq} does not contain amino acids"
    @staticmethod
    def prep_sequence(seq):
        assert type(seq) == list, "Input must be a list of strings"
        sequences = [" ".join(list(re.sub(r"[UZOB*_]", "X", sequence))) for sequence in seq]
        return sequences
    @staticmethod
    def prep_run_model(sequences, model_name):
        tokenizer = AutoTokenizer.from_pretrained(model_name)
        model = AutoModel.from_pretrained(model_name, output_attentions=True)
        inputs = tokenizer(sequences, return_tensors="pt")
        input_ids = inputs["input_ids"]
        tokens = [tokenizer.decode(id_) for id_ in input_ids[0]]
        outputs = model(input_ids)
        return outputs, tokens, model

    def tidy(self, seq:list, model:str, method = "outlier"):
        self.check_input(seq)
        sequences = self.prep_sequence(seq)
        outputs, tokens, model = self.prep_run_model(sequences, model_name=model)
        FindClusters = Find_head_clusters(outputs, model)
        FindClusters.remove_outliers_and_average()
        average_attention_map = FindClusters.non_outlier_attention_maps
        return average_attention_map

class Find_head_clusters:
    def __init__(self, output, model):
        self.model = model
        self.attentions = self.get_attentions(output)
        self.cluster_assignments = None  # Stores cluster assignments
        self.non_outlier_attention_maps = None

    def get_attentions(self, output):
        return output.attentions

    @staticmethod
    def flatten_attentions(attentions):
        flattened_attention_maps = []
        # Extract and flatten the attention maps, excluding CLS and SEP tokens
        for layer_attention in attentions:
            # layer_attention shape is [1, 16, seq_len, seq_len]
            # we remove the batch dimension and the first/last tokens (CLS and SEP)
            for head_attention in layer_attention.squeeze(0):
                # Exclude the first and last row/column corresponding to CLS and SEP
                attention_without_cls_sep = head_attention[1:-1, 1:-1]
                flattened_attention_maps.append(attention_without_cls_sep.detach().flatten().numpy())
        # Convert the list to a numpy array for further processing
        flattened_attention_maps_array = np.array(flattened_attention_maps)
        return flattened_attention_maps_array

    def pca_and_cluster_heads(self, n_components=4):
        flattened_attention_maps_array = self.flatten_attentions(self.attentions)
        pca = PCA(n_components=n_components)
        principal_components = pca.fit_transform(flattened_attention_maps_array)
        silhouette_scores = self.calc_silhouette_score(principal_components)
        n_clusters = silhouette_scores.index(max(silhouette_scores)) + 2
        kmeans = KMeans(n_clusters=n_clusters, random_state=42)
        self.cluster_assignments = kmeans.fit_predict(principal_components)

    def remove_outliers_and_average(self):
        # Flatten and calculate variance for each head's attention map
        flattened_attention_maps_array = self.flatten_attentions(self.attentions)
        attention_variances = np.var(flattened_attention_maps_array, axis=1)

        # Identify outliers using the IQR method
        Q1, Q3 = np.percentile(attention_variances, [25, 75])
        IQR = Q3 - Q1
        outlier_indices = np.where((attention_variances < (Q1 - 1.5 * IQR)) |
                                   (attention_variances > (Q3 + 1.5 * IQR)))[0]

        # Calculate non-outlier indices, ensuring they are within the valid range
        non_outlier_indices = [i for i in range(flattened_attention_maps_array.shape[0]) if i not in outlier_indices]

        # Collect and average non-outlier attention maps, excluding [CLS] and [SEP]
        non_outlier_attention_maps = []
        for i in non_outlier_indices:
            layer_index = i // self.model.config.num_attention_heads  # Correct division to get layer index
            head_index = i % self.model.config.num_attention_heads  # Correct modulo to get head index
            # Validate indices are within the expected range
            if layer_index < self.model.config.num_hidden_layers and head_index < self.model.config.num_attention_heads:
                attention_map = self.attentions[layer_index][0, head_index, 1:-1, 1:-1].detach().numpy() # this is only for the first sequence batch = 0
                non_outlier_attention_maps.append(attention_map)

        # Average the non-outlier attention maps, if any are collected
        if non_outlier_attention_maps:
            self.non_outlier_attention_maps = np.mean(non_outlier_attention_maps, axis=0)
        else:
            print("No non-outlier attention maps found. Check outlier detection criteria.")

    @staticmethod
    def calc_silhouette_score(principal_components, lower_range = 2, higher_range = 10):
        assert lower_range >= 2
        K = range(lower_range, higher_range)
        from sklearn.metrics import silhouette_score
        silhouette_scores = []
        for k in K:  # Assuming K is defined as in the Elbow Method
            kmeanModel = KMeans(n_clusters=k, random_state=42)
            labels = kmeanModel.fit_predict(principal_components)
            silhouette_scores.append(silhouette_score(principal_components, labels))
        return silhouette_scores

    def get_cluster_attention_maps(self):
        # Initialize a dict to hold attention maps for each cluster
        cluster_attention_maps = {cluster: [] for cluster in set(self.cluster_assignments)}

        # Iterate through each attention map and its cluster assignment
        for i, cluster in enumerate(self.cluster_assignments):
            # Retrieve the original attention map (unflattened)
            layer = i // 16  # Assuming 16 heads per layer
            head = i % 16
            attention_map = self.attentions[layer][0, head].detach().numpy()
            cluster_attention_maps[cluster].append(attention_map)

        # Average the attention maps within each cluster
        for cluster in cluster_attention_maps:
            cluster_attention_maps[cluster] = np.mean(cluster_attention_maps[cluster], axis=0)

        return cluster_attention_maps


class AverageAttentionHeatmap:
    def __init__(self, sequences:list, model:str,cmap = "Purples", annotate_cells = True,  fmt = ".2f", ax = None, font_settings=  {}):
        PrepData = PrepareData()
        self.matrix = PrepData.tidy(sequences, model)
        self.ax = ax
        if self.ax:
            self.plot_heatmap(len(sequences[0]), cmap, annotate_cells, fmt)
            if font_settings != {}:
                self.add_labels(font_settings)
                self.add_header(font_settings)

    def plot_heatmap(self, len_sequence, cmap, annotate_cells, fmt):
        annot_kws = {"size": annot_kws}
        pos_range = list(range(1, len_sequence ))
        sns.heatmap(self.matrix, xticklabels=pos_range, yticklabels=pos_range,
                    cmap = cmap,ax = self.ax, annot_kws=annot_kws, annot=annotate_cells, fmt = fmt)
        
    def add_labels(self, font_settings):
        self.ax.set_xlabel("Total sampled sequences", **font_settings)
        self.ax.set_ylabel("Total Unique Sequences", **font_settings)
        
    def add_header(self, font_settings):
        original_fontsize = font_settings["fontsize"]
        font_settings["fontsize"] = 22
        title = "Average Attention for chosen sequences in the model."
        title = "\n".join(wrap(title, 40))
        self.ax.set_title(title,pad = 12, **font_settings)
        font_settings["fontsize"] = original_fontsize



#import pandas as pd

#umap_data = pd.read_csv(r"C:\Users\nilsh\my_projects\ExpoSeq\tmp_test\embedding_umap_test.csv")
#reduced_umap = umap_data[umap_data["cluster_id"] == 1]
#reduced_umap['length'] = reduced_umap['sequences'].apply(len)
#most_common_length = reduced_umap['length'].mode().iloc[0]
#filtered_df = reduced_umap[reduced_umap['length'] == most_common_length]
#sequences = filtered_df["sequences"].to_list()
#Start = AverageAttentionHeatmap(sequences, "alchemab/antiberta2")
#import matplotlib
#matplotlib.use("Qt5Agg")
#Start.plot_heatmap()
