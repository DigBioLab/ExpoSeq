import editdistance
from ..tidy_data.tidy_sample_cluster import cleaning
import networkx as nx
import math
import seaborn as sns
import community
import matplotlib.pyplot as plt


class ClusterExperiment():
    def __init__(self, sequencing_report, ax, region_string, summed_clonefraction, max_num_reads, edge_color = "peachpuff", max_weight_lines = 100):
        report = sequencing_report
        self.region_string = region_string
        report = cleaning(report, summed_clonefraction=summed_clonefraction, max_num_reads=max_num_reads)
        self.G = nx.Graph()
        distances = self.nodes_and_edges(report, max_weight_lines)
        nodesize = self.customizations(report)
        sample_colors, node_colors = self.create_colors(report)
        self.draw_plot(nodesize, node_colors, ax, edge_color)
        self.legend(report, sample_colors)
        
    def nodes_and_edges(self, report, max_weight_lines = 100):
        distances = []
        for sample, group in report.groupby('Experiment'):
            sequences = group[self.region_string].tolist()
            for i in range(len(sequences)):
                seq1 = sequences[i]
                # Add node for seq1
                self.G.add_node(seq1, sample=sample)

                for j in range(i + 1, len(sequences)):
                    seq2 = sequences[j]
                    distance = editdistance.distance(seq1, seq2)


                    self.G.add_node(seq2, sample=sample)
                    distances.append((sample, seq1, seq2, distance))
                    # Add nodes for seq2 and connect them with an edge weighted by the Levenshtein distance
                    if distance != 0:
                        if distance > 5:
                            self.G.add_edge(seq1, seq2, weight = int(max_weight_lines * 1/(distance**2)))
                        
                    else:
                        self.G.add_edge(seq1, seq2, weight= max_weight_lines)
        return distances
    
    def customizations(self, report):
        nodesize = []
        label_numbers = {}
        label_sequences = {}


        for index, g in enumerate(self.G):
            for i in report[self.region_string]:
                if i == g:
                    node = math.sqrt(report.loc[report[self.region_string] == i, "cloneFraction"].values[0])
                    nodesize.append(node)
                    if node * 500 > 1 / 100 * 500:
                        label_numbers[g] = index
                        label_sequences[g] = g
                    else:
                        label_numbers[g] = ""
                        label_sequences[g] = ""
                    break
        for y, n in enumerate(nodesize):
            nodesize[y] = int(n * 500)
        max_fraction_node = max(label_numbers, key=label_numbers.get)
     #   label_numbers = {node: index if node == max_fraction_node else "" for node, index in label_numbers.items()}
      #  label_sequences = {node: seq if node == max_fraction_node else "" for node, seq in label_sequences.items()}

        return nodesize

    def create_colors(self, report):
        num_samples = len(report['Experiment'].unique())

        # Generate a color palette with num_samples distinct colors
        palette = sns.color_palette("hsv", num_samples)

        # Create a dictionary to map samples to colors
        sample_colors = {sample: color for sample, color in zip(report['Experiment'].unique(), palette)}
        node_colors = [sample_colors[self.G.nodes[n]['sample']] for n in self.G.nodes]
        return sample_colors, node_colors
        
    def draw_plot(self, nodesize, node_colors, ax, edge_color):
        G = self.G
        partition = community.best_partition(G)

        weight = 0.1
        
        pos = nx.spring_layout(G)
        nx.draw_networkx_nodes(self.G, pos, node_size=nodesize, node_color = node_colors, ax = ax)
        nx.draw_networkx_edges(self.G, pos, alpha=0.1, edge_color=edge_color, ax = ax)
        
    def legend(self, report, sample_colors):
        unique_samples = report['Experiment'].unique()

                # Create a list of colors corresponding to the samples
        colors = [sample_colors[sample] for sample in unique_samples]

        # Create a legend mapping sample names to colors
        legend = {sample: color for sample, color in zip(unique_samples, colors)}

        # Add a legend to the plot
        plt.legend(handles=[plt.Line2D([0], [0], marker='o', color='w', label=sample,
                                    markersize=10, markerfacecolor=color)
                            for sample, color in legend.items()],
                loc =  'center left', bbox_to_anchor =  (1, 0.5))
