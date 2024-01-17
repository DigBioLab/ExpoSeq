import networkx as nx
import math
from textwrap import wrap
import numpy as np
from ..tidy_data.tidy_protbert_embedding import TransformerBased
import matplotlib.pyplot as plt
import community

class Network_Embedding:
    def __init__(self, ax, sequencing_report, list_experiments, batch_size,font_settings,cmap = plt.cm.RdYlBu, nodesize = None, threshold_distance = 1) -> None:
        Transformer = TransformerBased(choice = "protbert")
        sequences,sequences_filtered,  selected_rows = Transformer.filter_sequences(sequencing_report, batch_size=batch_size, experiments = list_experiments, region_of_interest="aaSeqCDR3", )
        sequences_filtered = sequences_filtered.to_list()
        sequences_list = Transformer.embedding_per_seq(sequences, normalize = True)
        
        self.G = nx.Graph()
        self.add_edges(sequences_list, sequences_filtered, threshold_distance)
        if nodesize == None:
            nodesize = self.add_nodesize(self.G, selected_rows)
        else:
            nodesize = 30
        self.generate_plot(self.G, nodesize, ax, cmap)
        title = "\n".join(wrap("t-SNE embedding of for given samples", 40))
        self.ax.set_title(title, pad= 12, **font_settings)

    def add_edges(self, sequences_list, sequences_filtered, threshold_distance):
        num_sequences = sequences_list.shape[0]
        for seq_1 in range(num_sequences):
            for seq_2 in range(seq_1 + 1, num_sequences):
                p = sequences_list[seq_1]
                q = sequences_list[seq_2]
                distance = np.sum(abs(sequences_list[seq_1] - sequences_list[seq_2]))
                print(distance)
                if distance < threshold_distance/num_sequences:
                    string1 = sequences_filtered[seq_1]
                    string2 = sequences_filtered[seq_2]
                    self.G.add_edge(string1, string2, weight=distance)
    @staticmethod
    def add_nodesize(G, selected_rows):
        nodesize = []
        label_numbers = {}
        label_sequences = {}
        for index, g in enumerate(G):
            
            for i in selected_rows["aaSeqCDR3"]:
                
                if i == g:
                    node = math.sqrt(selected_rows.loc[selected_rows["aaSeqCDR3"] == i, "readFraction"].values[0])
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
        return nodesize
    
    @staticmethod     
    def generate_plot(G, nodesize, ax, cmap):
        partition = community.best_partition(G)
        pos = nx.spring_layout(G)
        nx.draw_networkx_nodes(G, pos,node_size=nodesize, cmap=cmap, node_color=list(partition.values()), ax=ax)
        nx.draw_networkx_edges(G, pos, alpha=0.3, ax=ax)




