import networkx as nx
import math
from textwrap import wrap
import numpy as np
from .tidy_protbert_embedding import TransformerBased
import matplotlib.pyplot as plt
import community.community_louvain as community
from sklearn.cluster import KMeans

class PrepareData:
    @staticmethod
    def input_check(model, list_experiments, batch_size):
        model_types_all = ["Rostlab/ProstT5_fp16", "Rostlab/prot_t5_xl_uniref50", "Rostlab/prot_t5_base_mt_uniref50", "Rostlab/prot_bert_bfd_membrane", "Rostlab/prot_t5_xxl_uniref50", "Rostlab/ProstT5", "Rostlab/prot_t5_xl_half_uniref50-enc", "Rostlab/prot_bert_bfd_ss3", "Rostlab/prot_bert_bfd_localization", "Rostlab/prot_electra_generator_bfd", "Rostlab/prot_t5_xl_bfd", "Rostlab/prot_bert", "Rostlab/prot_xlnet", "Rostlab/prot_bert_bfd", "Rostlab/prot_t5_xxl_bfd"]
        assert model in model_types_all, f"Please enter a valid model name which are\n{model_types_all}. You can find nearly all of the models at: https://huggingface.co/Rostlab"
        assert type(list_experiments) == list, "You have to give a list with the samples you want to analyze"
        assert type(batch_size) == int, "You have to give an integer as input for the batch_size"

        
    def cleaning(self, sequencing_report, model_choice, batch_size, list_experiments, region_of_interest, ):
        self.input_check(model_choice, list_experiments, batch_size)
        Transformer = TransformerBased(choice = model_choice)
        sequences,sequences_filtered,  selected_rows = Transformer.filter_sequences(sequencing_report, 
                                                                                    batch_size=batch_size,
                                                                                    experiments = list_experiments,
                                                                                    region_of_interest=region_of_interest,binding_data = None )
        sequences_filtered = sequences_filtered.to_list()
        sequences_list = Transformer.embedding_per_seq(sequences,
                                                       normalize = True)
        return sequences_list, sequences_filtered, selected_rows
    

class Network_Embedding:
    def __init__(self, ax, sequencing_report, list_experiments,model_choice, batch_size,font_settings,cmap = plt.cm.RdYlBu, nodesize = None, threshold_distance = 1, region_of_interest="aaSeqCDR3") -> None:
        self.ax = ax
        sequences_list, sequences_filtered, selected_rows = PrepareData().cleaning(sequencing_report, model_choice, batch_size, list_experiments, region_of_interest)
        self.G = nx.Graph()
        self.add_edges(sequences_list, sequences_filtered, threshold_distance)
        if nodesize == None:
            nodesize = self.add_nodesize(selected_rows, region_of_interest)
        else:
            nodesize = 500
        self.generate_plot(
                           nodesize = nodesize,
                           cmap = cmap)
        title = "\n".join(wrap("Network embedding for given samples", 40))
        self.ax.set_title(title, pad= 12, **font_settings)

    @staticmethod
    def k_means_clustering(distance_matrix):
        """
        This function follows the elbow methods to determine the number of clusters
        Args:
            distance_matrix (_type_): _description_
        """
        sse = [] # sum of squared errors
        sequence_size = distance_matrix.shape[0]
        max_cluster_number = sequence_size//3
        list_k = list(range(1, max_cluster_number))

        for k in list_k:
            km = KMeans(n_clusters=k)
            km.fit(distance_matrix)
            sse.append(km.inertia_)
        differences = np.diff(sse)
        
    def add_edges(self, sequences_list, sequences_filtered, threshold_distance):
        num_sequences = sequences_list.shape[0]
        distance_matrix = np.zeros((num_sequences, num_sequences))
        for seq_1 in range(num_sequences):
            for seq_2 in range(seq_1 + 1, num_sequences):
                p = sequences_list[seq_1]
                q = sequences_list[seq_2]
                distance = np.sqrt(np.sum(abs(sequences_list[seq_1] - sequences_list[seq_2])**2)) # euclidian distance of embedding coordinates for the 1024 dimensions
                if distance < threshold_distance/num_sequences:
                    string1 = sequences_filtered[seq_1]
                    string2 = sequences_filtered[seq_2]
                    self.G.add_edge(string1, string2, weight=distance)
                    
    def add_nodesize(self, selected_rows, region_of_interest):
        nodesize = []
        label_numbers = {}
        label_sequences = {}
        for index, g in enumerate(self.G):
            
            for i in selected_rows[region_of_interest]:
                
                if i == g:
                    node = math.sqrt(selected_rows.loc[selected_rows[region_of_interest] == i, "cloneFraction"].values[0])
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
    
    def generate_plot(self, nodesize, cmap):
        partition = community.best_partition(self.G)
        pos = nx.spring_layout(self.G)
        nx.draw_networkx_nodes(self.G, pos, cmap=cmap, node_color=list(partition.values()))
        nx.draw_networkx_edges(self.G, pos, alpha=1)



