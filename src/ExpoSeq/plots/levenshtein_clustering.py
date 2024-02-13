import matplotlib.pyplot as plt
import networkx as nx
import math
import pandas as pd
import community
from .layout_finder import best_layout
from textwrap import wrap
import editdistance




class PrepareData:
    @staticmethod
    def check_input(samples, max_ld, min_ld, batch_size, binding_data):
        assert type(samples) == list, "You have to give a string as input for the sample"
        assert type(max_ld) == int, "You have to give an integer as input for the maximum levenshtein distance"
        assert type(min_ld) == int, "You have to give an integer as input for the minimum levenshtein distance"
        assert type(batch_size) == int
        assert batch_size > 1
        assert max_ld > min_ld
        assert max_ld >= 1
        assert type(antigens) == list or type(antigens) == None
        assert antigens in binding_data.columns.tolist() or type(binding_data) == None # for case when class is used to solely cluster ngs sequences
        
    
    @staticmethod
    def clean_data(sequencing_report,samples, batch_size, experiment_column, binding_data, antigens, region_of_interest):
        sample_report = sequencing_report.loc[sequencing_report[experiment_column].isin(samples)]
        sample_report = sample_report.groupby(experiment_column).head(batch_size)
        sample_report = sample_report.drop_duplicates(subset = [region_of_interest])
        if binding_data is not None:
            merged_columns = [region_of_interest] + antigens
            binding_data = binding_data[merged_columns]
            mix = sample_report.merge(binding_data, on = region_of_interest, how = "outer")
            sample_report = mix.fillna(0)
        return sample_report
        
    def calc_edges(self, sequencing_report,samples, batch_size, max_ld, min_ld, region_string,
                   binding_data = None, antigens = None, experiment_column = "Experiment"):
        self.check_input()
        report = self.clean_data(sequencing_report, samples, batch_size, experiment_column, binding_data, antigens, region_string)
        G = nx.Graph()
        aa = list(report[region_string])
        num_strings = len(aa)
        string_index = [*range(0, num_strings, 1)]
        max_ld = max_ld
        min_ld = min_ld
        for seq in aa:
            G.add_node(seq)
        for i in range(len(aa)):
            string1 = aa[i]
            index_1 = i
            for index in range(index_1):
                string2 = aa[index]
                distance = editdistance.distance(string1, string2)
                if distance>= min_ld and distance <= max_ld: # sequences are not added to Graph if tehy dont meet the threshold requirements
                    G.add_edge(string1, string2)
        G_deg = G.degree()
        to_remove = [n for (n, deg) in G_deg if deg == 0]
        G.remove_nodes_from(to_remove)
        return G, report



class LevenshteinClustering:
    def __init__(self, sequencing_report, samples, ax = None, region_string = "aaSeqCDR3", max_ld = 2, min_ld =0, 
                 batch_size = 200,label_type = "numbers", font_settings = {}, binding_data = None, antigens = None):
        self.ax = ax
        self.G, report = PrepareData().calc_edges(sequencing_report, samples, batch_size,
                                                  max_ld, min_ld, region_string, binding_data, antigens)
        nodesize, label_numbers, label_sequences = self.calculate_nodesize(report, region_string, self.G)
        label_sequences = self.create_sequence_labels(self.G)
        partition = self.get_partition(self.G)
        if ax != None:
            partition = self.create_network(label_type, label_numbers, label_sequences,
                                            nodesize, partition, region_string, report, binding_data, antigens)
            if font_settings != {}:
                self.add_header(font_settings, samples)
            
        self.cluster_report = self.generate_report(partition)
    
    
    @staticmethod
    def calculate_nodesize(sample_report, region_string, G, cf_column_name = "cloneFraction"):
        nodesize = []
        label_numbers = {}
        label_sequences = {}
        for index, g in enumerate(G):
            for i in sample_report[region_string]:
                if i == g:
                    node = math.sqrt(sample_report.loc[sample_report[region_string] == i, cf_column_name].values[0])
                    if node == 0:
                        node = 0.5 # this is for the once which do not have a clone count, so the binding data. otherwise they are not visible 
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

        return nodesize, label_numbers, label_sequences

    @staticmethod
    def create_sequence_labels(G):
                ## labels
        label_sequences = {}
        node_ids = list(G.nodes())
        n = 0
        for index, g in enumerate(G):
            print(g)
            if n == 9:
                label_sequences[g] = g
                n = 0
            else:
                label_sequences[g] = ""
            n += 1
        return label_sequences
    
    @staticmethod
    def get_partition(G):
        partition = community.best_partition(G)
        return partition
    
    
    def map_binding(self, binding_data, region_of_interest, report, antigens):
        node_colors = {}
        assert region_of_interest in binding_data.columns, "You do not have binding data for the corresponding region of interest"
        new_df = binding_data.loc[:, antigens]     
        maximum_binding_value = new_df.max().max()
        minimum_binding_value = new_df.min().min()
        cmap = plt.colormaps.get_cmap("viridis")
        norm = plt.Normalize(minimum_binding_value, maximum_binding_value)
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm) # for colorbar
        for i, g in enumerate(self.G):
            sequence_row = report[report[region_of_interest] == g]
            right_binding_value = max([sequence_row[column].iloc[0] for column in antigens]) # if you have multiple antigens this will find the right value to map
            sequence = sequence_row[region_of_interest].iloc[0]
            node_colors[sequence] = right_binding_value
        node_colors_norm_heatmap = [cmap(norm(node_colors[node])) if node_colors[node] != 0 else 'gray' for node in self.G.nodes()] # all nodes without binding values are gray
        return node_colors_norm_heatmap, sm 
        
        
    def create_network(self, label_type, label_numbers, label_sequences, nodesize, partition, region_of_interest, report, binding_data, antigens):
        if binding_data is not None:
            node_colors, sm = self.map_binding(binding_data, region_of_interest, report, antigens)
        else:
            node_colors = list(partition.values())
        pos = nx.spring_layout(self.G)
        nx.draw_networkx_nodes(self.G, pos, node_size=nodesize, cmap=plt.cm.RdYlBu, node_color=node_colors, ax=self.ax)
        nx.draw_networkx_edges(self.G, pos, alpha=0.3, ax=self.ax)
        if label_type == "numbers":
            label = label_numbers 
        else:
            label = label_sequences 
        
        nx.draw_networkx_labels(self.G, pos,
                        labels=label,
                        alpha=0.5,
                        font_color="grey",
                        horizontalalignment="right",
                        verticalalignment="bottom", ax=self.ax)
        if binding_data is not None:
            plt.colorbar(sm, ax = self.ax)
        return partition
        
    def add_header(self, font_settings, samples):
        original_fontsize = font_settings["fontsize"]
        font_settings["fontsize"] = 22
        self.ax.set_title(f"Connected components of {', '.join(samples)}", **font_settings)
        font_settings["fontsize"] = original_fontsize
        self.ax.set_axis_off()
        
    @staticmethod
    def generate_report(partition):
        cluster_report = pd.DataFrame([partition]).T
        cluster_report['Sequences'] = cluster_report.index
        cluster_report.reset_index(drop=True, inplace=True)
        cluster_report.rename(columns={0: 'Cluster No.'}, inplace=True)
        return cluster_report
        
        

import pandas as pd
file = r"C:\Users\nilsh\OneDrive\Desktop\DTU\NGS_pipeline\data\Binding_data\Chris_main_df.csv"
binding_data = pd.read_csv(file)
antigens = ["Ecarpholin_S_DDEL_5ug/mL", "Myotoxin_II_DDEL_5ug/mL"]
filename = r"C:\Users\nilsh\my_projects\ExpoSeq\my_experiments\test_all_samples\sequencing_report.csv"
sequencing_report = pd.read_csv(filename)
sequencing_report["cloneFraction"] = sequencing_report["readFraction"]
fig = plt.figure(1)
ax = fig.gca()
samples = ["Library_1_F11_2", "test_1_F10"]
font_settings = {'fontfamily': 'serif', 'fontsize': '18', 'fontstyle': 'normal', 'fontweight': 'bold'}
LevenshteinClustering(sequencing_report, samples, ax, font_settings=font_settings, batch_size=800, binding_data=binding_data, antigens = antigens)
plt.show()






