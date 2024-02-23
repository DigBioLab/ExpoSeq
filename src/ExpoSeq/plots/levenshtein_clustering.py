import matplotlib.pyplot as plt
import networkx as nx
import math
import pandas as pd
import community.community_louvain as community
import editdistance




class PrepareData:
    @staticmethod
    def check_input(samples, max_ld, min_ld, batch_size, binding_data, antigens):
        assert type(samples) == list, "You have to give a string as input for the sample"
        assert type(max_ld) == int, "You have to give an integer as input for the maximum levenshtein distance"
        assert type(min_ld) == int, "You have to give an integer as input for the minimum levenshtein distance"
        assert type(batch_size) == int
        assert batch_size > 1
        assert max_ld > min_ld
        assert max_ld >= 1
        if binding_data is not None:
            assert type(antigens) == list
            for antigen in antigens:
                assert antigen in binding_data.columns.tolist()
        
    
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
        
    def calc_edges(self, sequencing_report,samples, batch_size = 1000, max_ld = 2, min_ld = 0, region_string = "aaSeqCDR3",
                   binding_data = None, antigens = None, experiment_column = "Experiment"):
        """ Create the Graph object for networkx which can be used to visualize the clusters 

        Args:
            sequencing_report (pd.DataFrame):  report which contains all sequencing information from all samples after the processing with mixcr
            samples (list): List containing the names of the samples which should be used for the final plot. 
            batch_size (_type_): Equals to the number of sequences which are drawn from the sequencing report for the embedding and the final plot. Defaults to 1000.
            max_ld (int): Is the maximum allowed difference between two sequences based on levenshtein distance to be still in one cluster. Default is 2.
            min_ld (int): Is the minimum allowed levenshtein distance between two sequence to be still in one cluster. Default is 0. 
            region_string (str): A string which indicates the column name from which the amino acid sequences should be taken from.
            binding_data (pd.DataFrame, optional): Dataframe which contains the sequences and the binding values to the antigens. Defaults to None.
            antigens (list, optional): list containing the names of the antigens which are the headers of the binding data. From these antigens the binding data will be taken for the plot. Defaults to None.
            experiment_column (str, optional): Name of the column which contains the sample names in the sequencing report. Defaults to "Experiment".

        Returns:
            G: Graph object for networkx containing the nodes and edges. 
            report: tidied report as input for the threshold algorithm based on levenshtein distance
        """
        self.check_input(samples, max_ld, min_ld, batch_size, binding_data, antigens)
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
        to_remove = [n for (n, deg) in G_deg if deg == 0] # degree = 0 means that these nodes are not connected to any other nodes. Since they have their own cluster they can be removed
        G.remove_nodes_from(to_remove)
        sequences_clustered = list(G.nodes)
        report = report[report[region_string].isin(sequences_clustered)] # remove those  which are not in self.G
        return G, report



class LevenshteinClustering:
    def __init__(self, sequencing_report, samples, ax = None, region_string = "aaSeqCDR3", max_ld = 2, min_ld =0, 
                 batch_size = 200,label_type = "numbers", font_settings = {}, binding_data = None, antigens = None, prefered_cmap="viridis"):
        self.ax = ax
        self.G, report = PrepareData().calc_edges(sequencing_report, samples, batch_size,
                                                  max_ld, min_ld, region_string, binding_data, antigens)
        nodesize, label_numbers, label_sequences = self.calculate_nodesize(report, region_string, self.G)
        label_sequences = self.create_sequence_labels(self.G)
        partition = self.get_partition(self.G)
        if ax != None:
            partition = self.create_network(label_type, label_numbers, label_sequences,
                                            nodesize, partition, region_string, report, binding_data, antigens, prefered_cmap)
            if font_settings != {}:
                self.add_header(font_settings, samples)
            
        self.cluster_report = self.generate_report(partition)
    
    
    @staticmethod
    def calculate_nodesize(sample_report, region_string, G, cf_column_name = "cloneFraction"):
        """ Calculates the nodesize and sets the labels for the corresponding sequences. The label can be either a number which indicates the id of the sequence in the report (sample_report) or a sequence.

        Args:
            sample_report (pd.DataFrame): is the report as output from PrepareData.tidy()
            region_string (str): A string which indicates the column name from which the amino acid sequences should be taken from.
            G (networkx.Graph): Graph object from networkx which contains the nodes and the edges for the plot.
            cf_column_name (str, optional):  Name of the column which contains the clone fraction in the sequencing report. Defaults to "cloneFraction".

        Returns:
            nodesize: size of the nodes in order
            label_numbers: dictionary containing the ids for the corresponding sequences
            label_sequencse: dicitinoary containing the sequences as label for the corresponding sequences
        """
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
            if n == 9:
                label_sequences[g] = g
                n = 0
            else:
                label_sequences[g] = ""
            n += 1
        return label_sequences
    
    @staticmethod
    def get_partition(G):
        """This function is used to partition the graph G into communities using the Louvain method, which is a popular algorithm for community detection. It returns a dictionary where the keys are the nodes of the graph G, and the values are the community assignments for each node.

        Args:
            G (networkx.Graph): Graph object from networkx which contains the nodes and the edges for the plot.

        Returns:
            _type_: _description_
        """
        partition = community.best_partition(G)
        return partition
    
    @staticmethod
    def map_binding(G, binding_data, region_of_interest, report, antigens, prefered_cmap, partition):
        """_summary_

        Args:
            G (networkx.Graph): Graph object from networkx which contains the nodes and the edges for the plot.
            binding_data (pd.DataFrame, optional): Dataframe which contains the sequences and the binding values to the antigens. Defaults to None.
            region_of_interest (str): A string which indicates the column name from which the amino acid sequences should be taken from.
            report: tidied report from PrepareData().tidy()
            antigens (list, optional): list containing the names of the antigens which are the headers of the binding data. From these antigens the binding data will be taken for the plot. Defaults to None.
            prefered_cmap (str): cmap for the colorbar
            partition (dict): output from LevenshteinClustering.get_partition()
        """
        node_colors = {}
        assert region_of_interest in binding_data.columns, "You do not have binding data for the corresponding region of interest"
       # new_df = binding_data.loc[:, antigens]     
        new_df = report.loc[:, antigens]
        maximum_binding_value = new_df.max().max()
        minimum_binding_value = new_df.min().min()
        if maximum_binding_value == 0:
            print("no binding values mapped")
            node_colors_norm_heatmap = list(partition.values())
            sm = None
        else:
            cmap = plt.colormaps.get_cmap(prefered_cmap)
            norm = plt.Normalize(minimum_binding_value, maximum_binding_value)
            sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm) # for colorbar
            for i, g in enumerate(G):
                sequence_row = report[report[region_of_interest] == g]
                right_binding_value = max([sequence_row[column].iloc[0] for column in antigens]) # if you have multiple antigens this will find the right value to map
                sequence = sequence_row[region_of_interest].iloc[0]
                node_colors[sequence] = right_binding_value
            node_colors_norm_heatmap = [cmap(norm(node_colors[node])) if node_colors[node] != 0 else 'gray' for node in G.nodes()] # all nodes without binding values are gray
        return node_colors_norm_heatmap, sm , maximum_binding_value
        
        
    def create_network(self, label_type, label_numbers, label_sequences, nodesize, partition, region_of_interest, report, binding_data, antigens, prefered_cmap):
        if binding_data is not None:
            node_colors, sm, maximum_binding_value = self.map_binding(self.G, binding_data, region_of_interest, report, antigens, prefered_cmap, partition)
        else:
            maximum_binding_value = 0
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
        if maximum_binding_value > 0:
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
        
        

#import pandas as pd
#file = r"C:\Users\nilsh\OneDrive\Desktop\DTU\NGS_pipeline\data\Binding_data\Chris_main_df.csv"
#binding_data = pd.read_csv(file)
#antigens = ["Ecarpholin_S_DDEL_5ug/mL", "Myotoxin_II_DDEL_5ug/mL"]
#filename = r"C:\Users\nilsh\my_projects\ExpoSeq\my_experiments\test_all_samples\sequencing_report.csv"
#sequencing_report = pd.read_csv(filename)
#sequencing_report["cloneFraction"] = sequencing_report["readFraction"]
#fig = plt.figure(1)
#ax = fig.gca()
#samples = ["Library_1_F11_2", "test_1_F10"]
#font_settings = {'fontfamily': 'serif', 'fontsize': '18', 'fontstyle': 'normal', 'fontweight': 'bold'}
#LevenshteinClustering(sequencing_report, samples, ax, font_settings=font_settings, batch_size=800, binding_data=binding_data, antigens = antigens)
#plt.show()






