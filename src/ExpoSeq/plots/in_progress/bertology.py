from tidy_protbert_embedding import TransformerBased
import re
import torch
import matplotlib.pyplot as plt
from contents.contact_prediction import ContactPredictionHead





class SetupModel:
    def __init__(self, model_path:str):
        self.Transformer = TransformerBased(choice = model_path, custom_model = True)
        self.tokenizer = self.Transformer.tokenizer
        self.model = self.Transformer.model
    
    def get_attention(self, sequence:list):
        sequences = [" ".join(list(re.sub(r"[UZOB*_]", "X", sequence))) for sequence in sequence]
        encoded_input = self.tokenizer(sequences, return_tensors='pt', padding=True)
        self.model.eval()
        tokens = encoded_input['input_ids']
        with torch.no_grad():
            output = self.model(**encoded_input, output_attentions=True)
        attentions = torch.stack(output['attentions'], dim=1)
                
        # Remove [CLS] and [SEP] token attentions
        # Assuming [CLS] at index 0 and [SEP] at the last index of each sequence
        seq_len = attentions.shape[3]  # Assuming all sequences are padded to the same length
        # Remove the first and last token (typically [CLS] and [SEP])
        attentions = attentions[:, :, :, 1:seq_len-1, 1:seq_len-1]
        return attentions


#Setup = SetupModel(model_path = r"C:\Users\nilsh\my_projects\ExpoSeq\models\nanobody_model")
#attentions = Setup.get_attention(sequence = ["GDIAGLNNMGWYRQAPGKQRELVAVQARGGNTNYTDSVKGRFTISRNNAGNTVYLQMNNLKSEDTAVYYCYATVGNWYTSGYYVDDYWGQGTQVTVSS_"])
#attention_shape = attentions.shape

import numpy as np
    

    
from Bio import AlignIO
from Bio import SeqIO
from Bio.Align.Applications import MuscleCommandline
import editdistance

class MSACluster:
    @staticmethod
    def make_seq_record_of_fasta(fasta_file:str) -> list:
        seq_records = []

        # Open the FASTA file and read each record
        with open(fasta_file, 'r') as fasta_file:
            for seq_record in SeqIO.parse(fasta_file, 'fasta'):
                seq_records.append(seq_record)
        return seq_records
    
    def run_msa(self, fasta_path, out = "aligned_sequences.fasta") -> list[str]:
        cline = MuscleCommandline(r"C:\Users\nilsh\Downloads\muscle3.8.31_i86win32.exe", input = fasta_path, out = out)
        cline()
        alignment = AlignIO.read(out, 'fasta')
        sequence_strings = []
        for record in alignment:
            # Convert the Seq object to a string and append to the list
            sequence_strings.append(str(record.seq))
        return sequence_strings
    
    @staticmethod
    def calculate_distance(seq1, seq2):
        """Calculate the edit distance between two sequences, including gap penalties."""
        return editdistance.distance(seq1, seq2)
    
    def distance_on_gaps(self, sequence_strings:list, max_distance = 3) -> list[list]:
        clusters = []
        # Iterate over each sequence
        for seq in sequence_strings:
            # Try to find a cluster for the sequence
            found_cluster = False
            for cluster in clusters:
                # Check if the sequence is close enough to any element in the cluster
                if any(editdistance.distance(seq, member) <= max_distance for member in cluster):
                    cluster.append(seq)
                    found_cluster = True
                    break
            # If no suitable cluster is found, start a new cluster
            if not found_cluster:
                clusters.append([seq])
        return clusters
        
    @staticmethod
    def find_variable_positions(clusters:list[list]) -> list[tuple]:
        """Finds the variable positions in the cluster of sequences. 

        Args:
            clusters (list[list]): Clusters found by the distance_on_gaps function. Number of clusters equals the number of items in the first list. 
            Clusters contain multiple sequence aligned sequences, so the sequences have all the same length.

        Returns:
            list[tuple]: Number of tuples equals the number of clusters. First item is the list of sequences and second item is the list of variable positions
        """
        variable_positions = []
        for cluster in clusters:
            if len(cluster) == 1:  # If there's only one sequence, no variability
                variable_positions.append((cluster, []))
                continue
            
            # Initialize a set to keep track of varying positions
            varying_indexes = set()
            
            # Compare each sequence with every other sequence in the cluster
            reference_sequence = cluster[0]
            sequence_length = len(reference_sequence)
            
            for i in range(sequence_length):
                reference_char = reference_sequence[i]
                for sequence in cluster[1:]:
                    if sequence[i] != reference_char:
                        varying_indexes.add(i)
                        break
            
            # Store the result as a tuple of the cluster and the sorted list of varying positions
            variable_positions.append((cluster, sorted(varying_indexes)))
        
        return variable_positions
    
    def find_fixed_positions(self, clusters:list[list]) -> list[tuple]:
        variable_positions_sequences: list[tuple] = self.find_variable_positions(clusters)
        
        for index, cluster in enumerate(variable_positions_sequences):
            example_sequence = cluster[0][0]
            variable_positions_cluster = cluster[1]
            all_positions = list(range(len(example_sequence)))
            fixed_positions = [pos for pos in all_positions if pos not in variable_positions_cluster]
            variable_positions_sequences[index] = (cluster[0], fixed_positions)
        return variable_positions_sequences
    
    @staticmethod
    def translate_positions(cluster):
        gapped_sequence = cluster[0][0]
        positions = cluster[1]
        ungapped_sequence = gapped_sequence.replace('-', '')  # Remove gaps
        mapping = []  # This will map gapped indices to ungapped indices
        ungapped_index = 0
        
        for index, char in enumerate(gapped_sequence):
            if char != '-':
                mapping.append(ungapped_index)
                ungapped_index += 1
            else:
                mapping.append(None)  # No corresponding index in ungapped sequence

        translated_positions = [mapping[pos] for pos in positions if mapping[pos] is not None]
        return ungapped_sequence, translated_positions
            

aligned_sequences = MSACluster().run_msa(r"c:\Users\nilsh\Downloads\aligned_sequences.fasta")
assert type(aligned_sequences) == list, "The aligned sequences should be a list of strings"
assert len(aligned_sequences) == 18, "The number of aligned sequences should be 18"
clusters = MSACluster().distance_on_gaps(aligned_sequences, max_distance = 3)
assert len(clusters) == 7, "The number of clusters should be 3"
variable_positions = MSACluster.find_variable_positions(clusters)
assert type(variable_positions[0]) == tuple
assert variable_positions[3][1] == [5, 7, 13, 14, 18, 24, 25], f"The variable positions are {variable_positions[3][1]}"
assert len(variable_positions[0][0]) == 1, "there is only one sequence in this cluster"
assert len(variable_positions[0][1]) == 0, "there is only one sequence in this cluster, so you cannot find variable positions"
fixed_positions = MSACluster().find_fixed_positions(clusters)
all_pos = list(range(len(fixed_positions[3][0]))) # all positions in the first sequence of the 4th cluster
assert not any(elem in fixed_positions[3][1] for elem in [5, 7, 13, 14, 18, 24, 25]), "The fixed positions should not be the variable positions"
# property function f is often symmetric: which means that i,j and j,i return 1 or 0, respectively. 
# The asymmetric case when i,j = 1 and j,i = 0 would happen if your attention is direction dependent. Direction should not be important for protein structures

# we will start with symmetric properties - so you parse only one list of residues
class Bertology:
    def __init__(self, residues:list, sequence:str, function = "binding_site", **kwargs):
        """_summary_

        Args:
            residues (list): This is a list of indexes of residues in the sequence.
        """
        self.residues = residues
        self.decision_func = getattr(self, f"_f_{function}")
        self.sequence = sequence
        
    def _f_binding_site(self, i, j):
        """Symmetric function. It returns 1 if i or j is in residues, otherwise 0. can be used if only the position of the residues is of interest (interesting for cdr3)

        Args:
            i (_type_): x value in specific attention head
            j (_type_): y value in specific attention head

        Returns:
            _type_: 
        """
        if i in self.residues or j in self.residues:
            return 1
        else:
            return 0


    # a weighted decision function could be interesting which depends on multiple parameters and which sum of weights is 1
    
    def _f_hydrophob(self, i, j, ):
        kyte_doolittle = {
            'R': -4.5, 'K': -3.9, 'D': -3.5, 'E': -3.5, 'N': -3.5, 'Q': -3.5,
            'H': -3.2, 'S': -0.8, 'T': -0.7, 'G': -0.4, 'A': 1.8, 'M': 1.9,
            'C': 2.5, 'Y': -1.3, 'W': -0.9, 'P': -1.6, 'V': 4.2, 'I': 4.5,
            'L': 3.8, 'F': 2.8, "_": -4.5, "*": -4.5
        }

        min_value = min(kyte_doolittle.values())
        max_value = max(kyte_doolittle.values())

        normalized_kyte_doolittle = {key: (value - min_value) / (max_value - min_value) for key, value in kyte_doolittle.items()}
        # doesnt work well, probably because hydrobhobicity is a very universal characteristic which may appear in each head
        aa_i = self.sequence[i]
        aa_j = self.sequence[j]
        hydrophobicity_vector_mean = normalized_kyte_doolittle[aa_i] * normalized_kyte_doolittle[aa_j]
        return hydrophobicity_vector_mean
    

    # should be reprogrammed with matrix multiplication, because it is very slow currently
    def compute_pa_f(self, attentions):
        assert len(attentions.shape) == 5, "The attentions tensor must have shape [batch_size, layers, heads, max_seq_len, max_seq_len]"
        assert attentions.shape[3] == len(self.sequence), "The sequence length of the attentions tensor must match the length of the sequence. You should remove the CLS and SEP token"
        numerator = torch.zeros((attentions.shape[0], attentions.shape[1], attentions.shape[2]), device=attentions.device) # shape [batch_size, layers, heads] holds the sum of each sequence, layer and head
        denominator = torch.zeros_like(numerator)
        for sequence in range(attentions.shape[0]):
            for layer in range(attentions.shape[1]):
                # Iterate through each head in the current layer
                for head in range(attentions.shape[2]):
                    # Iterate through the rows of the attention matrix for the current head
                    for i in range(0, attentions.shape[3]):
                        # Iterate through the columns of the attention matrix for the current token
                        for j in range(0,attentions.shape[4]):
                            # Check if the attention weight is greater than zero
                            alpha_ij = attentions[sequence, layer, head, i, j]
                            if alpha_ij > 0:
                                # Apply function f to the indices and multiply by the attention weight
                                numerator[sequence, layer, head] += self.decision_func(i, j) * alpha_ij
                                denominator[sequence, layer, head] += alpha_ij
        pa_f = numerator / denominator # should be the mean of the attention weights per head
        assert pa_f.shape[0] == attentions.shape[0], "The batch size of the pa_f tensor should match the attentions tensor"
        assert pa_f.shape[1] == attentions.shape[1], "The number of layers of the pa_f tensor should match the attentions tensor"
        assert pa_f.shape[2] == attentions.shape[2], "The number of heads of the pa_f tensor should match the attentions tensor"
        return pa_f
        
import seaborn as sns

class PlotAttention:
    def __init__(self, matrix, cmap = "blues", figure_no = 1):
        self.matrix = matrix
        self.create_fig(figure_no)
        self.cmap = cmap

    def create_fig(self, figure_no):
        self.fig = plt.figure(figure_no)
        self.ax = self.fig.gca()
        
    def plot_mean_head(self, sequence:str, residues:list):
        assert len(self.matrix.shape) == 3, "The input matrix should three dimensions. the first one is the batch size."
        layers_head_tensor = self.matrix[0, :, :]
        layers_head_numpy = layers_head_tensor.cpu().numpy()
        sns.heatmap(layers_head_numpy, ax = self.ax, cmap = self.cmap)
        self.ax.set_xlabel("Heads")
        self.ax.set_ylabel("Layers")
        self.ax.set_title(f"Mean Attention per Head for {sequence} and residues {residues}")


    def plot_residue_residue(self, sequence:str, attentions, no_heads_average = 5):
        """This function creates a heatmap of the mean attention weights for the top n heads given the residue settings of pa_f.

        Args:
            sequence (str): Sequence which was the input of the model
            attentions (tensor): Tensor with shape [batch_size, layers, heads, max_seq_len, max_seq_len]. This contains the attention weights of the model.
            no_heads_average (int, optional): Here you choose how many of the top attention heads for the given constraints you would like to choose to calculate the average from. Defaults to 5.
        """
        assert len(self.matrix.shape) == 3, "The input matrix should three dimensions. the first one is the batch size."
        assert type(sequence) == str, "The sequence should be a string"
        layers_head_tensor = self.matrix[0, :, :]
        layers_head_numpy = layers_head_tensor.cpu().numpy()
        sorted_indices = np.argsort(layers_head_numpy, axis=None)[::-1]
        top_indices_flat = sorted_indices[:no_heads_average] 
        top_indices = np.unravel_index(top_indices_flat, layers_head_numpy.shape)
        heads_top = attentions[0, top_indices[0], top_indices[1], :, :].cpu().numpy()
        assert len(heads_top.shape) == 3, "The shape of the heads_top should be 3"
        assert heads_top.shape[0] == no_heads_average, "The number of heads should be the same as the no_heads_average"
        heads_top_mean = np.mean(heads_top, axis = 0) # you average over first dimension which is the number of heads
        sns.heatmap(heads_top_mean, cmap = "Blues", ax = self.ax)
        self.ax.set_xticks(np.arange(len(sequence)) + 0.5)  # Centering the labels
        self.ax.set_xticklabels(list(sequence))  # Setting labels to letters
        self.ax.set_yticks(np.arange(len(sequence)) + 0.5)  # Centering the labels
        self.ax.set_yticklabels(list(sequence))  # Setting labels to letters
        self.ax.set_title(f"Mean Attention for top {no_heads_average} for {sequence}.")

specific_cluster = fixed_positions[3]

sequence, adjusted_positions = MSACluster.translate_positions(specific_cluster)
Setup = SetupModel(model_path = r"C:\Users\nilsh\my_projects\ExpoSeq\models\nanobody_model")
attentions = Setup.get_attention(sequence = [sequence])

Berto = Bertology(residues = adjusted_positions, sequence = sequence, function = "binding_site")
pa_f = Berto.compute_pa_f(attentions)
import matplotlib.pyplot as plt
#Plotter = PlotAttention(matrix = pa_f, figure_no = 1)
#Plotter.plot_mean_head(sequence = sequence, residues = adjusted_positions)
Plotter = PlotAttention(matrix = pa_f, figure_no = 2)
Plotter.plot_residue_residue(sequence = sequence, attentions = attentions, no_heads_average=5)