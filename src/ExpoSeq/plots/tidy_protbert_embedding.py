import torch
import re
import pandas as pd
# + pip install sentencepiece
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
import numpy as np 
import umap 

class TransformerBased:

    def __init__(self, choice='Rostlab/prot_t5_xl_half_uniref50-enc', truncation = False, max_length = None):
        self.truncation = truncation
        self.max_length = max_length
        self.choice = choice
        self.prep_model()

    def prep_model(self):
        model_types = ["facebook/esm2_t6_8M_UR50D", "Rostlab/ProstT5_fp16", "Rostlab/prot_t5_xl_uniref50", "Rostlab/prot_t5_base_mt_uniref50", "Rostlab/prot_bert_bfd_membrane", "Rostlab/prot_t5_xxl_uniref50",
                       "Rostlab/ProstT5", "Rostlab/prot_t5_xl_half_uniref50-enc", "Rostlab/prot_bert_bfd_ss3", "Rostlab/prot_bert_bfd_localization",
                       "Rostlab/prot_t5_xl_bfd", "Rostlab/prot_bert", "Rostlab/prot_xlnet", "Rostlab/prot_bert_bfd", "Rostlab/prot_t5_xxl_bfd"]
        
        t5 = ["Rostlab/ProstT5", "Rostlab/ProstT5_fp16", "Rostlab/prot_t5_xl_uniref50", "Rostlab/prot_t5_xxl_uniref50", "Rostlab/prot_t5_xl_half_uniref50-enc", "Rostlab/prot_t5_xl_bfd", "Rostlab/prot_t5_xxl_bfd"]
        bertis = ["Rostlab/prot_bert_bfd_ss3", "Rostlab/prot_bert_bfd_localization", "Rostlab/prot_bert", "Rostlab/prot_xlnet", "Rostlab/prot_bert_bfd"]
        esm = ["facebook/esm2_t6_8M_UR50D"]
        assert self.choice in model_types, f"Your choice has to be a model from https://huggingface.co/Rostlab and has to be one of {model_types}"
        if self.choice in t5:
            from transformers import T5Tokenizer, T5EncoderModel
            self.tokenizer = T5Tokenizer.from_pretrained(self.choice,
                                                    do_lower_case=False,)

            self.model = T5EncoderModel.from_pretrained(self.choice)
        elif self.choice in bertis:
            from transformers import BertModel, BertTokenizer
            self.tokenizer = BertTokenizer.from_pretrained(self.choice, do_lower_case=False )
            self.model = BertModel.from_pretrained(self.choice)
        elif self.choice in esm:
            from transformers import AutoTokenizer, EsmModel
            self.tokenizer = AutoTokenizer.from_pretrained(self.choice)
            self.model = EsmModel.from_pretrained(self.choice)
        
        
        

    
    def filter_sequences(self, sequencing_report, batch_size, experiments,binding_data,
                         region_of_interest = "aaSeqCDR3", cf_column_name = "cloneFraction", sample_column_name = "Experiment"):
        report_batch = sequencing_report.groupby(sample_column_name).head(batch_size)
        selected_rows = report_batch.loc[report_batch[sample_column_name].isin(experiments)]
        selected_rows = selected_rows.drop_duplicates(subset = [region_of_interest])
        if binding_data is not None:
            mix = selected_rows.merge(binding_data, on = region_of_interest, how = "outer")
            selected_rows = mix.fillna(0)
        max_fraction = max(selected_rows[cf_column_name])
        selected_rows.loc[selected_rows[cf_column_name] == 0.0, cf_column_name] = max_fraction
        selected_rows = selected_rows.sort_values(by=cf_column_name, ascending=False)
        sequences_filtered = selected_rows[region_of_interest]
        sequences = [" ".join(list(re.sub(r"[UZOB*]", "X", sequence))) for sequence in sequences_filtered]
        return sequences,sequences_filtered, selected_rows
        
    def prepare_sequences(self,sequences, device = "cpu"):
        ids = self.tokenizer.batch_encode_plus(sequences, add_special_tokens=True, truncation = True, max_length = 1430)
        input_ids = torch.tensor(ids['input_ids'])
        attention_mask = torch.tensor(ids['attention_mask'])
        return attention_mask, input_ids
    
    def get_result(self, sequences):        
        if not self.truncation:    
            encoded_input = self.tokenizer(sequences, return_tensors='pt', padding=True)
        else:
            if self.max_length == None:
                max_length = 1440
            else:
                max_length = self.max_length
            encoded_input = self.tokenizer(sequences, return_tensors = "pt", truncation = True, max_length = max_length)
        embedding_repr = self.model(**encoded_input)
        return embedding_repr
    
    
    def embedding_per_seq(self, sequences, normalize = False):
        sequences_list = []
        for seq in sequences:
            embeddings = self.get_result(seq)
            last_hidden_state = embeddings.last_hidden_state
            maximum_length = last_hidden_state.shape[1]
           
            avg_seq = np.squeeze(last_hidden_state, axis=0)
            
            avg_seq = last_hidden_state.mean(dim = 1) # take average for each feature from all amino acids
            if normalize == True:
                #avg_seq = (avg_seq - torch.min(avg_seq)) / (torch.max(avg_seq) - torch.min(avg_seq))
                min_value = torch.min(avg_seq)
                shifted_distribution = avg_seq - min_value
                sum_shifted = torch.sum(shifted_distribution)
                normalized_distribution = shifted_distribution / sum_shifted
                avg_seq = normalized_distribution
            sequences_list.append(avg_seq.cpu().detach().numpy()[0])
        
        sequences_list = np.array(sequences_list)
        return sequences_list
    
    def do_pca(self, sequences, batch_size, pca_components):
        """
        output: X: principal components of the embedding which has the shape: x = batch_size, y = principal component
        """
      #  assert batch_size > pca_components, "Your batch size (Number of sequences) needs to be bigger than the pca components. Further you should decrease the pca components more, otherwise the dimension reduction does not make really sense."
        # batch size, sequence length, features

        # loop over each sequence and make two dimensional and take average for all amino acids: output size: (batch_size, 1024)
        sequences_list = self.embedding_per_seq(sequences)
         # shape y = 1024, x = batch_size (No. sequences)
        pca = PCA(n_components=pca_components)
        pca.fit(sequences_list)
        X = pca.transform(sequences_list)
        
        print("Explained variance after reducing to " + str(pca_components) + " dimensions:" + str(np.sum(pca.explained_variance_ratio_).tolist()))

        return X
    
    @staticmethod
    def do_tsne(X, perplexity, iterations_tsne):
        """
        input: vector with shape: x = n_samples and y: n_features
        """
        tsne = TSNE(n_components=2,
            verbose=0,
            perplexity=perplexity,
            n_iter=iterations_tsne)
        tsne_results = tsne.fit_transform(X)
        tsne_results = pd.DataFrame(tsne_results,
                                    columns = [["tsne1", "tsne2"]])
        return tsne_results
    @staticmethod
    def do_umap(X, n_neighbors = 15,min_dist = 0.2, random_seed = 42, densmap = True, n_components = 2, y = None):
        """_summary_

        Args:
            X (array): an array with the PC's from PCA
            n_neighbors (int): Larger values will result in more global structure being preserved at the loss of detailed local structure. In general this parameter should be between 5 to 50.
            random_seed: Set a certrain seed for reproducing your results
            densmap: This parameter allows you to visualize points more densily which are also more dense in all dimensions to each other. You can have an idea about this here: https://umap-learn.readthedocs.io/en/latest/densmap_demo.html

        Returns:
            pd.DataFrame: return a pandas dataframe with the two umap_dimensions
        """
        reducer = umap.UMAP(random_state = random_seed, n_neighbors = n_neighbors, densmap = densmap, n_components = n_components, min_dist = min_dist)
        if y == None:
            reduced_dim = reducer.fit_transform(X)
        else:
            reduced_dim = reducer.fit_transform(X, y = y)
        assert reduced_dim.shape[1] == 2
        results = pd.DataFrame(reduced_dim, columns = [["UMAP_1", "UMAP_2"]])
        return results
