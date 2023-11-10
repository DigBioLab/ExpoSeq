import torch
import re
import pandas as pd
# + pip install sentencepiece
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
from sklearn.manifold import TSNE
import numpy as np 

class TransformerBased:

    def __init__(self, choice='Rostlab/prot_t5_xl_half_uniref50-enc'):
        
        self.choice = choice
        self.prep_model()

    def prep_model(self):
        model_types = ["Rostlab/ProstT5_fp16", "Rostlab/prot_t5_xl_uniref50", "Rostlab/prot_t5_base_mt_uniref50", "Rostlab/prot_bert_bfd_membrane", "Rostlab/prot_t5_xxl_uniref50", "Rostlab/ProstT5", "Rostlab/prot_t5_xl_half_uniref50-enc", "Rostlab/prot_bert_bfd_ss3", "Rostlab/prot_bert_bfd_localization", "Rostlab/prot_electra_generator_bfd", "Rostlab/prot_t5_xl_bfd", "Rostlab/prot_bert", "Rostlab/prot_xlnet", "Rostlab/prot_bert_bfd", "Rostlab/prot_t5_xxl_bfd"]
        assert self.choice in model_types, f"Your choice has to be a model from https://huggingface.co/Rostlab and has to be one of {model_types}"
        if "t5" in self.choice.lower():
            from transformers import T5Tokenizer, T5EncoderModel
            self.tokenizer = T5Tokenizer.from_pretrained(self.choice,
                                                    do_lower_case=False)

            self.model = T5EncoderModel.from_pretrained(self.choice)
        else:
            from transformers import BertModel, BertTokenizer
            self.tokenizer = BertTokenizer.from_pretrained(self.choice, do_lower_case=False )
            self.model = BertModel.from_pretrained(self.choice)
        

    
    def filter_sequences(self, sequencing_report, batch_size, experiments,binding_data, region_of_interest = "aaSeqCDR3"):
        report_batch = sequencing_report.groupby("Experiment").head(batch_size)
        selected_rows = report_batch.loc[report_batch["Experiment"].isin(experiments)]
        if binding_data is not None:
            #mix = selected_rows.merge(binding_data, on = "aaSeqCDR3", how = "outer")
            mix = pd.concat([selected_rows, binding_data])
            selected_rows = mix.fillna(0)
        selected_rows["cloneFraction"] = selected_rows["cloneFraction"].replace(0.0, max(selected_rows["cloneFraction"])) # all the added sequences from binding data get highest clone fraction to visualize them, otherwise they will not appear since the plot is fraction sensitive
        selected_rows = selected_rows.sort_values(by='cloneFraction', ascending=False)
        selected_rows.drop_duplicates(subset='aaSeqCDR3', keep="first", inplace=True) # remove duplicates due to adding of binding data
        sequences_filtered = selected_rows[region_of_interest]
        sequences = [" ".join(list(re.sub(r"[UZOB*]", "X", sequence))) for sequence in sequences_filtered]
        return sequences,sequences_filtered, selected_rows
        
    def prepare_sequences(self,sequences, device = "cpu"):
        ids = self.tokenizer.batch_encode_plus(sequences, add_special_tokens=True, padding="longest")
        input_ids = torch.tensor(ids['input_ids'])
        attention_mask = torch.tensor(ids['attention_mask'])
        return attention_mask, input_ids
    
    def get_result(self, sequences):
        if self.choice == "T5":
            attention_mask, input_ids = self.prepare_sequences(sequences, self.tokenizer)
            with torch.no_grad():
                embedding_repr = self.model(input_ids=input_ids, attention_mask=attention_mask)
        else:
            
            encoded_input = self.tokenizer(sequences, return_tensors='pt', padding=True, truncation=True)
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
    
