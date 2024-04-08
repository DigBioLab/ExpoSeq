from .tidy_protbert_embedding import TransformerBased
from bertviz.neuron_view import show
from bertviz import model_view
import re 
import seaborn as sns
import matplotlib.pyplot as plt
from captum.attr import visualization as viz
from captum.attr import IntegratedGradients, LayerConductance, LayerIntegratedGradients, LayerActivation
from captum.attr import configure_interpretable_embedding_layer, remove_interpretable_embedding_layer

import torch



class EmbeddingPrep:
    def __init__(self, model ,tokenizer, device) -> None:
        self.model = model
        self.tokenizer = tokenizer
        self.device = device    
    @staticmethod
    def create_interpretable_embedding(model):
        interpretable_embedding = configure_interpretable_embedding_layer(model, 'bert.embeddings.word_embeddings')
        return interpretable_embedding
    
    @staticmethod
    def construct_input_ref_pair(tokenizer, input, device = "cpu"):
        input = list("CAKDIGGGTRYYYYGMDVW")

        encoded_input = tokenizer.encode(input, add_special_tokens=False)
        ref_token_id = tokenizer.pad_token_id # A token used for generating token reference
        sep_token_id = tokenizer.sep_token_id # A token used as a separator between question and text and it is also added to the end of the text.
        cls_token_id = tokenizer.cls_token_id
        # construct input token ids
        input_ids = [cls_token_id] + encoded_input + [sep_token_id] 
        # construct reference token ids 
        ref_input_ids = [cls_token_id] + [ref_token_id] * len(encoded_input) + [sep_token_id]
        input_ids = torch.tensor([input_ids], device=device)
        ref_input_ids = torch.tensor([ref_input_ids], device=device)
        return input_ids, ref_input_ids
    
    def construct_whole_bert_embeddings(self, input_ids, ref_input_ids):

        output = self.model(input_ids)
        input_embeddings = output[0]

        ref_output = self.model(ref_input_ids)
        ref_input_embeddings = ref_output[0]
        
        return input_embeddings, ref_input_embeddings
    
    
    @staticmethod
    def summarize_attributions(attributions):
        if torch.__version__ >= '1.7.0':
            norm_fn = torch.linalg.norm
        else:
            norm_fn = torch.norm
        attributions = attributions.sum(dim=-1).squeeze(0)
        attributions = attributions / norm_fn(attributions)
        return attributions
    
    def peptide_forward_func(self, inputs, token_type_ids=None, position_ids=None, attention_mask=None):
        pred = self.model(inputs_embeds=inputs, token_type_ids=token_type_ids,
                    position_ids=position_ids, attention_mask=attention_mask)[0]

        print(pred.shape)
        print(inputs.shape)
        return pred
    
    def construct_input_ref_token_type_pair(self, input_ids, sep_ind=0):
        seq_len = input_ids.size(1)
        token_type_ids = torch.tensor([[0 if i <= sep_ind else 1 for i in range(seq_len)]], device=self.device)
        ref_token_type_ids = torch.zeros_like(token_type_ids, device=self.device)# * -1
        return token_type_ids, ref_token_type_ids
    
    def construct_input_ref_pos_id_pair(self, input_ids):
        seq_length = input_ids.size(1)
        position_ids = torch.arange(seq_length, dtype=torch.long, device=self.device)
        # we could potentially also use random permutation with `torch.randperm(seq_length, device=device)`
        ref_position_ids = torch.zeros(seq_length, dtype=torch.long, device=self.device)

        position_ids = position_ids.unsqueeze(0).expand_as(input_ids)
        ref_position_ids = ref_position_ids.unsqueeze(0).expand_as(input_ids)
        return position_ids, ref_position_ids
    
    @staticmethod
    def construct_attention_mask(input_ids):
        return torch.ones_like(input_ids)
    
    def layer_over_token(self, input:str) :
        layer_attrs_start = []
        layer_attrs_end = []
        layer_attn_mat_start = []
        layer_attn_mat_end = []
        
        input_ids, ref_input_ids = self.construct_input_ref_pair(self.tokenizer, input)
        input_embeddings, ref_input_embeddings = self.construct_whole_bert_embeddings(input_ids, ref_input_ids)
        token_type_ids, _= self.construct_input_ref_token_type_pair(input_ids)
        position_ids, _ = self.construct_input_ref_pos_id_pair(input_ids)
        attention_mask = self.construct_attention_mask(input_ids)


        for i in range(self.model.config.num_hidden_layers):
            print(input_embeddings.shape)
            lc = LayerConductance(self.peptide_forward_func, self.model.encoder.layer[i])

            layer_attributions_start = lc.attribute(inputs=input_embeddings,
                                                    baselines=ref_input_embeddings, 
                                                    target = 1,
                                                    )
           # layer_attributions_end = lc.attribute(inputs=input_embeddings,
            #                                      baselines=ref_input_embeddings,
             #                                     additional_forward_args=(token_type_ids,
              #                                                             position_ids, 
               #                                                            attention_mask, 1))
            
            layer_attrs_start.append(self.summarize_attributions(layer_attributions_start[0]))
      #      layer_attrs_end.append(self.summarize_attributions(layer_attributions_end[0]))

            layer_attn_mat_start.append(layer_attributions_start[1])
         #   layer_attn_mat_end.append(layer_attributions_end[1])# layer x seq_len
        layer_attrs_start = torch.stack(layer_attrs_start)
        # layer x seq_len
       # layer_attrs_end = torch.stack(layer_attrs_end)

        # layer x batch x head x seq_len x seq_len
        layer_attn_mat_start = torch.stack(layer_attn_mat_start)
        # layer x batch x head x seq_len x seq_len
        layer_attn_mat_end = torch.stack(layer_attn_mat_end)

        return layer_attrs_start, layer_attrs_end, input_ids

class PrepareData:
    
    def __init__(self, sequence, model = "Rostlab/prot_bert"):
        ModelManager = TransformerBased(model)
        self.sequence = self.prepare_sequence(sequence)
        self.model = ModelManager.model
        self.tokenizer = ModelManager.tokenizer
        
    @staticmethod
    def prepare_sequence(sequence):
        seq = [sequence]
        assert type(seq) == list, "Input must be a list of strings"
        sequences = [" ".join(list(re.sub(r"[UZOB*_]", "X", sequence))) for sequence in seq]
        return sequences
        
    def token_to_ids(self, encoded_input, ):
        input_ids = EmbeddingPrep.construct_input_ref_pair(self.tokenizer, encoded_input)
        indices = input_ids[0].detach().tolist()
        all_tokens = self.tokenizer.convert_ids_to_tokens(indices)
        return all_tokens
    
    def tokenize_input(self, ):
        encoded_input = self.tokenizer(self.sequence, return_tensors = "pt", truncation = True, max_length = 1000)
        embedding_repr = self.model(**encoded_input)
        return embedding_repr

    def prep_layer_over_token(self):
        PrepEmbed = EmbeddingPrep(self.model, self.tokenizer, device = "cpu")
        encoded_input = self.tokenize_input()
        layer_attrs_start, layer_attrs_end = PrepEmbed.layer_over_token(self.sequence)
        all_tokens = self.token_to_ids(encoded_input)
        return layer_attrs_start, layer_attrs_end, all_tokens
    
    
    def cleaning(self,  layer = 8):
        model_type = "bert"
        embedding_repr = self.tokenize_input()
        seq = [self.sequence]
        xlabels = list(seq[0])
        attention = embedding_repr[-1]  # Retrieve attention from model outputs
        layer = layer -1
        # Extract attention matrices for each layer and head
        attention_matrices = attention[layer].detach().numpy()
        return attention_matrices
    
PrepData = PrepareData(sequence = "CAKDIGGGTRYYYYGMDVW")      
layer_attrs_start, layer_attrs_end, all_tokens = PrepData.prep_layer_over_token()  
fig, ax = plt.subplots(figsize=(15,5))
xticklabels=all_tokens
yticklabels=list(range(1,13))
ax = sns.heatmap(layer_attrs_start.cpu().detach().numpy(), xticklabels=xticklabels, yticklabels=yticklabels, linewidth=0.2)
plt.xlabel('Tokens')
plt.ylabel('Layers')
plt.show()
    
#attention_matrices[0, 14, 1:-1, 1:-1] # seq no., attention heads = 16, seqlen, seqlen

#sns.heatmap(, xticklabels=xlabels, yticklabels=xlabels)
#plt.show()
#print(attention_matrices.shape)