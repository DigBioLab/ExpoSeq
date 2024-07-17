from tidy_protbert_embedding import TransformerBased
import re
import torch
import matplotlib.pyplot as plt
from contents.contact_prediction import ContactPredictionHead

Transformer = TransformerBased(choice = r"C:\Users\nilsh\my_projects\ExpoSeq\models\nanobody_model", custom_model = True)

tokenizer = Transformer.tokenizer

model = Transformer.model

sequence = [""]
sequences = [" ".join(list(re.sub(r"[UZOB*_]", "X", sequence))) for sequence in sequence]
encoded_input = tokenizer(sequences, return_tensors='pt', padding=True)
model.eval()
tokens = encoded_input['input_ids']
with torch.no_grad():
    output = model(**encoded_input, output_attentions=True)

# Extract attentions from the model output
attentions = torch.stack(output['attentions'], dim=1)  # B x L x H x T x T
batch_size, layers, heads, seqlen, _ = attentions.shape
# Create a padding mask based on tokens equaling the padding index
attentions = attentions.view(batch_size, layers * heads, seqlen, seqlen)
padding_mask = tokens.eq(model.embeddings.padding_idx)

from contents.contact_prediction import symmetrize, apc
attentions = attentions.to(
    model.contact_head.regression.weight.device
)  # attentions always float32, may need to convert to float16
attentions = apc(symmetrize(attentions))
attentions = attentions.permute(0, 2, 3, 1)

contact_map = model.contact_head.activation(model.contact_head.regression(attentions).squeeze(3))
# Apply the padding mask to the attentions if necessary
#if padding_mask is not None:
 #   attention_mask = 1 - padding_mask.type_as(attentions)
 #   attention_mask = attention_mask.unsqueeze(1) * attention_mask.unsqueeze(2)
 #   attentions *= attention_mask[:, None, None, :, :]

# Assuming your contact_head expects tokens and attentions, and computes contacts
#contacts = model.contact_head(tokens, attentions)
#print(contacts.shape)

import matplotlib.pyplot as plt
import numpy as np

# Assuming `output_tensor` is your tensor with shape [1, 102, 102]
# Convert the tensor to a numpy array and remove the first dimension since it's single
heatmap_data = contact_map.detach().squeeze().numpy()
import seaborn as sns
sns.heatmap(heatmap_data)
# Create the plot
