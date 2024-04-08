from antiberty import AntiBERTyRunner
import seaborn
import numpy as np
import matplotlib.pyplot as plt
from iglabel import IMGT
antiberty = AntiBERTyRunner()

sequences = [
    "CATVRWERGVFDPW",
 #   "DVVMTQTPFSLPVSLGDQASISCRSSQSLVHSNGNTYLHWYLQKPGQSPKLLIYKVSNRFSGVPDRFSGSGSGTDFTLKISRVEAEDLGVYFCSQSTHVPYTFGGGTKLEIK",
]
embeddings = antiberty.embed(sequences)

print(embeddings[0].shape)

embeddings, attentions = antiberty.embed(sequences, return_attention=True)
list_attentions = []
for seq in attentions:
    seq = seq.numpy()
    list_attentions.append(seq)

example = list_attentions[0]
a, b= IMGT(sequences = [sequences[0]], regions = ["CDR3"])
imgt_labels = list(a.values())[0]

seaborn.heatmap(example[-1, 0, :, :], xticklabels=imgt_labels, yticklabels=imgt_labels)
plt.show()
class RunAntibert:
    def __init__(self, sequences) -> None:
        self.embeddings, self.attentions = self.run(sequences)
    def run(sequences):
        embeddings, attentions = antiberty.embed(sequences, return_attention=True)
        list_attentions = []
        for seq in attentions:
            seq = seq.numpy()
            list_attentions.append(seq)
        return embeddings, list_attentions
    
    
    
    
class PlotAttention:
    def __init__(self, attentions, sequences) -> None:
        self.attentions = attentions
        self.sequences = sequences
        self.plot()
    def plot(self):
        seaborn.heatmap(self.attentions[0], xticklabels=self.sequences[0], yticklabels=self.sequences[0])
        seaborn.heatmap(self.attentions[1], xticklabels=self.sequences[1], yticklabels=self.sequences[1])
    