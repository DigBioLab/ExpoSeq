# oriented on tije.co/post/seqlogo/from_multiple_sequence_alignment/
import pandas as pd
import seqlogo

from python_scripts.tidy_data.interpret_data import mapFunc
from python_scripts.genetic_dogma import genetic_dogma
import matplotlib.pyplot as plt



def cleaning(sample_name, report):
    sample = report[report["Experiment"] == sample_name]
    local_report = sample[["Experiment", "cloneFraction", "aaSeqCDR3"]]
    sequences = local_report.aaSeqCDR3
    #local_report = local_report.sort_values("cloneFraction").groupby("Experiment", as_index = False).head(3) # filters out three highest values of each group
    aminoacids = "ACDEFGHIKLMNPQRSTVWY"
    longest_sequence = sequences.str.len().max() #raise error if 0
    compDict = {aa: longest_sequence*[0] for aa in aminoacids}
    for seq in sequences:
        for aa_position in range(len(seq)):
            aminoacid = seq[aa_position]
            if aminoacid == '*':
                pass
            else:
                compDict[aminoacid][aa_position] += 1
    aa_distribution = pd.DataFrame.from_dict(compDict)
    #aa_distribution = aa_distribution.divide(aa_distribution.shape[1])
    aa_distribution = aa_distribution.divide(aa_distribution.sum(axis = 1), axis = 0)
    return aa_distribution


