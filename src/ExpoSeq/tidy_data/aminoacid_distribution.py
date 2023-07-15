from collections import Counter

import matplotlib.pyplot as plt
import numpy as np
from pandas import DataFrame, concat


sequencing_report = tidy_data.sequencing_report
third_round = ['EDE_K', 'ABA_K', 'ACA_K', 'FGF_L']
second_round = ['ED-_K', 'AB-_K', 'AC-_K', 'FG-_L']
first_round = ['E--_K', 'A--_K', 'F--_L']
rounds = [first_round, second_round, third_round]
width = 1/(len(rounds))
labels = ["First Round", "Second Round", "Third Round"]
n = 0
colors = ["orange", "royalblue", "limegreen"] # has to be changed
# works only for three panning rounds otherwise barplot will be confused
for round in rounds:
    sub_report = sequencing_report.loc[sequencing_report["Experiment"].isin(round)]
    sequences = list(sub_report["aaSeqCDR3"])
    counts = np.array(sub_report["cloneCount"])
    aminoacids = {"A": 0, "C": 0, "D": 0, "E": 0, "F": 0, "G": 0, "H": 0, "I": 0, "K": 0, "L": 0, "M": 0, "N": 0, "P": 0, "Q": 0, "R": 0, "S": 0, "T": 0,  "V": 0, "W": 0, "Y": 0}
    aa_count = DataFrame([aminoacids], columns = aminoacids.keys())

    for i in range(len(sequences)):
        aa = Counter(sequences[i])
        occurences = np.asarray(list(aa.values())) * counts[i]
        keys = aa.keys()
        local_dic = dict(zip(keys, occurences))
        local_dic = DataFrame([local_dic])
        aa_count = concat([aa_count, local_dic], axis =0).reset_index()
        aa_count.drop(columns = aa_count.columns[0], axis=1, inplace=True)
    summed_aa = aa_count.sum()
    summed_aa = np.array(summed_aa)
    normalized_aa = summed_aa/np.sum(summed_aa)
    aminoacids_keys = aa_count.columns
    X = np.arange(normalized_aa.shape[0])
    plt.bar(X + (0.3 - 0.3*n), normalized_aa, label = labels[n], width = 0.3, color = colors[n])
    n += 1
X = np.arange(normalized_aa.shape[0])
plt.xticks(X, aminoacids_keys)
plot_style = openParams('plot_style.txt')
plt.ylabel("relative Abundance", **plot_style)
plt.xlabel("Amino Acids", **plot_style)
params_legend = openParams("USQ_plot_legend_params.txt")
plt.legend(**params_legend)






sequences = sequencing_report_all['aaSeqCDR3'].apply(lambda x: list(x))
a = sequences.explode()