import numpy as np


length = []
for i, j in batch:
    length.append(len(i))
length = example["aaSeqCDR3"].str.len()
unique_length, counts_length = np.unique(np.array(length), return_counts = True)
plt.bar(unique_length, counts_length)
p_values = counts_length/np.sum(counts_length)
#normalized on one:
plt.bar(unique_length, counts_length/np.max(counts_length))
#p-values
plt.bar(unique_length, p_values)