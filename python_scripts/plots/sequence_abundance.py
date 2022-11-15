import numpy as np
from sklearn.preprocessing import StandardScaler
import pandas as pd
from difflib import SequenceMatcher
import itertools
import re
from random import sample
from sklearn.model_selection import train_test_split


def create_batch_on_length(batch_size, sequencing_report_all, lib = "Library_1_F2_2"):
    report_lib = sequencing_report_all[sequencing_report_all["Experiment"] == "Library_1_F2_2"]
    #number_positions = report_lib["lengthOfCDR3"].max()/3
    report_lib = sequencing_report_all[sequencing_report_all["Experiment"] == "Library_1_F2_2"]
    length = report_lib["aaSeqCDR3"].str.len()
    report_lib = report_lib.assign(length_aminoacid = length.values)
    grouped_fraction = report_lib.groupby("length_aminoacid")["cloneFraction"].sum()
    unique_length = np.unique(np.array(length))
    p_values = grouped_fraction/np.sum(np.array(grouped_fraction))
    zipped_values = list(zip(unique_length, p_values))
    unique_counts_dic = dict(zipped_values)
    example["p_values_length"] = report_lib["length_aminoacid"].map(unique_counts_dic)
    random_sequences = np.random.choice(list(report_lib["aaSeqCDR3"]), batch_size, list(report_lib["length_aminoacid"]))
    random_cloneFraction = report_lib[report_lib["aaSeqCDR3"].isin(list(random_sequences))]["cloneFraction"]
    random_cloneFraction = report_lib.shape[0]/random_sequences.shape[0] * np.array(random_cloneFraction)
    batch = list(zip(random_sequences, random_cloneFraction))

    return report_lib, batch

kmers = [3, 4, 5]


last_three = example["aaSeqCDR3"].str.strip().str[-3:]
unique, counts = np.unique(np.array(last_three), return_counts = True)

sorted_unique = unqiue_last_three[index_counts[::-1]]
sorted_counts = counts[index_counts[::-1]]

cols = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "T", "U", "V", "W", "Y"]
fraction = np.array(example["cloneFraction"])
sequences = sequencing_report_all["aaSeqCDR3"]
aa_pos_global = pd.DataFrame(columns = cols, index = range(31))


# first calculate the patterns with length 3,4,5 but as an example try first length 3
# then calculate the probability of having the pattern on the all siffrent positions per nucleotide
# include an option where certain patterns can be ignored since they might be too short or similar
# so create a distribution for each pattern which shows at which length the pattern can appear
# you can create for each pattern and sequence length a two dimensional distribution to find the probability for having the specific pattern in the given but also each other sequence with different length
# then plot a barplot or smth which shows you on the x axis the position and inside the bars the proportion of having the pattern there
#

# create pattern finder
#delete minus 6 from each sequence because all start and end aa are nearly the same
#pattern_length = 4
#threshold = 6
#local_frame = example[example["aaSeqCDR3"].str.len() > (pattern_length + threshold)]
#sequences_local = local_frame["aaSeqCDR3"]
#fraction = local_frame["cloneFraction"]
#batch_sequences, rest_seq, batch_fraction, rest_fraction = train_test_split(sequences_local, fraction, train_size = 1000, shuffle = True)
#del rest_seq, rest_fraction
#batch = list(zip(batch_sequences, batch_fraction))
# look only for patterns > 3
report_lib, batch = create_batch_on_length(batch_size = 200,
                                           sequencing_report_all = sequencing_report_all,
                                           )
pattern_collection = pd.DataFrame(columns=range(31))

n = 0
for peptide, fraction_one in batch:
    for peptide_two_index in range(n, len(batch)):
        peptide_two = batch[peptide_two_index][0]
        fraction_two = batch[peptide_two_index][1]
        match = SequenceMatcher(None, peptide, peptide_two).get_matching_blocks()
        for block in match:
            size_pattern = block.size
            if size_pattern > 2:
                index_pattern_one = block[0]
                index_pattern_two = block[1]
                patterns_found = pattern_collection.index

                local_pattern = peptide[index_pattern_one:(index_pattern_one+size_pattern)]
                sec_ind_one = index_pattern_one + size_pattern - 1
                sec_ind_two = index_pattern_two + size_pattern - 1


                if local_pattern in patterns_found:

                    pattern_collection.loc[local_pattern, index_pattern_one:sec_ind_one] += 1
                    pattern_collection.loc[local_pattern, index_pattern_two:sec_ind_two] += 1
                else:
                    pattern_collection.loc[local_pattern, index_pattern_one:sec_ind_one] = 1
                    pattern_collection = pattern_collection.fillna(0)
                    pattern_collection.loc[local_pattern, index_pattern_two:sec_ind_two] = 1

    n += 1


example_pattern = "DYW"
report_lib, batch = create_batch_on_length(batch_size = 2000,
                                           sequencing_report_all = sequencing_report_all,
                                           )
batch = pd.DataFrame(batch, columns = ["aaSeqCDR3", "fractions"])
length = batch["aaSeqCDR3"].str.len()
batch = batch.assign(length_aminoacid=length.values)
fractions_new = report_lib.shape[0]/batch.shape[0] * np.array(batch.fractions)
batch["fractions"] = fractions_new

kmers = 3
col_number = pattern_collection.shape[1]
patterns_to_search = []
for i in range(col_number-kmers):
    k_mers_frame = pattern_collection.iloc[:, i:i+kmers]
    summed_cols = k_mers_frame.sum(axis = 1)
    largest_patterns = summed_cols.nlargest(n = 3)
    patterns_to_search.append(list(largest_patterns.index))
patterns_to_search = list(itertools.chain.from_iterable(patterns_to_search))
patterns_to_search = np.unique(np.array(patterns_to_search))







    match = SequenceMatcher
match = SequenceMatcher(None, string1, string2).get_matching_blocks()
for i in match:
    size_pattern = i.size
    if size_pattern == pattern_length:


for seq in sequences:
    aa_pos = pd.DataFrame(columns=cols, index=range(31))
    seq_list = list(seq)
    length_sequence = len(seq_list)
    factor
    for i in range(len(seq_list)):
        aminoacid = seq_list[i]
        aa_pos.loc[i, aminoacid] = (1*(i+1)*fraction[n])/(length_sequence/2)
    aa_pos_global = pd.concat([aa_pos_global, aa_pos])
    n += 1

##distribution:
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