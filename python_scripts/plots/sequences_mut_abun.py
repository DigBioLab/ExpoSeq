

batch = sequencing_report[sequencing_report["Experiment"] == "Library_1_F2_2"]
clone_fractions = batch["clonesFraction"].to_list()
sequences = batch["aaSeqCDR3"]
sequences = sequences.to_list()
mutated_sequences = []
mutated_fractions = []

for i in range(len(sequences)):
    seq = sequences[i]
    for aa in range(len(seq)):
        new_seq = seq[:aa] + "*" + seq[(aa + 1):]
        mutated_sequences.append(new_seq)
        mutated_fractions.append(clone_fractions[i])

data = {"mutated_fractions": mutated_fractions, "mutated_sequences": mutated_sequences}
frame = pd.DataFrame(data)



for index, row in df.iterrows():
    for i, aa in enumerate(row['sequence']):
        for new_aa in mutation_list[i]:
            new_sequence = row['sequence'][:i] + new_aa + row['sequence'][i + 1:]
            new_clone_fraction = row['clone_fraction']
            mutated_sequences.append((new_sequence, new_clone_fraction))