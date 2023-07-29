import pandas as pd


def cleaning(sample_name, report, chosen_seq_length):
    sample = report[report["Experiment"] == sample_name]
    local_report = sample[["Experiment", "clonesFraction", "aaSeqCDR3"]]
    sequences = local_report.aaSeqCDR3
    #local_report = local_report.sort_values("cloneFraction").groupby("Experiment", as_index = False).head(3) # filters out three highest values of each group
    aminoacids = "ACDEFGHIKLMNPQRSTVWY"
    if chosen_seq_length == "median":
        chosen_sequence = int(sequences.str.len().median())#raise error if 0
    if type(chosen_seq_length) == int:
        chosen_sequence = int(chosen_seq_length)
    if chosen_seq_length == "max":
        chosen_sequence = int(sequences.str.len().max())
    if chosen_seq_length == "mean":
        chosen_sequence = int(sequences.str.len().mean())
    compDict = {aa: chosen_sequence*[0] for aa in aminoacids}
    sequences = local_report[local_report.aaSeqCDR3.astype(str).str.len() == chosen_sequence]["aaSeqCDR3"]
    length_filtered_seqs = sequences.shape[0]
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
    return aa_distribution, chosen_sequence, length_filtered_seqs


