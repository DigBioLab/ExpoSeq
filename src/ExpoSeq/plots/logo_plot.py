import pandas as pd
import numpy as np
def calculate_entropy(probs):
    """Calculate Shannon entropy."""
    return -np.sum([p * np.log2(p) if p > 0 else 0 for p in probs])
def calculate_bit_values(entropies):
  """Calculates the bit values for a sequence logo, given a list of Shannon entropies."""

  bit_values = []
  for entropy in entropies:
    bit_value = 2 - entropy
    if bit_value < 0:
      bit_value = 0
    bit_values.append(bit_value)
  return bit_values
def cleaning(sample_name, report, chosen_seq_length, region_string, method):
    sample = report[report["Experiment"] == sample_name]
    local_report = sample[["Experiment", "cloneFraction", region_string]]
    sequences = local_report[region_string]
    #local_report = local_report.sort_values("cloneFraction").groupby("Experiment", as_index = False).head(3) # filters out three highest values of each group
    aminoacids = "ACDEFGHIKLMNPQRSTVWY"

    compDict = {aa: chosen_seq_length*[0] for aa in aminoacids}
    sequences = local_report[local_report.aaSeqCDR3.astype(str).str.len() == chosen_seq_length]["aaSeqCDR3"]
    length_filtered_seqs = sequences.shape[0]
    for seq in sequences:
        for aa_position in range(len(seq)):
            aminoacid = seq[aa_position]
            if aminoacid == '*':
                pass
            else:
                compDict[aminoacid][aa_position] += 1
    if method == "bits":
        frequencies = {aa: [count / length_filtered_seqs for count in compDict[aa]] for aa in aminoacids}
        # Calculate Shannon entropy for each position
        entropies = [calculate_entropy([frequencies[aa][i] for aa in aminoacids]) for i in range(chosen_seq_length)]
        # Calculate bits for sequence logo for each amino acid
        bits_dict = {aa: [2 - entropies[i] if frequencies[aa][i] > 0 else 0 for i in range(chosen_seq_length)] for aa in aminoacids}
        aa_distribution = pd.DataFrame.from_dict(bits_dict)
    else:
        aa_distribution = pd.DataFrame.from_dict(compDict)
        aa_distribution = aa_distribution.divide(aa_distribution.sum(axis=1), axis=0)
    return aa_distribution, chosen_seq_length, length_filtered_seqs


