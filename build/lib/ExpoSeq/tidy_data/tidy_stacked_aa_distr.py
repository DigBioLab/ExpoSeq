import pandas as pd

def cleaning(sequencing_report, sample,region, region_string):

    sample = sequencing_report[sequencing_report["Experiment"] == sample]
    local_report = sample[["Experiment", "cloneFraction", region_string]]
    aminoacids = "ACDEFGHIKLMNPQRSTVWY"

    sequences = local_report[local_report[region_string].astype(str).str.len() >= region[1]][region_string]
    max_length = local_report[region_string].str.len().max()
    if not region[1] <= max_length:
        raise ValueError("you upper region limit is above the longest sequence. That is not possible. Please reduce it.")
    compDict = {aa: max_length*[0] for aa in aminoacids}
    for seq in sequences:
        for aa_position in range(len(seq)):
            aminoacid = seq[aa_position]
            if aminoacid == '*':
                pass
            else:
                compDict[aminoacid][aa_position] += 1
    aa_distribution = pd.DataFrame.from_dict(compDict)
    aa_distribution = aa_distribution.divide(aa_distribution.sum(axis=1), axis=0)
    aa_distribution = aa_distribution[(aa_distribution.index >= region[0]) & (aa_distribution.index <= region[1])]
    return aa_distribution