import pandas as pd

def cleaning(sequencing_report, sample, region, protein):
    if protein == True:
        sequence_type = "aaSeqCDR3"
    else:
        sequence_type = "nSeqCDR3"
    sample = sequencing_report[sequencing_report["Experiment"] == sample]
    local_report = sample[["Experiment", "clonesFraction", sequence_type]]
    sequences = local_report.aaSeqCDR3
    aminoacids = "ACDEFGHIKLMNPQRSTVWY"

    sequences = local_report[local_report["aaSeqCDR3"].astype(str).str.len() >= region[1]]["aaSeqCDR3"]
    max_length = local_report["aaSeqCDR3"].str.len().max()
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