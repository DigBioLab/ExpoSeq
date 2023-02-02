import editdistance
import pandas as pd

# create an empty dictionary to store the summarized clone fractions


# display the summarized dataframe
def tidy_levenshtein_fractions(sequencing_report,max_levenshtein_distance,length_filter, batch = 3000):
    sequencing_report = sequencing_report[sequencing_report['aaSeqCDR3'].str.len() >= length_filter]
    sequences = sequencing_report["aaSeqCDR3"].head(3000).to_list()
    clone_fractions = sequencing_report["clonesFraction"].head(3000).to_list()
    max_distance = max_levenshtein_distance
    summary = {}
    total_clone_fraction = 0
    # iterate over the sequences in the dataframe
    for i in range(len(sequences)):
        sequence = sequences[i]
        clone_fraction = clone_fractions[i]
        total_clone_fraction += clone_fraction
        # check if the sequence is already in the summary dictionary
        added = False
        for key in summary.keys():
            if editdistance.distance(sequence, key) <= max_distance:
                summary[key] += clone_fraction
                added = True
                break
        # if the sequence is not in the summary dictionary, add it
        if not added:
            summary[sequence] = clone_fraction
    # convert the summary dictionary to a pandas dataframe
    grouped = pd.DataFrame({'sequence': list(summary.keys()), 'clone_fraction': list(summary.values())})
    grouped['clone_fraction'] /= total_clone_fraction
    grouped = grouped.sort_values(by = "clone_fraction", ascending = False)
    return grouped

