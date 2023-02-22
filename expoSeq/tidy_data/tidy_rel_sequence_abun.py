import editdistance
import pandas as pd

# create an empty dictionary to store the summarized clone fractions


# display the summarized dataframe
def cleaning(sequencing_report,max_levenshtein_distance, samples, length_filter, batch = 3000):
    data = {}
    for sample in samples:
        report = sequencing_report[sequencing_report["Experiment"] == sample]
        report = report[report['aaSeqCDR3'].str.len() >= length_filter]
        sequences = report["aaSeqCDR3"].head(batch).to_list()
        clone_fractions = report["clonesFraction"].head(batch).to_list()
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
        count_data = {}
        # convert the summary dictionary to a pandas dataframe
        summary = sorted(summary.items(), key=lambda x: x[1], reverse=True)
        first_20 = dict(summary[:20])
        intermediate_df = pd.DataFrame.from_dict(first_20, orient='index')
        data[sample] = intermediate_df
    all_samples = pd.concat(data, axis=1)
    all_samples = all_samples.fillna(0)
    all_samples = all_samples.rename(columns={col: col[0].replace('/sample', "") for col in all_samples.columns})
    all_samples.columns = samples
    return all_samples





