from ExpoSeq.settings.collecting_all_arguments import ExpoSeqArgs, prep_args
from ExpoSeq.plots.hist_lvst_dist import PrepareData
import pandas as pd


def call_args():
    Args = ExpoSeqArgs()
    Args.add_save_csv()
    Args.add_sequencing_report()
    Args.add_samples()
    Args.add_batch_size()
    Args.add_levenshtein_distance()
    Args.add_region_plots()
    return Args

if __name__ == "__main__":
    args = call_args()
    parser = prep_args(args)
    sequencing_report = pd.read_csv(parser.sequencing_report)
    # plot prepare
    PrepData = PrepareData()
    linked, aa_clustered = PrepData.tidy(
        sequencing_report, parser.samples, parser.batch_size, parser.region_plots, parser.levenshtein_distance
    )
    assert linked.shape[0] > 0, "No sequences found for the given sample"
    ## A warning or something which indicates that there are too many nodes what be useful, because the plot becomes messy otherwise.
    distance = linked[:, 2]
    cluster_1 = linked[:, 0]
    cluster_2 = linked[:, 1] # there are higher numbers for cluster id than sequences avaialbale because new clusteres are formed by mergin existing clusters during the hierarchical clustering process
    
    pd.DataFrame({"distance": distance, "cluster_id1": cluster_1, "cluster_id2": cluster_2}).to_csv(
        parser.save_csv
    )
    
