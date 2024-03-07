from ExpoSeq.settings.collecting_all_arguments import ExpoSeqArgs, prep_args
from ExpoSeq.plots.length_distribution import LengthDistributionSingle
import pandas as pd


def call_args():
    Args = ExpoSeqArgs()
    Args.add_save_csv()
    Args.add_sequencing_report()
    Args.add_single_sample()
    Args.add_region_plots()
    return Args

if __name__ == "__main__":
    args = call_args()
    parser = prep_args(args)
    sequencing_report = pd.read_csv(parser.sequencing_report)
    # plot prepare
    unique_length, counts_length = LengthDistributionSingle.tidy(
        sequencing_report, parser.single_sample, parser.region_plots
    )
    pd.DataFrame({"unique_length": unique_length, "counts_length": counts_length}).to_csv(
        parser.save_csv
    )
