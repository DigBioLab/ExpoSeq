from ExpoSeq.settings.collecting_all_arguments import ExpoSeqArgs, prep_args
from ExpoSeq.plots.logo_plot import PrepareData
import pandas as pd


def call_args():
    Args = ExpoSeqArgs()
    Args.add_save_csv()
    Args.add_sequencing_report()
    Args.add_chosen_seq_length()
    Args.add_samples()
    Args.add_region_plots()
    Args.add_method_logo()
    return Args


if __name__ == "__main__":
    args = call_args()
    parser = prep_args(args)
    sequencing_report = pd.read_csv(parser.sequencing_report)
    # plot prepare
    PrepData = PrepareData()
    aa_distribution = PrepData.cleaning(
        parser.samples,
        sequencing_report,
        parser.chosen_seq_length,
        parser.region_plots,
        parser.method_logo,
    )
    aa_distribution.to_csv(
        parser.save_csv
    )
