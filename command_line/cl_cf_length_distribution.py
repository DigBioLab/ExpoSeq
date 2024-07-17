from ExpoSeq.settings.collecting_all_arguments import ExpoSeqArgs, prep_args
from ExpoSeq.plots.length_clone_fraction import PrepareData
import pandas as pd


def call_args():
    Args = ExpoSeqArgs()
    Args.add_save_csv()
    Args.add_sequencing_report()
    Args.add_single_sample()
    Args.add_no_sequences_to_viz()
    Args.add_region_plots()
    return Args

if __name__ == '__main__':
    args = call_args()
    parser = prep_args(args)
    sequencing_report = pd.read_csv(parser.sequencing_report)
    report, colors = PrepareData().tidy(sequencing_report, parser.region_plots, parser.no_sequences_to_viz)
    report = report[report["Experiment"] == parser.single_sample]
    report.to_csv(parser.save_csv, index = False)