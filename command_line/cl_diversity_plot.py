from ExpoSeq.settings.collecting_all_arguments import ExpoSeqArgs, prep_args
from ExpoSeq.plots.diversity_plot import PrepareData
import pandas as pd


def call_args():
    Args = ExpoSeqArgs()
    Args.add_save_csv()
    Args.add_sequencing_report()
    Args.add_region_plots()
    Args.method_diversity()
    return Args

if __name__ == '__main__':
    args = call_args()
    parser = prep_args(args)


    sequencing_report = pd.read_csv(parser.sequencing_report)
    PrepData = PrepareData()
    values, unique_experiments = PrepData.cleaning(
        sequencing_report, parser.region_plots, parser.method
    )
    data_diversity = pd.DataFrame({"y_axis": values, "x_axis": unique_experiments})
    data_diversity.to_csv(parser.save_csv)