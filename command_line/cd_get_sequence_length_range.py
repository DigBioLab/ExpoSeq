
from ExpoSeq.settings.collecting_all_arguments import ExpoSeqArgs, prep_args
import pandas as pd


def call_args():
    Args = ExpoSeqArgs()
    Args.add_sequencing_report()
    Args.add_region_plots()
    Args.add_single_sample()
    return Args

    
    
if __name__ == "__main__":
    args = call_args()
    parser = prep_args(args)
    report = pd.read_csv(parser.sequencing_report)
    filtered = report[report["Experiment"] == parser.single_sample]
    max_len = filtered[parser.region_plots].str.len().max()
    min_len = filtered[parser.region_plots].str.len().min()
    print(f"{max_len} {min_len}")