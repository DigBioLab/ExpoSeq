
from ExpoSeq.settings.collecting_all_arguments import ExpoSeqArgs, prep_args
import pandas as pd


def call_args():
    Args = ExpoSeqArgs()
    Args.add_sequencing_report()
    Args.add_region_plots()
    return Args


def get_main_information() -> list:

    sequencing_report = pd.read_csv(parser.sequencing_report)  
    samples = sequencing_report["Experiment"].unique().tolist()
    return samples
    
    
if __name__ == "__main__":
    args = call_args()
    parser = prep_args(args)
    samples = get_main_information()
    # Convert list to a space-separated string
    samples_str = ' '.join(samples)
    print(samples_str)