from ExpoSeq.settings.collecting_all_arguments import ExpoSeqArgs, prep_args
from ExpoSeq.plots.clone_fraction import VisFrac
import pandas as pd


def call_args():
    Args = ExpoSeqArgs()
    Args.add_save_csv()
    Args.add_sequencing_report()
    Args.add_single_sample()
    Args.add_limit_seq()
    Args.add_fraction()
    
    return Args


if __name__ == '__main__':
    args = call_args()
    parser = prep_args(args)
    sequencing_report = pd.read_csv(parser.sequencing_report)
    
    top_95_percent = VisFrac.get_top_fraction(sequencing_report, 
                             parser.single_sample,
                             parser.limit_seq,
                             parser.fraction
                             )
    cfs = top_95_percent["cloneFraction"].values.tolist()
    cf_id = list(range(1, len(cfs) + 1))
    pd.DataFrame([cf_id, cfs]).T.to_csv(parser.save_csv)
