from ..settings.reports import SequencingReport 
from .collecting_all_arguments import ExpoSeqArgs, prep_args
import pandas as pd


 # If we want to let the user decide in each plot individually which region one wants to analyse then we need to choose the output from the loop.
    # Otherwise, the user can decide in the highest layer and the region for the plots will always be chosen automatically based on that input
def call_args():
    Args = ExpoSeqArgs()
    Args.add_tsv_dir() 
    Args.add_region_of_interest() 
    Args.add_length_threshold() 
    Args.add_min_read_count_threshold()
    Args.add_remove_gaps()
    Args.add_remove_errors()
    Args.add_save_csv()
    return Args
    
args = call_args()
parser = prep_args(args)
    
all_sequences_report = pd.DataFrame([])

Report = SequencingReport(sequencing_report=parser.tsv_dir) # first extract original sequencing_report (raw data)
Report.prepare_seq_report(region_string = parser.region, 
                          length_threshold=parser.length_threshold,
                          min_read_count=parser.min_read_count,
                          remove_gaps = parser.remove_gaps,
                          remove_errors=parser.remove_errors
                          )
orig_report_cols = set(Report.origin_seq_report.columns.tolist())
possible_regions = set(["aaSeqFR1", "aaSeqCDR1", "aaSeqFR2", "aaSeqCDR2", "aaSeqFR3", "aaSeqCDR3", "aaSeqFR4"])

avail_regions = list(orig_report_cols.intersection(possible_regions)) # get available aaSeq cols from your raw data
for region in avail_regions: # loop through all available regions and merge the tidied data in one report -> you will have Na's then
    region = region.replace("aaSeq", "")
    RegionReport = SequencingReport(sequencing_report=parser.tsv_dir)
    RegionReport.prepare_seq_report(region_string = parser.region, 
                          length_threshold=parser.length_threshold,
                          min_read_count=parser.min_read_count,
                          remove_gaps = parser.remove_gaps,
                          remove_errors=parser.remove_errors
                          )
    all_sequences_report = pd.concat([all_sequences_report, RegionReport.sequencing_report])

Report.sequencing_report.to_csv(parser.save_csv)


