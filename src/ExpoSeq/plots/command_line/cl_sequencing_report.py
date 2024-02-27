from ExpoSeq.settings.reports import SequencingReport
from .collecting_all_arguments import ExpoSeqArgs, TestArgs

def call_args():
    Args = ExpoSeqArgs()
    Args.tsv_dir() 
    Args.region_of_interest() 
    Args.length_threshold() 
    Args.min_read_count_threshold()
    Args.remove_gaps()
    Args.remove_errors()
    Args.save_csv()
    return Args
    
args = call_args()
args.parser.parse_args()
tests = args.chosen_tests
TestInput = TestArgs(args)
for test in tests: # tests basically input values for chosen arguments. they should be all tested in the pipeline but this just prevents wrong inputs on a lower level.
    getattr(TestInput, test)
    
Report = SequencingReport(sequencing_report=args.tsv_dir)
Report.prepare_seq_report(region_string = args.region, 
                          length_threshold=args.length_threshold,
                          min_read_count=args.min_read_count,
                          remove_gaps = args.remove_gaps,
                          remove_errors=args.remove_errors
                          )
seq_report = Report.sequencing_report
seq_report.to_csv(args.save_csv)


