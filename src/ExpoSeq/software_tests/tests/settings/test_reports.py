from src.ExpoSeq.settings.reports import SequencingReport, ManageImportFiles
import pandas as pd
import os


def test_ManageTsv():
    dir_tsv = r"src\ExpoSeq\software_tests\test_files\test_show\tables_mixcr"
    TSVManager = ManageImportFiles()
    seq_report = TSVManager.merge_tsvs(dir_tsv)
    assert isinstance(seq_report, pd.DataFrame) == True
    cols_seq_report = seq_report.columns.to_list()
    assert "Experiment" in cols_seq_report
    assert "cloneFraction" in cols_seq_report
    assert "readCount" in cols_seq_report
    assert "cloneId" in cols_seq_report

    

def test_sequencingreportclass():
    report_path = os.path.join(r"src/ExpoSeq/software_tests/test_files/test_show/sequencing_report.csv")
    divisble_by = 3
    length_threshold = 3
    min_read_count = 10
    example_report = pd.read_csv(report_path)
    Report = SequencingReport(example_report)
    Report.prepare_seq_report(region_string = "CDR3",length_threshold = length_threshold, min_read_count=min_read_count, remove_gaps=False)
    test = Report.sequencing_report    
    assert Report.sequencing_report["aaSeqCDR3"].str.contains("_").any(), "There are no sequence errors in test file"
    Report.prepare_seq_report(region_string = "CDR3", length_threshold = length_threshold, min_read_count=min_read_count, remove_gaps=True)
    test = Report.sequencing_report
    assert not test["aaSeqCDR3"].str.contains("[*_]").any(), "There are still sequence errors in test file"
    experiment_names = test["Experiment"].unique().tolist()
    for exp in experiment_names:
        exp_specific = test[test["Experiment"] == exp]
        assert exp_specific["aaSeqCDR3"].shape[0] == exp_specific["aaSeqCDR3"].unique().shape[0], "Removing duplicates did not work"
    assert "aaSeqCDR3" in test.columns.tolist(), "Column not in sequencing report"
    assert "nSeqCDR3" in test.columns.tolist(), "Column not in sequencing report"
    assert test["aaSeqCDR3"].str.len().min() > length_threshold, "Sequence length filter does not work"
    assert test["readCount"].min() > min_read_count, "Read count filter does not work"
    tolerance = 0.000001
    for exp in experiment_names:
        sub_test = test[test["Experiment"] == exp]
        assert abs(sub_test["cloneFraction"].sum() - 1) < tolerance, "cloneFraction does not sum up to 1" # most important test, otherwise data is not tidied correctly
    unitest_dir = os.path.join("src", "ExpoSeq", "software_tests", "test_files","my_experiments", "unitest")
    if not os.path.isdir(unitest_dir):
        os.mkdir(unitest_dir)
    Report.check_sample_name(unitest_dir, "not_existent")
    assert not os.path.isfile(os.path.join(unitest_dir, "not_existend", "sequencing_report.csv"))
    Report.prepare_seq_report(region_string="targetSequences",  length_threshold = length_threshold, min_read_count=min_read_count)
    assert "nSeqtargetSequences" in Report.sequencing_report.columns.tolist()
    assert "aaSeqtargetSequences" in Report.sequencing_report.columns.tolist()
    assert "aaSeqCDR3" not in Report.sequencing_report.columns.tolist(), f"{Report.sequencing_report.columns.tolist()}"
    assert "aaSeqCDR3" in Report.origin_seq_report.columns.tolist()
    assert "targetSequences" in Report.origin_seq_report.columns.tolist()
    ### Platforma tests
    dir_tsv = r"src\ExpoSeq\software_tests\test_files\test_show\tables_mixcr"
    Report = SequencingReport(dir_tsv)
    Report.prepare_seq_report(region_string="targetSequences",  length_threshold = length_threshold, min_read_count=min_read_count, remove_gaps = False, remove_errors = True)
    
    
    


