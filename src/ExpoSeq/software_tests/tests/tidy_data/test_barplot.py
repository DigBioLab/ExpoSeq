
from glob import glob
import os
from src.ExpoSeq.augment_data.loop_collect_reports import load_alignment_reports
import pandas as pd
from src.ExpoSeq.tidy_data.barplot import cleaning_data


def test_barplot():
    align_repo_path = os.path.join(r"src/ExpoSeq/software_tests/test_files/alignment_reports/*")
    alignment_reports = glob(align_repo_path)
    all_alignment_reports = load_alignment_reports(alignment_reports)
    sequecing_report_path =  os.path.join(r"src/ExpoSeq/software_tests/test_files/test_show/sequencing_report.csv")
    sequencing_report_all = pd.read_csv(sequecing_report_path)
    data = cleaning_data(all_alignment_reports, sequencing_report_all)
    unique_experiments = sequencing_report_all["Experiment"].unique().tolist()
    assert len(data["Experiment"].tolist()) == len(unique_experiments), "All alignment reports has more or less amount of experiments than sequencing report."
    assert data.columns.tolist() == ["Experiment", "Aligned_Reads", "tot_sequenced_reads"], "Columns are not correct."
    for i in range(len(data["Aligned_Reads"].tolist())):
        assert data["Aligned_Reads"].tolist()[i] <= data["tot_sequenced_reads"].tolist()[i], "Aligned reads are more than total sequenced reads."
