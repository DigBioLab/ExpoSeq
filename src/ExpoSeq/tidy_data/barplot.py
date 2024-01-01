
from pandas import DataFrame

def cleaning_data(all_alignment_reports, sequencing_report_all):
    all_alignment_reports = all_alignment_reports.sort_values("Input file(s)").reset_index(drop=True)
    unique_experiments = sequencing_report_all.sort_values("Experiment")["Experiment"].unique()
    unique_experiments = list(unique_experiments)
    tot_sequencing_reads = all_alignment_reports["Total sequencing reads"]
    aligned_reads = all_alignment_reports["Successfully aligned reads"].str.split("(", expand = True).iloc[:, 0]
    data = DataFrame([unique_experiments, aligned_reads, tot_sequencing_reads]).T
    data.columns = ["Experiment", "Aligned_Reads", "tot_sequenced_reads"]
    return data

