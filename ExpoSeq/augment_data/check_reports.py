
def check_completeness(all_alignment_reports, sequencing_report_all):
    all_alignment_reports = all_alignment_reports.sort_values("Input file(s)").reset_index(drop=True)
    unique_experiments = sequencing_report_all.sort_values("Experiment")["Experiment"].unique()
    unique_experiments = list(unique_experiments)
    not_found_strings = []
    for i in all_alignment_reports.iloc[:, 0]:
        found = False
        for j in unique_experiments:
            if j in i:
                found = True
        if found == False:
            print("Your data is imcomplete for the following sample: " + j)
            not_found_strings.append(j)


    for string in not_found_strings:
        all_alignment_reports = all_alignment_reports[~all_alignment_reports.iloc[:, 0].str.contains(string, na=False)]

    return all_alignment_reports