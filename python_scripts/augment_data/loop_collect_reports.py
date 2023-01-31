import pandas as pd
import re
from python_scripts.augment_data.structure_files import find_exp_name


def load_mixed_files(filenames):
    all_alignment_reports = pd.DataFrame([])
    sequencing_report = pd.DataFrame([])
    for i in filenames:
        if "AlignmentReport" in i:
            report = pd.read_table(i)
            splitted_report = report.iloc[:, 0].str.split(":",
                                                          expand=True)
            splitted_report = splitted_report.iloc[:, 0:2].T
            transposed_report = splitted_report.rename(columns=splitted_report.iloc[0]).drop(splitted_report.index[0])
            all_alignment_reports = pd.concat([all_alignment_reports, transposed_report])
        else:
            try:
                local_intermediate_report = pd.read_table(file)
                sequencing_report = pd.concat([sequencing_report, local_intermediate_report])
            except:
                print(i + "could not be read.")
    try:
        del all_alignment_reports["index"]
    except:
        pass
    return sequencing_report, all_alignment_reports





def load_alignment_reports(filenames):
    all_alignment_reports = pd.DataFrame()
    for i in filenames:
        if "AlignmentReport" in i:
            report = pd.read_table(i)
            splitted_report = report.iloc[:, 0].str.split(":",
                                                          expand = True)
            splitted_report = splitted_report.iloc[:, 0:2].T
            transposed_report = splitted_report.rename(columns=splitted_report.iloc[0]).drop(splitted_report.index[0])
            all_alignment_reports = pd.concat([all_alignment_reports, transposed_report])
    all_alignment_reports = all_alignment_reports.reset_index()
    try:
        del all_alignment_reports["index"]
    except:
        pass
    return all_alignment_reports