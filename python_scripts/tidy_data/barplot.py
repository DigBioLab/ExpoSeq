import re
from python_scripts.tidy_data.tidy_alignment_reports import collect_boxplot_data
import pandas as pd
from python_scripts.test_data.rename_labels import Library_2_to_panning



def cleaning_data(common_vars, local_pattern_more_digits):
    grouped_filenames = common_vars.grouped_filenames
    boxplot_data_frame = pd.DataFrame()
    report_name = "fastq_AlignmentReport"
    max_files = len(grouped_filenames)
    limit = int(input(f"how many do you want to plot? Max No. is {max_files}"))
    current_status = 0
     # not in class because it reduces flexibility
    not_used = []
    for cluster in grouped_filenames:
        if current_status == limit:
            break
        else:
            for file in cluster:
                filename = re.search(local_pattern_more_digits,
                                    file).group()
                if report_name in file:
                    alignment_frame = collect_boxplot_data(file,
                                                           local_pattern_more_digits)
                    boxplot_data_frame = pd.concat([boxplot_data_frame, alignment_frame],
                                                    ignore_index = False)
                    file_found = True
                    current_status += 1
                    break
                else:
                    if file == cluster[len(cluster)-1]:
                        not_used.append(filename)
                    pass
    experiment = boxplot_data_frame["Experiment"].map(Library_2_to_panning)
    boxplot_data_frame["Experiment"] = experiment
    return boxplot_data_frame
