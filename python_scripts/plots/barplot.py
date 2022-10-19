import re
from python_scripts.tidy_data.tidy_alignment_reports import collect_boxplot_data
import pandas as pd
import matplotlib.pyplot as plt



number_hierarchy = local_pattern_more_digits.count("\d")



## should be in a class where I can choose the funciton for the corresponding plot

# def boxplot_files
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
                if file == cluster[len(cluster)]:
                    not_used.append(filename)
                pass
