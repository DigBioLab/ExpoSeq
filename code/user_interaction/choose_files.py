import re
from tidy_alignment_reports import collect_boxplot_data
import pandas as pd

number_hierarchy = local_pattern_more_digits.count("\d")



## should be in a class where I can choose the funciton for the corresponding plot

# def boxplot_files
boxplot_data_frame = pd.DataFrame()
alignment_report = "fastq_AlignmentReport"
max_files = len(grouped_filenames)
limit = input(f"how many do you want to plot? Max No. is {max_files}")
current_status = 0
while current_status != limit: # not in class because it reduces flexibility
    for cluster in grouped_filenames:
        while file_found == False:
            for file in cluster:
                file_found = False
                filename = re.search(local_pattern_more_digits,
                                file).group()
                if alignment_report in file:
                    alignment_frame = collect_boxplot_data()
                    boxplot_data_frame = pd.concat([boxplot_data_frame, alignment_frame],
                                               ignore_index = False)
                    file_found = True
                    current_status += 1



