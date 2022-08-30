from python_scripts.augment_data.loop_collect_reports import collect_intermediate_files
#from python_scripts.main import local_pattern_more_digits, grouped_filenames
from python_scripts.tidy_data import read_extract_data
from difflib import SequenceMatcher
import numpy as np

sequencing_report_all = collect_intermediate_files(grouped_filenames,
                                                   local_pattern_more_digits)
tidy_data = read_extract_data.read_intermediate_reports(filename = sequencing_report_all,
                                                        input_pattern = local_pattern_more_digits)
tidy_data.filter_rows_not_divisible(divisor = 3,
                                    column = 'lengthOfCDR3') # if app insert a test and error report if someone fails to insert correct column name
tidy_data.filter_rows_on_min(column = 'cloneCount',
                             min_count = 1)
sequencing_report_all = tidy_data.sequencing_report
unique_experiments = tidy_data.get_individuals("Experiment")
raw_sequences = np.array(sequencing_report_all["nSeqCDR3"])
experiments = np.array(sequencing_report_all["Experiment"])
max_experiments = unique_experiments.shape[0]
heatmap_axis = input("How many experiments do you want to plot? Max. No. {max_experiments}")
heatmap_absolute = np.zeros([heatmap_axis, heatmap_axis])
matches_collected = []
sum_matches = 0
for nt_sequence in raw_sequences:
    for i in range(0, len(raw_sequences)):
        current_experiment = experiments[i]
        if current_experiment != experiments[i+1]:
            number_matches = sum(a == b for a, b in zip(nt_sequence,
                                                        raw_sequences[i]))
            sum_matches += number_matches
            sum_matches = 0

        else:
            number_matches = sum(a==b for a, b in zip(nt_sequence,
                                                  raw_sequences[i]))
            sum_matches += number_matches


