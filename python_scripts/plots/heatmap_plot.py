from python_scripts.augment_data.loop_collect_reports import collect_intermediate_files
#from python_scripts.main import local_pattern_more_digits, grouped_filenames
from python_scripts.tidy_data import read_extract_data
import numpy as np
import seaborn as sns

sequencing_report_all = collect_intermediate_files(grouped_filenames,
                                                   local_pattern_more_digits)
tidy_data = read_extract_data.read_intermediate_reports(filename = sequencing_report_all,
                                                        input_pattern = local_pattern_more_digits)
tidy_data.filter_rows_not_divisible(divisor = 3,
                                    column = 'lengthOfCDR3')
# if app insert a test and error report if someone fails to insert correct column name
tidy_data.filter_rows_on_min(column = 'cloneCount',
                             min_count = 1)
tidy_data.nest_data(nest_by = "Experiment",
                    column_to_list = ['nSeqCDR3',
                                      'cloneCount'])

sequencing_report_all = tidy_data.sequencing_report
unique_experiments = tidy_data.get_individuals("Experiment")
max_experiments = unique_experiments.shape[0]
heatmap_axis = input("How many experiments do you want to plot? Max. No. {max_experiments}")
heatmap_axis = int(heatmap_axis)
heatmap_absolute = np.zeros([heatmap_axis, heatmap_axis])
heatmap_absolute = heatmap_absolute.astype(np.uint32)
matches_collected = []
sum_matches = 0

nested_seq = sequencing_report_all.nSeqCDR3


start = 0
index_x = 0
lib_matches = 0
for first_cluster in range(0, 37):
    for second_cluster in range(start, 37):
        for seq in nested_seq[first_cluster]:
            for seq2 in nested_seq[second_cluster]:
                if seq == seq2:
                    lib_matches += 1
                else:
                    pass
        heatmap_absolute[start, second_cluster] = lib_matches
        heatmap_absolute[second_cluster, start] = lib_matches

       # index_x += 1
        lib_matches = 0
    start += 1


nested_clone_counts = sequencing_report_all.cloneCount
array_clone_counts = np.array(nested_clone_counts)
nested_seq = sequencing_report_all.nSeqCDR3

start = 0
index_x = 0
lib_matches = 0
for first_cluster in range(0, 37):
    for second_cluster in range(start, 37):
        for seq in nested_seq[first_cluster]:
            for seq2 in nested_seq[second_cluster]:
                if seq == seq2:
                    lib_matches += 1
                else:
                    pass
        heatmap_absolute[start, second_cluster] = lib_matches
        heatmap_absolute[second_cluster, start] = lib_matches

       # index_x += 1
        lib_matches = 0
    start += 1




## normalizing
clone_counts = sequencing_report_all.cloneCount
def sum_clone_counnts
clone_counts = np.array(clone_counts)
summed_counts = []
for nest in clone_counts:
    sum = np.sum(nest)
    summed_counts.append(sum)

x = np.array(list(zip(summed_counts)))
y = np.array(summed_counts)
b = np.broadcast(x, y)
out = np.empty(b.shape)
out.flat = [u+v for (u,v) in b]

normalized_heatmap = heatmap_absolute/(out/2)
sns.heatmap(normalized_heatmap, cmap="Blues", annot=True, annot_kws={"size":4}, fmt='.2f')