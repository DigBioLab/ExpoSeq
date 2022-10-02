# USQ = unique sequences quality
from python_scripts.tidy_data import read_extract_data
from python_scripts.tidy_data.interpret_data import add_fraction
import matplotlib.pyplot as plt
import random
from collections import Counter
import warnings


divisible_by = 3
min_count = 1
tidy_data = read_extract_data.read_intermediate_reports(filename = sequencing_report_all,
                                                        input_pattern = local_pattern_more_digits)
tidy_data.filter_rows_not_divisible(divisor = divisible_by,
                                    column = 'lengthOfCDR3')
# if app insert a test and error report if someone fails to insert correct column name
tidy_data.filter_rows_on_min(column = 'cloneCount',
                            min_count = min_count)
sub_table = tidy_data.extract_substring_rows(lib_name = "Library_1") # change library here

sub_table = add_fraction(sequencing_report = sub_table)
tidy_data = read_extract_data.read_intermediate_reports(filename = sub_table,
                                                        input_pattern = local_pattern_more_digits)
tidy_data.nest_data(nest_by="Experiment",
                    column_to_list=['nSeqCDR3', 'clonefrac', 'cloneCount'],
                    unique=False)
sub_table = tidy_data.sequencing_report
x_axis = []
y_axis = []
for sample in range(len(sub_table)):
    sequences = list(sub_table.iloc[sample].nSeqCDR3)
    fraction = list(sub_table.iloc[sample].clonefrac)
    counter = 0
    temp_list = []
    unique_list = []
    total_list = []
    mean_count = (sum(sub_table.iloc[sample]["cloneCount"]) / 100)
    while counter < 100:
        temp = random.choices(sequences, weights=fraction, k = int(mean_count))  # randomly sample 100th of the sequences with the recalculated clonefractions as probablities
        temp_list.extend(temp)  ## extend the temp_list with the samples from each 'pull-out' of the 100 'pull-outs'
        unique = Counter(temp_list).keys()  ## count unique seequences from the increasing temp_list
        unique_list.append(unique)  ## for each run append the number of unique sequences to a new cell in this list - will be accumulative
        total_list.append(len(temp_list))  ## append the total sequences amount for plotting
        counter += 1
    plt.plot(total_list, [len(x) for x in unique_list], label = sub_table.iloc[sample].Experiment)
    x_axis.append(total_list)
    y_axis.append([len(x) for x in unique_list])

labels = sub_table.Experiment

if len(sub_table) > 7:
    warnings.warn("Risk of Overplotting: Many different colors are used. It may be hard to differ between the data")

for i in range(len(x_axis)):
    plt.plot(x_axis[i], y_axis[i])