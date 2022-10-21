from python_scripts.augment_data.global_func_variables import activate_common_variables
from python_scripts.augment_data.load_data import collectFiles

import numpy as np
from python_scripts.augment_data.global_func_variables import activate_common_variables
from python_scripts.tidy_data.trimming import trimming
from python_scripts.augment_data import structure_files
from python_scripts.plots import usq_plot, plt_heatmap, barplot
from python_scripts.tidy_data import tidy_heatmap

from python_scripts.plots.saveFig import saveFig
filenames = collectFiles()
filepattern = "Library_2_F3_2"
common_vars = activate_common_variables(filenames = filenames,
                                        filepattern = filepattern)
sequencing_report_all = common_vars.create_seq_report()
local_pattern, local_pattern_more_digits = structure_files.pattern_generator(filepattern)
sequencing_report_all, tidy_data = trimming(report = sequencing_report_all,
                                            pattern = local_pattern_more_digits,
                                            divisible_by = 3,
                                            min_count = 1)#


usq_plot.plot_USQ(sequencing_report_all, local_pattern_more_digits)
unique_sequences, unique_experiments = tidy_heatmap.cleaning_data(protein = False,
                                                                  sequencing_report_input = sequencing_report_all,
                                                                  tidy_data = tidy_data)
matrix = plt_heatmap.plot_heatmap(unique_sequences, unique_experiments)
barplot.barplot(common_vars, local_pattern_more_digits)

saveFig()