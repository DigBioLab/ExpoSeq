from python_scripts.augment_data.global_func_variables import activate_common_variables
from python_scripts.augment_data.load_data import collectFiles
from python_scripts.tidy_data.tidy_heatmap import heatmap_creator
import pyximport
#from python_scripts.plots.create_heatmap import loop_morosita_horn
import numpy as np
from python_scripts.augment_data.global_func_variables import activate_common_variables
from create_heatmap import loop_morosita_horn

filenames = collectFiles()
filepattern = "Library_2_F3_2"
common_vars = activate_common_variables(filenames = filenames, filepattern = filepattern)
sequencing_report_all = common_vars.create_seq_report()
create_heatmap = heatmap_creator(sequencing_report_all, filepattern)
sequencing_report_all, max_experiments = create_heatmap.tidy_for_dna()
nested_seq = np.array(sequencing_report_all.nSeqCDR3)
nested_clone_fraction = np.array(sequencing_report_all.cloneFraction)
heatmap_morosita_horn = loop_morosita_horn(nested_seq, nested_clone_fraction)



local_pattern_more_digits = common_vars.local_pattern_more_digits
