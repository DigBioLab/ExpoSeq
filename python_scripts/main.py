from python_scripts.augment_data.global_func_variables import activate_common_variables
from python_scripts.augment_data.load_data import collectFiles

filenames = collectFiles()
filepattern = "Library_2_F3_2"
common_vars = activate_common_variables(filenames = filenames, filepattern = filepattern)
sequencing_report_all = common_vars.create_seq_report()
local_pattern_more_digits = common_vars.local_pattern_more_digits



