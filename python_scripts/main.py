from augment_data import load_data
from augment_data.structure_files import input_pattern as pattern
from python_scripts.augment_data.structure_files import grouping_files
from python_scripts.augment_data.structure_files import pattern_generator
from python_scripts.augment_data.load_data import collectFiles

filenames = collectFiles()
filepattern = "Library_2_F3_2"
local_pattern, local_pattern_more_digits = pattern_generator(filepattern)
grouped_filenames = grouping_files(local_pattern,
                                   local_pattern_more_digits,
                                   filenames)



