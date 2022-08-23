from code.augment_data.load_data import collectFiles
from code.tidy_data.read_extract_data import pattern_generator, grouping_files

filenames = collectFiles()
filepattern = input("normal file pattern")
local_pattern, local_pattern_more_digits = pattern_generator(filepattern)
grouped_filenames = grouping_files(local_pattern,
                                   local_pattern_more_digits,
                                   filenames)



