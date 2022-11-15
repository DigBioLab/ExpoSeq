import re

def input_pattern():
    filepattern = input()
    return filepattern



def pattern_generator(filepattern): # works only if input is single digit
    filepattern_seperate = list(filepattern)
    local_pattern = ""
    local_pattern_more_digits = ""
    for letter in range(0, len(filepattern_seperate)):
        if filepattern_seperate[letter].isnumeric() != True:
            if filepattern_seperate[letter].isupper():
                local_pattern = local_pattern + '[A-Z]'
                local_pattern_more_digits = local_pattern_more_digits + '[A-Z]'
            elif filepattern_seperate[letter].islower():
                local_pattern = local_pattern + '[a-z]'
                local_pattern_more_digits = local_pattern_more_digits + '[a-z]'
            elif filepattern_seperate[letter].isspace():
                local_pattern = local_pattern + '\s'
                local_pattern_more_digits = local_pattern_more_digits + '\s'
            else:
                local_pattern = local_pattern + '[' + filepattern_seperate[letter] + ']'
                local_pattern_more_digits = local_pattern_more_digits + '[' + filepattern_seperate[letter] + ']'
        else:
            local_pattern_more_digits = local_pattern_more_digits + '\d+'
            local_pattern = local_pattern + '[1-9]'


    return local_pattern, local_pattern_more_digits





def grouping_files(local_pattern, local_pattern_more_digits, filenames): # will not work if you have three numbers in your filepattern on different places and your second number has two digits but your last number has one digit
    used_filenames = []
    grouped_filenames = []
    for index_filename in range(0, len(filenames)):
        subgroup = []
        if not bool(filenames[index_filename] in used_filenames):
            local_pattern_single = re.search(local_pattern,
                                             filenames[index_filename])
            local_pattern_multiple = re.search(local_pattern_more_digits,
                                               filenames[index_filename])
            if local_pattern_single is None:
               # local_pattern_multiple = re.search(local_pattern_more_digits,
                #                                   filenames[index_filename])
                local_pattern_single = local_pattern
                local_pattern_multiple = local_pattern_multiple.group()
            else:
                local_pattern_single = local_pattern_single.group()
                local_pattern_multiple = local_pattern_multiple.group()
            for filepath in filenames:
                if local_pattern_single == local_pattern_multiple:
                    if filepath not in used_filenames and local_pattern_single in filepath:
                        subgroup.append(filepath)
                        used_filenames.append(filepath)
                else:
                    if filepath not in used_filenames and local_pattern_multiple in filepath:
                        subgroup.append(filepath)
                        used_filenames.append(filepath)
            if bool(subgroup) != False:
                grouped_filenames.append(subgroup)
    return grouped_filenames


def find_exp_name(filename): # assumes that company appends uniqueCDR3_Exp_UniqueCDR3 to filenames so that they can be seperated by _
    list_filename = filename.split('\\')
    last_item_list = list_filename[-1]
    local_split = last_item_list.split("_")
    del local_split[len(local_split) - 3:]
    experiment_name = '_'.join(local_split)
    return experiment_name