import pandas as pd
from itertools import groupby
import re



def group_files(filenames):
    grouped_filenames = []
    used_filenames = []
    for index_filename in range(0, len(filenames)):
        substring = filenames[index_filename].split('.')[0]
        subgroup = []
        for filepath in filenames:
            if substring in filepath and filepath not in used_filenames:
                subgroup.append(filepath)
                used_filenames.append(filepath)
        if bool(subgroup) != False:
            grouped_filenames.append(subgroup)
    return grouped_filenames

def input_pattern():
    filepattern = input()
    return filepattern


def pattern_generator(filepattern):
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
            local_pattern_more_digits = local_pattern_more_digits + '[1-9]+'
            local_pattern = local_pattern + '[1-9]'


    return local_pattern, local_pattern_more_digits





def grouping_files(local_pattern, local_pattern_more_digits, filenames): # will not work if you have three numbers in your filepattern on different places and your second number has two digits but your last number has one digit
    used_filenames = []
    grouped_filenames = []
    for index_filename in range(0, len(filenames)):
        subgroup = []
        if not bool(filenames[index_filename] in used_filenames):
            try:
                local_pattern_single = re.search(local_pattern,
                                            filenames[index_filename]).group()
                local_pattern_multiple = re.search(local_pattern_more_digits,
                                                   filenames[index_filename]).group()
                for filepath in filenames:
                    if local_pattern_single != local_pattern_multiple:
                        if filepath not in used_filenames and local_pattern_single in filepath:
                            subgroup.append(filepath)
                            used_filenames.append(filepath)
                    else:
                        if filepath not in used_filenames and local_pattern_multiple in filepath:
                            subgroup.append(filepath)
                            used_filenames.append(filepath)
                if bool(subgroup) != False:
                    grouped_filenames.append(subgroup)
            except:
                pass
    return grouped_filenames












#def group_files(filenames):
 #   grouped_filenames = [list(group) for key, group in groupby(
  #                                      filenames,
   #                                     lambda string: string.split('.')[0])
    #                     ]
    #return grouped_filenames



#next steps: ideas:
# look for txt file in inner list
# extract sample name and save it in variable
# read txt file in pandas data table
# for the future: identify key variables for plotting and data analysis and tell whether they are present
# read other AlignmentReport and extract key variables and attach to table for each row
#

## read alignment report

if o.endswith('.fastq_AlignmentReport.txt')