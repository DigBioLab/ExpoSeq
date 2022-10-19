# USQ = unique sequences quality
from python_scripts.tidy_data import read_extract_data
from python_scripts.tidy_data.interpret_data import add_fraction

def cleaning_data(sequencing_report_all, local_pattern_more_digits):
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
    tidy_data.sequencing_report = sub_table
    #tidy_data = read_extract_data.read_intermediate_reports(filename = sub_table,
     #                                                       input_pattern = local_pattern_more_digits)
    tidy_data.nest_data(nest_by="Experiment",
                        column_to_list=['nSeqCDR3', 'clonefrac', 'cloneCount'],
                        unique=False)
    sub_table = tidy_data.sequencing_report
    return sub_table

