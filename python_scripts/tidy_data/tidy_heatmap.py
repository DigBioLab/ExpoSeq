from python_scripts.tidy_data import read_extract_data


def cleaning(local_pattern_more_digits, sequencing_report_all, min_count = 1, divisible_by = 3):

    tidy_data = read_extract_data.read_intermediate_reports(filename = sequencing_report_all,
                                                            input_pattern = local_pattern_more_digits)
    tidy_data.filter_rows_not_divisible(divisor = divisible_by,
                                        column = 'lengthOfCDR3')
# if app insert a test and error report if someone fails to insert correct column name
    tidy_data.filter_rows_on_min(column = 'cloneCount',
                                min_count = min_count)
    tidy_data.nest_data(nest_by = "Experiment",
                        column_to_list = ['nSeqCDR3', 'cloneFraction'],
                        unique = False)
    sequencing_report_all = tidy_data.sequencing_report
    unique_experiments = tidy_data.get_individuals("Experiment")
    max_experiments = unique_experiments.shape[0]
    return sequencing_report_all, max_experiments


