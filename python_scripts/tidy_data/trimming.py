from python_scripts.tidy_data.read_extract_data import read_intermediate_reports

def trimming(report, pattern, divisible_by = 3, min_count = 3):
    tidy_data = read_intermediate_reports(filename = report,
                                          input_pattern = pattern)
    tidy_data.filter_rows_not_divisible(divisor=divisible_by,
                                        column='lengthOfCDR3')
    # if app insert a test and error report if someone fails to insert correct column name

    tidy_data.filter_rows_on_min(column='cloneCount',
                                 min_count=min_count)
    sequencing_report_all = tidy_data.sequencing_report
    return sequencing_report_all, tidy_data