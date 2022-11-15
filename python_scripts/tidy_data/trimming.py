from python_scripts.tidy_data.read_extract_data import read_intermediate_reports
from numpy import array
def trimming(report, divisible_by = 3, min_count = 3, new_fraction = "cloneFraction"):
    tidy_data = read_intermediate_reports(filename = report)
    tidy_data.filter_rows_not_divisible(divisor=divisible_by,
                                        column='lengthOfCDR3')
    # if app insert a test and error report if someone fails to insert correct column name

    tidy_data.filter_rows_on_min(column='cloneCount',
                                 min_count=min_count)
    report = tidy_data.sequencing_report
    new_column = report['cloneCount'] / report.groupby('Experiment')['cloneCount'].transform('sum')
    report[new_fraction] = array(new_column)
    tidy_data.sequencing_report = report
    return report, tidy_data