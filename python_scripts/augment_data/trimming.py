from numpy import array
def trimming(sequencing_report, divisible_by = 3, min_count = 3, new_fraction = "cloneFraction"):
    sequencing_report = sequencing_report[(sequencing_report["lengthOfCDR3"] % divisible_by) == 0]
    sequencing_report = sequencing_report[(sequencing_report["readCount"] > min_count)]
    new_column = sequencing_report['readCount'] / sequencing_report.groupby('Experiment')['readCount'].transform('sum')
    sequencing_report[new_fraction] = array(new_column)
    return sequencing_report