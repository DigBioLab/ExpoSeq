import pandas as pd

def add_fraction(sequencing_report,
                group = 'Experiment',
                percentage_from = 'readCount',
                new_column_name = 'clonefrac'):
    pd.options.mode.chained_assignment = None # otherwise error will be raised
    sequencing_report[new_column_name] = sequencing_report.groupby([group])[percentage_from].transform(lambda x: x / x.sum())
        # raise error if name exists
    return sequencing_report


def mapFunc(sequencing_report, column, func, column_name):
    mapped_column = sequencing_report[column].apply(func)
    sequencing_report[column_name] = mapped_column
    return sequencing_report

