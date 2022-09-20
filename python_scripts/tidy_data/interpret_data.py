import pandas as pd

def add_fraction(sequencing_report,
                group = 'Experiment',
                percentage_from = 'cloneCount',
                new_column_name = 'clonefrac'):
    pd.options.mode.chained_assignment = None # otherwise error will be raised
    sequencing_report[new_column_name] = sequencing_report.groupby([group])[percentage_from].transform(lambda x: x / x.sum())
        # raise error if name exists

    return sequencing_report
