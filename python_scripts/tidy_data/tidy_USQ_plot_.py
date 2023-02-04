# USQ = unique sequences quality
from python_scripts.tidy_data.interpret_data import add_fraction
from test_data.rename_labels import Library_2_to_panning

def cleaning_data(sequencing_report, samples):
    sub_table = sequencing_report.loc[sequencing_report["Experiment"].isin(samples)]
    sub_table = add_fraction(sequencing_report = sub_table)
    apply_fun = lambda x: list(x)
    sub_table = sub_table.groupby("Experiment")[['nSeqCDR3', 'clonesFraction', 'readCount']].agg(apply_fun).reset_index()
    return sub_table

