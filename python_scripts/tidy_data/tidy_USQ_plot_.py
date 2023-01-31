# USQ = unique sequences quality
from python_scripts.tidy_data.interpret_data import add_fraction
from test_data.rename_labels import Library_2_to_panning

def cleaning_data(sequencing_report, lib_name):
    sub_table = sequencing_report[sequencing_report["Experiment"].str.contains(lib_name)]
    sub_table = add_fraction(sequencing_report = sub_table)
    apply_fun = lambda x: list(x)
    sub_table = sub_table.groupby("Experiment")[['nSeqCDR3', 'clonefrac', 'cloneCount']].agg(apply_fun).reset_index()
    experiment = sub_table["Experiment"].map(Library_2_to_panning) #delete in the end since it is very individual
    sub_table["Experiment"] = experiment
    return sub_table

