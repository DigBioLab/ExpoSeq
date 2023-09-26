
def cleaning_data(sequencing_report, samples):
    sub_table = sequencing_report.loc[sequencing_report["Experiment"].isin(samples)]
    apply_fun = lambda x: list(x)
    sub_table = sub_table.groupby("Experiment")[['nSeqCDR3', 'cloneFraction', 'readCount']].agg(apply_fun).reset_index()
    return sub_table

