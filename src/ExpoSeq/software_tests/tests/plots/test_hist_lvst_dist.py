import pandas as pd
from src.ExpoSeq.plots.hist_lvst_dist import LevenshteinDend, PrepareData
import matplotlib.pyplot as plt
import pytest

def test_LevenshteinDend():
    sequencing_report_path = r"src/ExpoSeq/software_tests/test_files/test_show/sequencing_report.csv"
    sequencing_report = pd.read_csv(sequencing_report_path)
    linked, clustered = PrepareData().tidy(sequencing_report, ["GeneMind_1"], 200, "aaSeqCDR3", 3)
    assert linked.shape[0] == len(clustered) - 1
    assert linked.shape[1] == 4
    LevenshteinDend(sequencing_report, "aaSeqCDR3", ["GeneMind_1"], batch_size = 200, max_cluster_dist = 3)
    fig = plt.figure(1, figsize = (12, 10))
    ax = fig.gca()
    LevenshteinDend(sequencing_report, "aaSeqCDR3", ["GeneMind_1"], batch_size = 200, max_cluster_dist = 3, ax = ax)

        
    


