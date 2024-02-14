import matplotlib.pyplot as plt
import pandas as pd
from src.ExpoSeq.plots.stacked_aa_distribution import StackedAADistribution
import numpy as np 

def test_stacked_aa_distribution():
    sequencing_report_path = r"C:\Users\nilsh\my_projects\ExpoSeq\src\ExpoSeq\software_tests\test_files\sequencing_report.csv"
    sequencing_report = pd.read_csv(sequencing_report_path)
    
    aa_distribution = StackedAADistribution.prepare_data(sequencing_report, "sample 0", [0, 5], "aaSeqCDR3")
    assert aa_distribution.shape[1] == 20
    assert aa_distribution.shape[0] == 6 # region from index 0 to 5
    for j in range(5):
        assert aa_distribution.iloc[j, :].sum() > 0.99
        assert aa_distribution.iloc[j, :].sum() <= 1.01
        
    fig = plt.figure(1, figsize = (12, 10))
    ax = fig.gca()
    StackedAADistribution(sequencing_report, "sample 0", [0, 5], "aaSeqCDR3", protein = True, font_settings = {"fontsize": 20}, ax = ax)

    