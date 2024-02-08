import pandas as pd
from src.ExpoSeq.plots.diversity_plot import DiversityPlot
import matplotlib.pyplot as plt

def test_diversity():
    
    sequencing_report_path = r"src/ExpoSeq/software_tests/test_files/test_show/sequencing_report.csv"
    sequencing_report = pd.read_csv(sequencing_report_path)
    sequencing_report["cloneFraction"] = sequencing_report["readFraction"]
    fig = plt.figure(1, figsize = (12, 10))
    ax = fig.gca()
    assert DiversityPlot(sequencing_report, region_of_interest="aaSeqCDR3",ax = ax, method = "InverseSimpson",)
    fig = plt.figure(1, figsize = (12, 10))
    ax = fig.gca()
    assert DiversityPlot(sequencing_report, ax= ax, method = "Shannon",  region_of_interest="aaSeqCDR3")
    font_settings = {'fontfamily': 'serif', 'fontsize': '18', 'fontstyle': 'normal', 'fontweight': 'bold'}
    fig = plt.figure(1, figsize = (12, 10))
    ax = fig.gca()
    assert DiversityPlot(sequencing_report, ax = ax, method = "InverseSimpson", font_settings=font_settings,  region_of_interest="aaSeqCDR3")
    fig = plt.figure(1, figsize = (12, 10))
    ax = fig.gca()
    assert DiversityPlot(sequencing_report, ax = ax, method = "Shannon", font_settings=font_settings,  region_of_interest="aaSeqCDR3")
    fig = plt.figure(1, figsize = (12, 10))
    ax = fig.gca()
    sequencing_report["aaSeqCDR1"] = sequencing_report["aaSeqCDR3"]
    sequencing_report["nSeqCDR1"] = sequencing_report["nSeqCDR3"]
    
    assert DiversityPlot(sequencing_report, ax = ax, method = "InverseSimpson", font_settings=font_settings,  region_of_interest="aaSeqCDR3"), "Plot does not work for other region"
    


    
