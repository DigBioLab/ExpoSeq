import pandas as pd
from src.ExpoSeq.plots.diversity_plot import DiversityPlot
import matplotlib.pyplot as plt

def test_diversity():
    
    sequencing_report_path = r"src/ExpoSeq/software_tests/test_files/test_show/sequencing_report.csv"
    sequencing_report = pd.read_csv(sequencing_report_path)
    fig = plt.figure(1, figsize = (12, 10))
    ax = fig.gca()
    DiversityPlot(sequencing_report, ax, method = "InverseSimpson")
    fig = plt.figure(1, figsize = (12, 10))
    ax = fig.gca()
    DiversityPlot(sequencing_report, ax, method = "Shannon")
    font_settings = {'fontfamily': 'serif', 'fontsize': '18', 'fontstyle': 'normal', 'fontweight': 'bold'}
    fig = plt.figure(1, figsize = (12, 10))
    ax = fig.gca()
    DiversityPlot(sequencing_report, ax, method = "InverseSimpson", font_settings=font_settings)
    fig = plt.figure(1, figsize = (12, 10))
    ax = fig.gca()
    DiversityPlot(sequencing_report, ax, method = "Shannon", font_settings=font_settings)


    
