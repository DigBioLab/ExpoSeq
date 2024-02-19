import pandas as pd
from ExpoSeq.plots.multiple_length_plot import LengthPlotMultiple
import matplotlib.pyplot as plt

def test_LengthPlotMultiple():
    sequencing_report = pd.read_csv(r"C:\Users\nilsh\my_projects\ExpoSeq\src\ExpoSeq\software_tests\test_files\sequencing_report.csv")
    region_of_interest = "aaSeqCDR3"
    fig = plt.figure(1)
    ax = fig.gca()
    font_settings = {'fontfamily': 'serif', 'fontsize': '18', 'fontstyle': 'normal', 'fontweight': 'bold'}
    LengthPlotMultiple(sequencing_report, region_of_interest, ax, font_settings)

    fig = plt.figure(1)
    ax = fig.gca()
    font_settings = {'fontfamily': 'serif', 'fontsize': '18', 'fontstyle': 'normal', 'fontweight': 'bold'}
    LengthPlotMultiple(sequencing_report, region_of_interest, ax, font_settings, plot_type="violin")