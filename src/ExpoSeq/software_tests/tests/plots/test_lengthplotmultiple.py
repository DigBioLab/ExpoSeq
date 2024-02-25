import pandas as pd
from src.ExpoSeq.plots.multiple_length_plot import LengthPlotMultiple
import matplotlib.pyplot as plt
from src.ExpoSeq.settings.reports import SequencingReport


def test_LengthPlotMultiple():
    sequencing_report_path = r"src/ExpoSeq/software_tests/test_files/test_show/sequencing_report.csv"
    sequencing_report = pd.read_csv(sequencing_report_path)
    fig = plt.figure(1)
    ax = fig.gca()
    font_settings = {'fontfamily': 'serif', 'fontsize': '18', 'fontstyle': 'normal', 'fontweight': 'bold'}
    LengthPlotMultiple(sequencing_report, "aaSeqCDR3", ax = ax, font_settings=font_settings)
    
    path = r"C:\Users\nilsh\my_projects\ExpoSeq\src\ExpoSeq\software_tests\test_files\my_experiments\multi_region\sequencing_report.csv"
    sequencing_report = pd.read_csv(path)
    Report = SequencingReport(sequencing_report)
    Report.prepare_seq_report("targetSequences", length_threshold = 3, min_read_count=0, remove_gaps=False )
    sequencing_report = Report.sequencing_report
    fig2 = plt.figure(2)
    ax2 = fig2.gca()
    LengthPlotMultiple(sequencing_report, "aaSeqtargetSequences", ax = ax2, font_settings=font_settings)
