from src.ExpoSeq.plots.logo_plot import LogoPlot
import pandas as pd
import matplotlib.pyplot as plt

def test_logo_plot():
    fig = plt.figure(1, constrained_layout=True)
    ax = fig.gca()
    sequencing_report_path = r"src/ExpoSeq/software_tests/test_files/test_show/sequencing_report.csv"
    font_settings = {}
    sequencing_report = pd.read_csv(sequencing_report_path)
    sequencing_report["cloneFraction"] = sequencing_report["readFraction"]
    logo_instance = LogoPlot(ax, sequencing_report, region_string = "aaSeqCDR3", sample = "GeneMind_1", chosen_seq_length=None,
             highlight_spec_position=False, font_settings=font_settings, method = "", color_scheme = "hydrophobicity")
    logo_instance2 = LogoPlot(ax, sequencing_report, region_string = "aaSeqCDR3",
                              sample = "GeneMind_1", chosen_seq_length=5, highlight_spec_position=False, font_settings = font_settings, method = "", color_scheme="hydrophobicity")
    assert logo_instance2.chosen_seq_length == 5
    logo_instance3 = LogoPlot(ax, sequencing_report, region_string = "aaSeqCDR3",
                              sample = "GeneMind_1", chosen_seq_length=1,
                              highlight_spec_position=False, font_settings = font_settings, method = "", color_scheme="hydrophobicity")
    assert logo_instance3.chosen_seq_length == 14
    fig = plt.figure(1, constrained_layout=True)
    ax = fig.gca()
    font_settings = {'fontfamily': 'serif', 'fontsize': '18', 'fontstyle': 'normal', 'fontweight': 'bold'}
    logo_instance4 = LogoPlot(ax, sequencing_report, region_string = "aaSeqCDR3",
                              sample = "GeneMind_1", chosen_seq_length=1,
                              highlight_spec_position=False, font_settings = font_settings, method = "", color_scheme="hydrophobicity")

    fig = plt.figure(1, constrained_layout=True)
    ax = fig.gca()
    logo_instance5 = LogoPlot(ax, sequencing_report, region_string = "aaSeqCDR3",
                              sample = "GeneMind_1", chosen_seq_length=1,
                              highlight_spec_position=True, font_settings = font_settings, method = "", color_scheme="hydrophobicity")
