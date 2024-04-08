from src.ExpoSeq.plots.length_clone_fraction import LengthSeqFraction, PrepareData
import matplotlib.pyplot as plt
import pandas as pd

def test_seq_fraction():
    PrepData = PrepareData()
    region = "aaSeqCDR3"
    sample = "GeneMind_1"
    sequencing_report = pd.read_csv(r"C:\Users\nilsh\my_projects\ExpoSeq\src\ExpoSeq\software_tests\test_files\test_show\sequencing_report.csv")
    sequencing_report["cloneFraction"] = sequencing_report["readFraction"]
    sequencing_report, color = PrepData.tidy(sequencing_report, region, no_sequences = 10)
    assert len(sequencing_report["Experiment"].unique().tolist()) > 1, "You must have more than one sample in this output."
    assert len(sequencing_report["color"].unique().tolist()) == 11, "There must be eleven different colors since you chose 10 sequences and another one is for gray."
    fig = plt.figure(1)
    ax = fig.gca()
    legend_settings = {'loc': 'center left', 'bbox_to_anchor': (1, 0.5), 'fontsize': 9, 'frameon': True, 'framealpha': 1, 'facecolor': 'white', 'mode': None, 'title_fontsize': 'small'}
    font_settings = {'fontfamily': 'serif', 'fontsize': '18', 'fontstyle': 'normal', 'fontweight': 'bold'}
    LengthSeqFraction(sequencing_report, sample, region, ax = ax, legend_params=legend_settings, font_settings=font_settings )
    # Group DataFrame by sequence length and calculate the sum of clone fractions for each group