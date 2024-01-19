from ExpoSeq.plots.rarefraction_curves import RarefractionCurves
import pandas as pd
import matplotlib.pyplot as plt

def test_rarefractioncurves():
    sequencing_report_path = r"src/ExpoSeq/software_tests/test_files/test_show/sequencing_report.csv"
    sequencing_report = pd.read_csv(sequencing_report_path)      
    fig = plt.figure(figsize=(10, 10))
    ax = fig.gca() 
    RarefractionCurves(sequencing_report, ["GeneMind_1", "GeneMind_2"], "aaSeqCDR3", ax = ax, font_settings = {"fontsize": 15}, legend_settings = {"fontsize": 15}, fraction_column="readFraction" )            
    RarefractionCurves(sequencing_report, ["GeneMind_1", "GeneMind_2"], "aaSeqCDR3", fraction_column="readFraction")
    plt.show() 

