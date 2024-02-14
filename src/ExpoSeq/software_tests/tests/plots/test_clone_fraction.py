from src.ExpoSeq.plots.clone_fraction import VisFrac
import pandas as pd
import matplotlib.pyplot as plt


def test_visualize_frac():      
        sequencing_report = pd.read_csv(r"C:\Users\nilsh\my_projects\ExpoSeq\src\ExpoSeq\software_tests\test_files\sequencing_report.csv")
        fig = plt.figure(1, figsize = (12, 10))
        ax = fig.gca()
        region_string = "aaSeqCDR3"
        print("Region String must be full sequence")
        assert VisFrac.get_top_fraction(sequencing_report, "sample 0",force_reducing=1000, fraction = 0.95).shape[1] > 1
        assert VisFrac.get_top_fraction(sequencing_report, "sample 0",force_reducing=3, fraction = 0.95).shape[0] == 3
        assert VisFrac.prepare_labels(sequencing_report, region_string,limit_seq_filter = 10, fraction = 0.1).iloc[-1] == ""
        Vis = VisFrac(sequencing_report, ax, "sample 0", visualize_sequences = False, prefered_cmap = "Reds", top_clone_fraction = 0.95, seqs_viz_fraction = 0.75, alpha_val = 0.6, pad_rectangles = False,force_reducing=1000, region_string = region_string, limit_seq_filter=10, font_settings = {"fontsize": 12})
        assert Vis.seqs == None
        assert VisFrac(sequencing_report, ax, "sample 0", visualize_sequences = True, prefered_cmap = "Reds", top_clone_fraction = 0.95, seqs_viz_fraction = 0.75, alpha_val = 0.6, pad_rectangles = False,force_reducing=None, region_string = region_string, limit_seq_filter=10, font_settings = {"fontsize": 12})
        assert VisFrac(sequencing_report, ax, "sample 0", visualize_sequences = True, prefered_cmap = "Reds", top_clone_fraction = 0.95, seqs_viz_fraction = 0.75, alpha_val = 0.6, pad_rectangles = True,force_reducing= 1000, region_string = region_string, limit_seq_filter=10, font_settings = {"fontsize": 12})
        
