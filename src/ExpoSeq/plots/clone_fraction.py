import squarify
import matplotlib.pyplot as plt
import matplotlib
import pandas as pd
from textwrap import wrap

class VisFrac:
    
    def __init__(self, sequencing_report, ax,sample, visualize_sequences,prefered_cmap, top_clone_fraction, seqs_viz_fraction,alpha_val, pad_rectangles, force_reducing, limit_seq_filter, font_settings, region_string):
        self.ax = ax
        top_95_percent = self.get_top_fraction(report = sequencing_report, sample = sample,force_reducing = force_reducing, fraction = top_clone_fraction)
        self.seqs = self.prepare_labels(top_95_percent, region_string, limit_seq_filter, seqs_viz_fraction)
        self.fracs = top_95_percent["cloneFraction"]
        self.make_plot(pad = pad_rectangles, 
                       visualize_seqs= visualize_sequences,
                       alpha_val = alpha_val,
                       prefered_cmap=prefered_cmap,
                       )
        self.no_axis()
        self.add_title(font_settings, sample)
            
    @staticmethod
    def get_top_fraction(report,sample, force_reducing, fraction = 0.95):
        report = report[report["Experiment"] == sample]
        index = report["cloneFraction"].cumsum() < fraction
        top_95_percent = report[index]
        if top_95_percent.shape[0] > 1000:
            if force_reducing != None:
                top_95_percent = top_95_percent.iloc[:force_reducing]
        top_95_percent = top_95_percent.copy()
        top_95_percent["cloneFraction"] = top_95_percent["cloneFraction"] / fraction # normalize 
        return top_95_percent
    
    @staticmethod
    def prepare_labels(report,region_string,limit_seq_filter, fraction = 0.75):
        index = report["cloneFraction"].cumsum() < fraction
        seqs = report[region_string]
        seq_length = len(seqs)
        filtered_seqs = seqs[index]
        filtered_length = len(filtered_seqs)
        if filtered_length > limit_seq_filter:
            filtered_seqs = filtered_seqs.iloc[:limit_seq_filter]
            filtered_length = len(filtered_seqs)
        empty_values = pd.Series([""] * int((seq_length - filtered_length)), dtype="str")
        new_seqs = pd.concat([filtered_seqs, empty_values], ignore_index = True)

            
        return new_seqs
    
    @staticmethod
    def prep_cmap(top_readF, prefered_cmap = "Reds"):
        cmap = matplotlib.colormaps.get_cmap(prefered_cmap)
        mini=min(top_readF)
        maxi=max(top_readF)
        norm = matplotlib.colors.Normalize(vmin=mini, vmax=maxi)
        colors = [cmap(norm(value)) for value in top_readF] 
        return colors   
    

    def no_axis(self):
        self.ax.set_axis_off()
    
    def add_title(self, font_settings, sample):
        title = "\n".join(wrap("Clone fraction of " + sample, 40))
        self.ax.set_title(title,pad = 12, **font_settings)
        
    def make_plot(self, pad,visualize_seqs, alpha_val, prefered_cmap, **kwargs):
        if visualize_seqs == None:
            self.seqs = None
        else:
            pass
        colors = self.prep_cmap(self.fracs, prefered_cmap)
        squarify.plot(self.fracs, pad = pad, label = self.seqs, alpha = alpha_val, color = colors, ax = self.ax, text_kwargs={'fontsize':5}, **kwargs)
        
        
class Test:
    def __init__(self) -> None:
        
        sequencing_report = self.get_report()
        self.get_fig()
        region_string = "aaSeqCDR3"
        print("Region String must be full sequence")
        assert VisFrac.get_top_fraction(sequencing_report, "sample 0",force_reducing=1000, fraction = 0.95).shape[1] > 1
        assert VisFrac.get_top_fraction(sequencing_report, "sample 0",force_reducing=3, fraction = 0.95).shape[0] == 3
        assert VisFrac.prepare_labels(sequencing_report, region_string, fraction = 0.1).iloc[-1] == ""
        Vis = VisFrac(sequencing_report, self.ax, "sample 0", visualize_sequences = None, prefered_cmap = "Reds", top_clone_fraction = 0.95, seqs_viz_fraction = 0.75, alpha_val = 0.6, pad_rectangles = False,force_reducing=1000, region_string = region_string)
        assert Vis.seqs == None
        assert VisFrac(sequencing_report, self.ax, "sample 0", visualize_sequences = True, prefered_cmap = "Reds", top_clone_fraction = 0.95, seqs_viz_fraction = 0.75, alpha_val = 0.6, pad_rectangles = False,force_reducing=None, region_string = region_string)
        assert VisFrac(sequencing_report, self.ax, "sample 0", visualize_sequences = True, prefered_cmap = "Reds", top_clone_fraction = 0.95, seqs_viz_fraction = 0.75, alpha_val = 0.6, pad_rectangles = True,force_reducing= 1000, region_string = region_string)
        
    def get_fig(self):
        self.fig = plt.figure(1, figsize = (12, 10))
        self.ax = self.fig.gca()
        
    def get_report(self):
        sequencing_report = pd.read_csv("my_experiments/unitest/sequencing_report.csv")
        return sequencing_report
        
        

