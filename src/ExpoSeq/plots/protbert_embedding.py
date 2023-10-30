from textwrap import wrap
from ..tidy_data.tidy_protbert_embedding import TransformerBased
import pandas as pd

class Plot_Embedding:
    def __init__(self,ax, sequencing_report, list_experiments,strands,add_clone_size, batch_size, pca_components, perplexity, iterations_tsne, font_settings, legend_settings):
        self.ax = ax
        Transformer = TransformerBased(choice = "protbert")

        sequences, selected_rows = Transformer.filter_sequences(sequencing_report, batch_size=batch_size, experiments = list_experiments, region_of_interest="aaSeqCDR3", )


        X = Transformer.do_pca(sequences, batch_size, pca_components)
        peptides = selected_rows["aaSeqCDR3"].to_list()
        self.clones = selected_rows["cloneFraction"]
        tsne_results = Transformer.do_tsne(X, perplexity, iterations_tsne)
        self.create_plot(tsne_results, selected_rows)
        self.add_legend(list_experiments, **legend_settings)
        self.ax.set_xlabel("t-SNE1", **font_settings) # add font_settings
        self.ax.set_ylabel("t-SNE2", **font_settings)
        title = "\n".join(wrap("t-SNE embedding of for given samples", 40))
        self.ax.set_title(title, pad= 12, **font_settings)
        if strands == True:
            self.add_seq_anotation(peptides, tsne_results)
        if add_clone_size != None:
            self.add_size(add_clone_size)
            
        
    def create_plot(self, tsne_results, selected_rows):
        experiments_batch = selected_rows["Experiment"]
        self.tsne_plot = self.ax.scatter(tsne_results.tsne1,
                                    tsne_results.tsne2,
                                    c = pd.factorize(experiments_batch)[0],
                                    )
        
    def add_size(self, add_clone_size):
        size_points = self.clones * add_clone_size
        self.tsne_plot.set_sizes(size_points)
        
    def add_legend(self, list_experiments, legend_settings):
        self.ax.legend(handles=self.tsne_plot.legend_elements()[0],
            labels=list_experiments,
            **legend_settings)
        

    def add_seq_anotation(self,peptides, tsne_results):
        x = tsne_results["tsne1"].values.tolist()
        y = tsne_results["tsne2"].values.tolist()
        for i in range(0, len(x), 10):
            self.ax.annotate(peptides[i],
                        (x[i][0], y[i][0]),
                        fontsize = 5,
                        )
            