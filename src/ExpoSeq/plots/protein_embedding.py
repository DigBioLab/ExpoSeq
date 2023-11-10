from textwrap import wrap
from ..tidy_data.tidy_protbert_embedding import TransformerBased
import pandas as pd
import matplotlib.pyplot as plt

class Plot_Embedding:
    def __init__(self,ax, sequencing_report,model_choice, list_experiments,strands,add_clone_size, batch_size, pca_components, perplexity, iterations_tsne, font_settings, legend_settings, region_of_interest, binding_data = None, colorbar_settings = None, antigens = None, toxin_names = None, extra_figure = False):
        self.ax = ax
        self.binding_data = binding_data
        self.filter_binding_data(region_of_interest, antigens)
        Transformer = TransformerBased(choice = model_choice)
        sequences, selected_rows, selected_rows = Transformer.filter_sequences(sequencing_report, batch_size, list_experiments, self.binding_data, region_of_interest=region_of_interest, )

        X = Transformer.do_pca(sequences, batch_size, pca_components)
        peptides = selected_rows["aaSeqCDR3"].to_list()
        self.clones = selected_rows["cloneFraction"]
        self.tsne_results = Transformer.do_tsne(X, perplexity, iterations_tsne)
        if self.binding_data is not None:
            self.create_binding_plot(selected_rows, antigens, region_of_interest, toxin_names, colorbar_settings)
            if extra_figure == True:
                self.create_second_bind_plot(font_settings)
                title = "\n".join(wrap(f"t-SNE embedding for {antigens}", 40))
                self.ax.set_title(title, pad= 12, **font_settings)
            
        else:
            self.create_plot(self.tsne_results, selected_rows)
            self.add_legend(list_experiments, legend_settings)
            title = "\n".join(wrap("t-SNE embedding for given samples", 40))
            self.ax.set_title(title, pad= 12, **font_settings)
            
        self.ax.set_xlabel("t-SNE1", **font_settings) # add font_settings
        self.ax.set_ylabel("t-SNE2", **font_settings)

        if strands == True:
            self.add_seq_anotation(peptides, self.tsne_results)
        if add_clone_size != None:
            self.add_size(add_clone_size)
            
        
    def create_plot(self, tsne_results, selected_rows):
        experiments_batch = selected_rows["Experiment"]

        self.tsne_plot = self.ax.scatter(tsne_results.tsne1,
                                    tsne_results.tsne2,
                                    c = pd.factorize(experiments_batch)[0],
                                    alpha = 0.5,
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
            
    def filter_binding_data(self, region_of_interest, antigens):
        if self.binding_data is not None:
            merged_columns = [region_of_interest] + antigens
            self.binding_data = self.binding_data[merged_columns]
        
    def return_binding_results(self,  selected_rows, antigens, region_of_interest):
        kds = selected_rows[antigens].max(axis = 1)
        ids = selected_rows[antigens].idxmax(axis = 1)
        aminoacids = selected_rows[region_of_interest].to_list()
        experiments_batch = selected_rows["Experiment"]
        unique_experiments_num = pd.factorize(experiments_batch)[0]
        self.tsne_results["experiments_factorized"] = list(unique_experiments_num)
        self.tsne_results["experiments_string"] = list(experiments_batch)
        self.tsne_results["binding"] = list(kds)
        self.tsne_results["sequences"] = list(aminoacids)
        self.tsne_results['sequence_id'] = pd.Series(range(self.tsne_results.shape[0]))
        return kds, ids
    
    def create_second_bind_plot(self, font_settings):

        fig2 = plt.figure()
        ax2 = fig2.gca()
        ax2.scatter(self.tsne_results.tsne1, self.tsne_results.tsne2, alpha = 0.0)
        ax2.set_xlabel('t-SNE 1', **font_settings)
        ax2.set_ylabel('t-SNE 2', **font_settings)
        n = 0
        for j, row in self.tsne_results.iterrows():
            if row["binding"] > 1:
                    ax2.text(row['tsne1'], row['tsne2'], row['sequence_id'], fontsize=10, weight = "bold")

            else:
                if n == 6:
                    n = 0
                    ax2.text(row['tsne1'], row['tsne2'], row['sequence_id'], fontsize=8)
            n += 1
    
    def create_binding_plot(self, selected_rows, antigens, region_of_interest, toxin_names, colorbar_settings):
        kds, ids = self.return_binding_results(selected_rows, antigens, region_of_interest)
        self.tsne_plot = self.ax.scatter(self.tsne_results.tsne1,
                            self.tsne_results.tsne2,
                            c = self.tsne_results.binding,
                            alpha = 1,
                            cmap = "magma"
                            )
        x_cor = list(self.tsne_results.tsne1.iloc[:, 0])
        y_cor = list(self.tsne_results.tsne2.iloc[:, 0])
        if toxin_names == True:
            for i, txt in enumerate(list(ids)):
                if list(kds)[i] > 0:
                    self.ax.annotate(txt, (x_cor[i], y_cor[i]))
        else:
            pass
        plt.colorbar(self.tsne_plot, **colorbar_settings)
        


