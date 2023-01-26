import matplotlib.pyplot as plt
import pandas as pd
from python_scripts.tidy_data.tidy_cluster_embedding import tidy_embed



def show_difference(sequencing_report, list_experiments,strands, batch_size, pca_components, perplexity, iterations_tsne, ax):
    tsne_results, aminoacids, experiments_batch = tidy_embed(sequencing_report,
                                                              batch_size,
                                                              list_experiments,
                                                              pca_components,
                                                              perplexity,
                                                              iterations_tsne)
    tsne_plot = ax.scatter(tsne_results.tsne1,
                           tsne_results.tsne2,
                           c = pd.factorize(experiments_batch)[0],
                            )
    ax.legend(handles=tsne_plot.legend_elements()[0],
               labels=list_experiments)
    x = tsne_results["tsne1"].values.tolist()
    y = tsne_results["tsne2"].values.tolist()
    if strands == True:
        for i in range(0, len(x), 10):
            ax.annotate(aminoacids[i],
                        (x[i][0], y[i][0]),
                        fontsize = 5,
                        )