import pandas as pd
from ..tidy_data.tidy_cluster_embedding import tidy_embed
from textwrap import wrap
def show_difference(sequencing_report, list_experiments,strands, batch_size, pca_components, perplexity, iterations_tsne,region_string, ax, legend_settings, font_settings):
    tsne_results, aminoacids, experiments_batch = tidy_embed(sequencing_report,
                                                              batch_size,
                                                              list_experiments,
                                                              pca_components,
                                                              perplexity,
                                                              iterations_tsne,
                                                              region_string)
    tsne_plot = ax.scatter(tsne_results.tsne1,
                           tsne_results.tsne2,
                           c = pd.factorize(experiments_batch)[0],
                            )
    ax.legend(handles=tsne_plot.legend_elements()[0],
               labels=list_experiments,
              **legend_settings)
    x = tsne_results["tsne1"].values.tolist()
    y = tsne_results["tsne2"].values.tolist()
    ax.set_xlabel("t-SNE1", **font_settings)
    ax.set_ylabel("t-SNE2", **font_settings)
    if strands == True:
        for i in range(0, len(x), 10):
            ax.annotate(aminoacids[i],
                        (x[i][0], y[i][0]),
                        fontsize = 5,
                        )
    original_fontsize = font_settings["fontsize"]
    font_settings["fontsize"] = 22
    title = "\n".join(wrap("t-SNE embedding of for given samples", 40))
    ax.set_title(title, pad= 12, **font_settings)
    font_settings["fontsize"] = original_fontsize