import pandas as pd
from ..tidy_data.tidy_cluster_embedding import tidy_embed
from textwrap import wrap
import pandas as pd
from sgt.SGT import SGT
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE

class ClusterEmbedding:
    def __init__(self) -> None:
        pass
    



    def tidy_embed(sequencing_report, batch_size, list_experiments, pca_components, perplexity ,iterations_tsne, region_string):
        report_batch = sequencing_report.groupby("Experiment").head(batch_size)
        selected_rows = report_batch.loc[report_batch["Experiment"].isin(list_experiments)]
        sgt = SGT(kappa=1,
                lengthsensitive=True)
        sequences = selected_rows[region_string].map(list)
        sequences_list = list(sequences)
        sgtembedding_df = sgt.fit_transform(corpus=sequences_list)
        pca = PCA(n_components=pca_components)
        pca.fit(sgtembedding_df)
        X = pca.transform(sgtembedding_df)
        print(np.sum(pca.explained_variance_ratio_))
        tsne = TSNE(n_components=2,
                    verbose=0,
                    perplexity=perplexity,
                    n_iter=iterations_tsne)
        tsne_results = tsne.fit_transform(X)
        tsne_results = pd.DataFrame(tsne_results,
                                    columns=[["tsne1", "tsne2"]])
        aminoacids = selected_rows["aaSeqCDR3"].to_list()
        experiments_batch = selected_rows.Experiment
        unique_experiments_num = pd.factorize(experiments_batch)[0]
        tsne_results["experiments_factorized"] = list(unique_experiments_num)
        tsne_results["experiments_string"] = list(experiments_batch)
        return tsne_results, aminoacids, experiments_batch
    
    
    

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
                           alpha = 0.5
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