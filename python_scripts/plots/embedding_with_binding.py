import matplotlib.pyplot as plt
import pandas as pd
from sgt.SGT import SGT
from sklearn.decomposition import PCA
import numpy as np
from sklearn.manifold import TSNE




def cluster_toxins_tsne(fig, ax,sequencing_report,sample, toxins, binding_data,font_settings, toxin_names, pca_components, perplexity, iterations_tsne, ):
    sequencing_report = sequencing_report[sequencing_report["Experiment"] == sample]
    batch = sequencing_report.groupby("Experiment").head(1000)
    mix = batch.merge(binding_data, on = "aaSeqCDR3", how = "left")
    mix = pd.concat([batch, binding_data])
    mix = mix.fillna(0)

    kds = mix[toxins].max(axis = 1)
    ids = mix[toxins].idxmax(axis = 1)
    corpus = pd.DataFrame([list(kds), mix.aaSeqCDR3.map(list)]).T
    sgt = SGT(kappa = 1, lengthsensitive = False)
    sequences = mix["aaSeqCDR3"].map(list)
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
    tsne_results = pd.DataFrame(tsne_results, columns = [["tsne1", "tsne2"]])

    x_cor = list(tsne_results.tsne1.iloc[:, 0])
    y_cor = list(tsne_results.tsne2.iloc[:, 0])

    tsne_results = pd.DataFrame(tsne_results,
                            columns=[["tsne1", "tsne2"]])
    aminoacids = mix["aaSeqCDR3"].to_list()
    experiments_batch = mix["Experiment"]
    unique_experiments_num = pd.factorize(experiments_batch)[0]
    tsne_results["experiments_factorized"] = list(unique_experiments_num)
    tsne_results["experiments_string"] = list(experiments_batch)
    tsne_results["binding"] = list(kds)
    tsne_plot = ax.scatter(tsne_results.tsne1, tsne_results.tsne2,c = tsne_results.binding)
    if toxin_names == True:
        for i, txt in enumerate(list(ids)):
            if list(kds)[i] > 0:
                ax.annotate(txt, (x_cor[i], y_cor[i]))
    else:
        pass
    plt.colorbar(tsne_plot)
    ax.set_xlabel("t-SNE1", **font_settings)
    ax.set_ylabel("t-SNE2", **font_settings)
    ax.set_title("Sequence Embedding on t-SNE Space for Sample " + sample)
    fig.tight_layout()