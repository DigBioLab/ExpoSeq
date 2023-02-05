import numpy as np
import pandas as pd
from sgt.SGT import SGT
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE


def tidy_embed(sequencing_report, batch_size, list_experiments, pca_components, perplexity ,iterations_tsne):
    report_batch = sequencing_report.groupby("Experiment").head(batch_size)
    selected_rows = report_batch.loc[report_batch["Experiment"].isin(list_experiments)]
    sgt = SGT(kappa=1,
              lengthsensitive=True)
    sequences = selected_rows["aaSeqCDR3"].map(list)
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