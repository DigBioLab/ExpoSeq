import matplotlib.pyplot as plt
import pandas as pd
from sgt.SGT import SGT
from sklearn.decomposition import PCA
import numpy as np
from sklearn.manifold import TSNE


def cluster_toxins_tsne(fig, sequencing_report, sample, toxins, binding_data, toxin_names, pca_components, perplexity, iterations_tsne, font_settings, colorbar_settings, extra_figure):
    report = sequencing_report[sequencing_report["Experiment"] == sample]
    batch = report.groupby("Experiment").head(1000)
    if extra_figure == True:
        ax = fig.gca()
        fig2 = plt.figure()
        ax2 = fig2.gca()
    else:
        ax = fig.add_subplot(221)
        ax2 = fig.add_subplot(222)
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
    print("Explained variance after reducing to " + str(pca_components) + " dimensions:" + str(np.sum(pca.explained_variance_ratio_).tolist()))
    tsne = TSNE(n_components=2,
                verbose=0,
                perplexity=perplexity,
                n_iter=iterations_tsne)
    tsne_results = tsne.fit_transform(X)
    tsne_results = pd.DataFrame(tsne_results,
                                columns = [["tsne1", "tsne2"]])
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
    tsne_results["sequences"] = list(aminoacids)
    tsne_results['sequence_id'] = pd.Series(range(1170))
    tsne_plot = ax.scatter(tsne_results.tsne1,
                           tsne_results.tsne2,

                           c = tsne_results.binding)
    if toxin_names == True:
        for i, txt in enumerate(list(ids)):
            if list(kds)[i] > 0:
                ax.annotate(txt, (x_cor[i], y_cor[i]))
    else:
        pass
    plt.colorbar(tsne_plot, **colorbar_settings)
    ax.set_xlabel("t-SNE1", **font_settings)
    ax.set_ylabel("t-SNE2", **font_settings)
    original_fontsize = font_settings["fontsize"]
    font_settings["fontsize"] = 22
   # ax.set_title("Sequence Embedding on t-SNE Space for Sample " + sample, pad= 12, **font_settings)
    font_settings["fontsize"] = original_fontsize
    fig.tight_layout()


    ax2.scatter(tsne_results.tsne1, tsne_results.tsne2, alpha = 0.0)
    ax2.set_xlabel('t-SNE 1', **font_settings)
    ax2.set_ylabel('t-SNE 2', **font_settings)
    n = 0
    for j, row in tsne_results.iterrows():
        if row["binding"] > 1:
                ax2.text(row['tsne1'], row['tsne2'], row['sequence_id'], fontsize=10, weight = "bold")

        else:
            if n == 6:
                n = 0
                ax2.text(row['tsne1'], row['tsne2'], row['sequence_id'], fontsize=8)
        n += 1


    return tsne_results


#from gensim.models import Word2Vec
#import numpy as np

# Define CDR3 sequences
#cdr3_sequences = ['CASSLGAEQYF', 'CASSLGAEQYF']

# Define ProtVec model parameters
#embedding_size = 100
#window_size = 15
#min_word_count = 1
#workers = 8

# Train ProtVec model on CDR3 sequences
#cdr3_tokens = [list(seq) for seq in cdr3_sequences]
#protvec_model = Word2Vec(cdr3_tokens, size=embedding_size, window=window_size, min_count=min_word_count, workers=workers, sg=1)

# Generate ProtVec embeddings for CDR3 sequences
#cdr3_embeddings = np.array([protvec_model.wv[token] for token in cdr3_tokens])

# Print embeddings
#print(cdr3_embeddings)














