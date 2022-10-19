import numpy as np
from scipy.stats import norm
from python_scripts.plots.plot_params.open_txtfiles import openParams
import matplotlib.pyplot as plt
from python_scripts.augment_data.read_fastq import read_fastq

def gc_plot():
    sequences, quality = read_fastq()
    params_hist = openParams("gc_content_hist.txt")
    params_legend = openParams("USQ_plot_legend_params.txt")
    tick_params = openParams("tick_params.txt")
    plot_style = openParams("plot_style.txt")
    gc_content = []
    sequences = np.array(sequences,
                         dtype = "object")
    unique_seqs, counts = np.unique(sequences,
                                    return_counts = True)
    for seq in unique_seqs:
        count_g = seq.count('G')
        count_c = seq.count('C')
        seq_length = len(seq)
        gc_ratio = (count_g + count_c)/ seq_length
        gc_content.append(gc_ratio)
    gc_content = np.array(gc_content)*100
    mean_gc = np.mean(gc_content)
    variance_square = np.std(gc_content)

    plt.hist(x = gc_content,weights = counts, **params_hist)
    x_norm = np.arange(np.min(gc_content), np.max(gc_content), 1)
    plt.plot(x_norm, norm.pdf(x_norm, mean_gc, variance_square))
    plt.xlabel("GC-Content in %", **plot_style)
    plt.ylabel("Normalized Frequency of GC-Content", **plot_style)
    plt.tick_params(axis = 'both', **tick_params)
    #for key, value in tick_params.items():
     #   plt.gca().spines[value].set_visible(False)
    #plt.legend(**params_legend)