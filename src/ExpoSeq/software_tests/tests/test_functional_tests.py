from src.ExpoSeq.pipeline import PlotManager
import pandas as pd

def test_pipeline():
    plot = PlotManager(experiment = "my_test_experiment", module_dir = "src/ExpoSeq/software_tests/test_files")

    plot.lengthDistribution_single()
    plot.lengthDistribution_multi()
    plot.print_antigens()
    plot.print_samples()
    plot.aa_distribution()
    plot.rarefraction_curves()
    plot.logoPlot_single()
    plot.logoPlot_multi()
    plot.rel_seq_abundance()
    plot.basic_cluster()
    plot.embedding_tsne(model = "sgt")
    plot.morosita_horn()
    plot.jaccard()
    plot.sorensen()
    plot.relative()
    plot.levenshtein_dendrogram()
    plot.alignment_quality()
    

