from src.ExpoSeq.pipeline import PlotManager
import pandas as pd

def test_pipeline():
    plot = PlotManager(experiment = "my_test_experiment", module_dir = "src/ExpoSeq/software_tests/test_files", allow_binding_data=False, test_version=True,)

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
    
    plot = PlotManager(experiment = "my_test_experiment", module_dir = "src/ExpoSeq/software_tests/test_files", allow_binding_data="src/ExpoSeq/software_tests/test_files/binding_data.csv", test_version=True)
  
   # plot.tsne_cluster_AG()
    plot.sample_diversity()
    plot.cluster_one_AG()
    plot.connect_samples()
    plot.dendro_bind()
    plot.embedding_network()
    plot.embedding_tsne(model = "sgt")
    
    plot = PlotManager(experiment = "my_test_experiment", module_dir = "src/ExpoSeq/software_tests/test_files", allow_binding_data="src/ExpoSeq/software_tests/test_files/binding_data.csv", test_version=True)
    
    sample_names = plot.experiments_list
    plot.lengthDistribution_single(sample = sample_names[0])
    samples_multi = [sample_names[0], sample_names[1]]
    plot.lengthDistribution_multi(samples = samples_multi)
    plot.aa_distribution(sample = sample_names[0], region = [3,7])
    plot.aa_distribution(sample = sample_names[0], region = [3,7], protein = False)
    plot.rarefraction_curves(samples = samples_multi)
    plot.logoPlot_single(sample = sample_names[0], highlight_specific_pos=4, chosen_seq_length = 10)
    plot.logoPlot_single(color_scheme = "hydrophobicity")
    plot.logoPlot_single(color_scheme = "charge")
    plot.logoPlot_single(color_scheme = "chemistry")
    plot.logoPlot_single(color_scheme = "NajafabadiEtAl2017")
    plot.logoPlot_single(color_scheme="dmslogo_charge")
    plot.logoPlot_single(color_scheme = "skylign_protein")
    plot.dendro_bind(sample = sample_names[0], max_cluster_dist = 3, ascending = False) 
    plot.sorensen(specific_experiments = samples_multi)
    plot.rel_seq_abundance(sample = sample_names[0], alpha_val = 0.5, top_clone_fraction = 0.5)
       

