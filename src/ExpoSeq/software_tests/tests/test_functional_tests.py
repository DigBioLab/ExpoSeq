from src.ExpoSeq.pipeline import PlotManager
import pandas as pd
import random 

def test_pipeline():
    
    plot = PlotManager(experiment = "multi_region", module_dir=r"C:\Users\nilsh\my_projects\ExpoSeq\src\ExpoSeq\software_tests\test_files", allow_binding_data=False, test_version=True)
    assert "targetSequences" in plot.avail_regions
    plot.change_region(region = "aaSeqCDR2")
    plot.change_region(region = "targetSequences")
    plot.change_filter(length_threshold_aa=2, min_read_count=2, remove_gaps = False, remove_errors = True)
    print(plot.region_string)
    seq_report = plot.sequencing_report
    assert "aaSeqtargetSequences" in seq_report.columns.to_list()
    assert "nSeqtargetSequences" in seq_report.columns.to_list()
    
    plot.print_antigens()
    plot.print_samples()
    assert plot.region_of_interest == "aaSeqtargetSequences"
    plot.discard_samples(["Pool1_R1_001"])
    seq_report = plot.sequencing_report
    assert plot.sequencing_report.shape[0] == 0
    
    plot = PlotManager(experiment = "test_show", module_dir = "src/ExpoSeq/software_tests/test_files", allow_binding_data=False, test_version=True,)

    plot.lengthDistribution_single()
    plot.lengthDistribution_multi()
    plot.print_antigens()
    plot.print_samples()
    plot.aa_distribution()
    plot.rarefraction_curves()
    plot.logoPlot()
    plot.logoPlot_multi()
    plot.rel_seq_abundance()
    plot.basic_cluster()
    plot.embedding_tsne(batch_size = 100)
    plot.morosita_horn()
    plot.jaccard()
    plot.sorensen()
    plot.relative()
    plot.levenshtein_dendrogram()
    plot.alignment_quality()
    plot.umap_sample_cluster(batch_size = 51)
    plot.umap_clustering_characteristic(batch_size = 51)
    
    plot = PlotManager(experiment = "test_show", module_dir = "src/ExpoSeq/software_tests/test_files", test_version=True, allow_binding_data=False)
    sequencing_report = plot.sequencing_report
    sample_0 = sequencing_report[sequencing_report["Experiment"] == "GeneMind_1"]
    aa = sample_0["aaSeqCDR3"]
    top_ten = aa.to_list()[:10]
    random_list = [random.randint(1, 3345832) for _ in range(10)]
    binding_data = pd.DataFrame(data = {"aaSeqCDR3": top_ten, "Antigen 0": random_list})
    plot.binding_data = binding_data
    plot.preferred_antigen= "Antigen 0"
    plot.dendro_bind(antigens = ["Antigen 0"])
    
   # plot.tsne_cluster_AG()
    plot.sample_diversity()
    plot.cluster_one_AG()
    plot.connect_samples()
    plot.embedding_network(batch_size = 100)
    
    plot = PlotManager(experiment = "test_show", module_dir = "src/ExpoSeq/software_tests/test_files", allow_binding_data="src/ExpoSeq/software_tests/test_files/binding_data.csv", test_version=True)
    plot.cluster_binding_data(batch_size = 100, iterations_tsne = 251)
    plot.ls_distance_binding()
    plot.tsne_cluster_AG(iterations_tsne = 251)
    sample_names = plot.experiments_list
    plot.lengthDistribution_single(sample = sample_names[0])
    samples_multi = [sample_names[0], sample_names[1]]
    plot.lengthDistribution_multi(samples = samples_multi)
    plot.aa_distribution(sample = sample_names[0], region = [3,7])
    plot.aa_distribution(sample = sample_names[0], region = [3,7], protein = False)
    plot.rarefraction_curves(samples = samples_multi)
    plot.logoPlot(sample = sample_names[0], highlight_specific_pos=4, chosen_seq_length = 10)
    plot.logoPlot(color_scheme = "hydrophobicity")
    plot.logoPlot(color_scheme = "charge")
    plot.logoPlot(color_scheme = "chemistry")
    plot.logoPlot(color_scheme = "NajafabadiEtAl2017")
    plot.logoPlot(color_scheme="dmslogo_charge")
    plot.logoPlot(color_scheme = "skylign_protein")
    plot.sorensen(specific_experiments = samples_multi)
    plot.rel_seq_abundance(sample = sample_names[0], alpha_val = 0.5, top_clone_fraction = 0.5)
    plot.cluster_binding_data_umap(batch_size = 51)
       

    