from src.ExpoSeq.plots.levenshtein_clustering import PrepareData, LevenshteinClustering
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np 

def test_LevenshteinClustering():
    sequencing_report_path = r"src/ExpoSeq/software_tests/test_files/test_show/sequencing_report.csv"
    sequencing_report = pd.read_csv(sequencing_report_path)
    sequencing_report["cloneFraction"] = sequencing_report["readFraction"]  
    samples = ["GeneMind_1"]
    batch_size = 1000
    max_ld = 2
    min_ld = 0
    region_string = "aaSeqCDR3"
    Prep = PrepareData()
    G, report = Prep.calc_edges(sequencing_report, samples, batch_size, max_ld, min_ld, region_string)
    assert report.shape[0] > len(list(G.nodes)), "Single nodes were not removed"
    assert len(list(G.edges)) > len(list(G.nodes)), "No cluster had been formed"
    max_ld = 1
    G2, report2 = Prep.calc_edges(sequencing_report, samples, batch_size, max_ld, min_ld, region_string)
    assert len(list(G2.edges)) < len(list(G.edges)), "There must be less clusters in G2 because of the lower levenshtein distance"
    assert isinstance(G, nx.Graph)
    assert isinstance(report, pd.DataFrame)
    fig = plt.figure(1)
    ax = fig.gca()
    font_settings = {'fontfamily': 'serif', 'fontsize': '18', 'fontstyle': 'normal', 'fontweight': 'bold'}
    LevenshteinClustering(sequencing_report, samples, ax, region_string, max_ld, min_ld, batch_size, font_settings)
    grouped = sequencing_report.groupby('Experiment')
    # Initialize an empty list to store the top sequences for each sample
    top_sequences = []
    top_n = 10
    binding_data = pd.DataFrame(columns=[ 'aaSeqCDR3', 'Antigen 1'])
    # Iterate over each group (sample) and select the top 10 sequences
    for sample, group_data in grouped:
        top_sequences_per_sample = group_data['aaSeqCDR3'].value_counts().nlargest(top_n).index.tolist()
        values = np.random.randint(low=1000, high=1000000, size=len(top_sequences_per_sample))  
        sample_values = pd.DataFrame({'aaSeqCDR3': top_sequences_per_sample, 'Antigen 1': values})
        top_sequences.extend(top_sequences_per_sample)
        binding_data = pd.concat([binding_data, sample_values])
        
    G3, report3 = Prep.calc_edges(sequencing_report, samples, batch_size, max_ld, min_ld, region_string, binding_data, antigens = ["Antigen 1"])
    assert "Antigen 1" in report3.columns.tolist()
    binding_data = pd.DataFrame(data = {"aaSeqCDR3": "AAAA", "Antigen 1": 1}, index = [1])
    G3, report4 = Prep.calc_edges(sequencing_report, samples, batch_size, max_ld, min_ld, region_string, binding_data, antigens = ["Antigen 1"])
    assert "AAAA" not in list(G3.nodes)    # because it does not cluster
    fig = plt.figure(1)
    ax = fig.gca()
    font_settings = {'fontfamily': 'serif', 'fontsize': '18', 'fontstyle': 'normal', 'fontweight': 'bold'}
    LevenshteinClustering(sequencing_report, samples, ax, region_string, max_ld, min_ld, batch_size, font_settings, binding_data = binding_data, antigens = ["Antigen 1"], )
    