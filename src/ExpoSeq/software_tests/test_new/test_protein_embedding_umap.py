import matplotlib.pyplot as plt
from src.ExpoSeq.plots.protein_embedding_umap import PrepareData, PlotEmbedding
import pandas as pd
import numpy as np


def test_plots():
    sequencing_report_path = r"src/ExpoSeq/software_tests/test_files/test_show/sequencing_report.csv"
    sequencing_report = pd.read_csv(sequencing_report_path)
    sequencing_report["cloneFraction"] = sequencing_report["readFraction"] 
    fig = plt.figure(1, figsize = (12, 10))
    ax = fig.gca()
    list_experiments = ["GeneMind_1"]
    font_settings = {'fontfamily': 'serif', 'fontsize': '18', 'fontstyle': 'normal', 'fontweight': 'bold'}
    legend_settings = {'loc': 'center left', 'bbox_to_anchor': (1, 0.5), 'fontsize': 9, 'frameon': True, 'framealpha': 1, 'facecolor': 'white', 'mode': None, 'title_fontsize': 'small'}
    
    #PrepData = PrepareData()
    #peptides, selected_rows, kds, ids = PrepData.tidy(sequencing_report, list_experiments, "aaSeqCDR3", batch_size = 80, characteristic = "hydrophobicity")

    Plot = PlotEmbedding(sequencing_report, ["GeneMind_1", "GeneMind_2"], "aaSeqCDR3", batch_size = 15, pca_components=10, ax = ax, strands = False, legend_settings=legend_settings, font_settings=font_settings)
    
    
def test_PlotEmbedding():
    sequencing_report_path = r"src/ExpoSeq/software_tests/test_files/test_show/sequencing_report.csv"
    sequencing_report = pd.read_csv(sequencing_report_path)
    sequencing_report["cloneFraction"] = sequencing_report["readFraction"] 
    fig = plt.figure(1, figsize = (12, 10))
    ax = fig.gca()
    list_experiments = ["GeneMind_1", ""]
    
    PrepData = PrepareData()
    peptides, selected_rows, kds, ids = PrepData.tidy(sequencing_report, list_experiments, "aaSeqCDR3", batch_size = 80)
    assert isinstance(PrepData.umap_results, pd.DataFrame)
    result_cols = PrepData.umap_results.columns
    assert "binding" not in result_cols
    assert "sequences" in result_cols
    assert "sequence_id" in result_cols
    assert "UMAP_1" in result_cols
    assert "UMAP_2" in result_cols
    umap_results = PrepData.umap_results
    assert umap_results.shape[0] == 80
    assert kds == None
    assert ids == None
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
        
    PrepBinding = PrepareData()
    peptides, selected_rows, kds, ids = PrepBinding.tidy(sequencing_report, list_experiments, region_of_interest="aaSeqCDR3", batch_size = 80, binding_data = binding_data, 
                                                    antigens=["Antigen 1"])
    umap_results = PrepBinding.umap_results
    result_cols = umap_results.columns
    assert "binding" in umap_results
    assert kds is not None
    assert ids is not None
    
    font_settings = {'fontfamily': 'serif', 'fontsize': '18', 'fontstyle': 'normal', 'fontweight': 'bold'}
    
    PlotEmbedding(sequencing_report=sequencing_report,
                model_choice="Rostlab/prot_bert", 
                list_experiments=["GeneMind_1", "GeneMind_2"],
                region_of_interest="aaSeqCDR3",
                strands = False, 
                add_clone_size=True,
                batch_size = 1000, 
                pca_components=70,
                n_neighbors=15,
                min_dist=0.3,
                random_seed=24,
                antigens = None, 
                ax = ax,
                font_settings = font_settings, 
)
    plt.show()