import pandas as pd
from src.ExpoSeq.plots.protein_embedding import PrepareData, PlotEmbedding
import matplotlib.pyplot as plt
import numpy as np 

def test_ProteinEmbedding():
    sequencing_report_path = r"src/ExpoSeq/software_tests/test_files/test_show/sequencing_report.csv"
    sequencing_report = pd.read_csv(sequencing_report_path)
    sequencing_report["cloneFraction"] = sequencing_report["readFraction"] 
    fig = plt.figure(1, figsize = (12, 10))
    ax = fig.gca()
    list_experiments = ["GeneMind_1"]
    PrepData = PrepareData()
    peptides, selected_rows, kds, ids = PrepData.tidy(sequencing_report, list_experiments, "aaSeqCDR3", batch_size = 80, iterations_tsne = 251)
    assert isinstance(PrepData.tsne_results, pd.DataFrame)
    result_cols = PrepData.tsne_results.columns
    assert "binding" not in result_cols
    assert "sequences" in result_cols
    assert "sequence_id" in result_cols
    assert "tsne1" in result_cols
    assert "tsne2" in result_cols
    tsne_results = PrepData.tsne_results
    assert tsne_results.shape[0] == 80
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
    peptides, selected_rows, kds, ids = PrepBinding.tidy(sequencing_report, list_experiments, region_of_interest="aaSeqCDR3", batch_size = 80, iterations_tsne = 251, binding_data = binding_data, 
                                                    antigens=["Antigen 1"])
    tsne_binding = PrepBinding.tsne_results
    result_cols = tsne_binding.columns
    assert "binding" in tsne_binding
    assert kds is not None
    assert ids is not None
    
    font_settings = {'fontfamily': 'serif', 'fontsize': '18', 'fontstyle': 'normal', 'fontweight': 'bold'}
    
    PlotEmbedding(sequencing_report=sequencing_report,
                model_choice="Rostlab/prot_bert", 
                list_experiments=["GeneMind_1"],
                region_of_interest="aaSeqCDR3",
                strands = False, 
                add_clone_size=None,
                batch_size = 100, 
                pca_components=70,
                perplexity=25, 
                iterations_tsne=251,
                antigens = None, 
                ax = ax,
                font_settings = font_settings)
    
    fig = plt.figure(1, figsize = (12, 10))
    ax = fig.gca()
    colorbar_settings = {'cmap': 'inferno', 'orientation': 'vertical', 'spacing': 'proportional', 'extend': 'neither'}

    PlotEmbedding(sequencing_report=sequencing_report,
                model_choice="Rostlab/prot_bert", 
                list_experiments=["GeneMind_1"],
                region_of_interest="aaSeqCDR3",
                strands = False, 
                add_clone_size=None,
                batch_size = 100, 
                pca_components=70,
                perplexity=25, 
                iterations_tsne=251,
                ax = ax,
                font_settings=font_settings,
                extra_figure=True,
                antigens = ["Antigen 1"],
                binding_data=binding_data,
                colorbar_settings=colorbar_settings
                )
    
    plt.show()