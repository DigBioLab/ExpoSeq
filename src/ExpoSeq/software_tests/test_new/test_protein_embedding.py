import pandas as pd
from src.ExpoSeq.plots.protein_embedding import PrepareData, PlotEmbedding
import matplotlib.pyplot as plt

def test_ProteintEmbedding():
    
    sequencing_report_path = r"src/ExpoSeq/software_tests/test_files/test_show/sequencing_report.csv"
    sequencing_report = pd.read_csv(sequencing_report_path)
    sequencing_report["cloneFraction"] = sequencing_report["readFraction"] 
    PlotEmbedding(sequencing_report=sequencing_report,
                model_choice="Rostlab/prot_bert", 
                list_experiments=["GeneMind_1"],
                region_of_interest="aaSeqCDR3",
                strands = False, 
                add_clone_size=None,
                batch_size = 200, 
                pca_components=70,
                perplexity=25, 
                iterations_tsne=1000,
                antigens = None)
    fig = plt.figure(1, figsize = (12, 10))
    ax = fig.gca()
    PlotEmbedding(sequencing_report=sequencing_report,
                model_choice="Rostlab/prot_bert", 
                list_experiments=["GeneMind_1"],
                region_of_interest="aaSeqCDR3",
                strands = False, 
                add_clone_size=None,
                batch_size = 200, 
                pca_components=70,
                perplexity=25, 
                iterations_tsne=1000,
                antigens = None,
                ax = ax,
                )