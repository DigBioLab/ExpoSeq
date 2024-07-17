from src.ExpoSeq.plots.protein_embedding_umap import PrepareData, PlotEmbedding
import pandas as pd
    
    
def test_different_models():
    list_experiments = ["GeneMind_1"]
    sequencing_report_path = r"src/ExpoSeq/software_tests/test_files/test_show/sequencing_report.csv"
    sequencing_report = pd.read_csv(sequencing_report_path)
    sequencing_report["cloneFraction"] = sequencing_report["readFraction"] 
    PrepData = PrepareData()
    peptides, selected_rows, kds, ids = PrepData.tidy(sequencing_report, list_experiments, "aaSeqCDR3", batch_size = 80, model_choice = r"Rostlab/prot_bert")
    PrepData = PrepareData()
    peptides, selected_rows, kds, ids = PrepData.tidy(sequencing_report, list_experiments, "aaSeqCDR3", batch_size = 80, model_choice = r"facebook/esm2_t6_8M_UR50D")
    PrepData = PrepareData()
    peptides, selected_rows, kds, ids = PrepData.tidy(sequencing_report, list_experiments, "aaSeqCDR3", batch_size = 80, model_choice = r"src\ExpoSeq\software_tests\test_files\nanobody_model")