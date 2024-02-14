import pandas as pd
from src.ExpoSeq.plots.tidy_protbert_embedding import TransformerBased

def test_transformerBased():
    sequencing_report_path = r"src/ExpoSeq/software_tests/test_files/test_show/sequencing_report.csv"
    sequencing_report = pd.read_csv(sequencing_report_path)
    sequencing_report["cloneFraction"] = sequencing_report["readFraction"] 
    Model = TransformerBased()
    sequences,sequences_filtered, selected_rows = Model.filter_sequences(sequencing_report, batch_size = 300, experiments = ["GeneMind_1"], binding_data = None)
    assert selected_rows.shape[0] == 300, "Shape of selected rows must equal to batch size"
    sequences,sequences_filtered, selected_rows = Model.filter_sequences(sequencing_report, batch_size = 3, experiments = ["GeneMind_1", "GeneMind_5"], binding_data = None)
    assert selected_rows.shape[0] == 6, "No duplicates for top 3 sequences for these samples, so selected rows must have 6 rows"
    sequences,sequences_filtered, selected_rows = Model.filter_sequences(sequencing_report, batch_size = 3, experiments = ["GeneMind_1", "GeneMind_2"], binding_data = None)
    assert selected_rows.shape[0] < 600, "Duplicates removed"
    example_seq = sequences[0]
    for i in range(1, len(example_seq), 2):
        assert example_seq[i] == " ", "Every second character in the sequence must be a space"
    binding_data = pd.DataFrame({"aaSeqCDR3": "CASSRLAGGTDTQYF", "antigen_1": 12345}, index= [1])
    sequences,sequences_filtered, selected_rows = Model.filter_sequences(sequencing_report, batch_size = 3, experiments = ["GeneMind_1", "GeneMind_5"], binding_data = binding_data)
    binding_values = selected_rows["antigen_1"]
    assert binding_values[0] == 12345
    binding_data = pd.DataFrame(data = {"aaSeqCDR3": ["CASSRLAGGTDTQYF", "ADDS"], "antigen_1": [12345,2]}, index = [1,2])
    sequences,sequences_filtered, selected_rows = Model.filter_sequences(sequencing_report, batch_size = 3, experiments = ["GeneMind_1", "GeneMind_5"], binding_data = binding_data)
    assert selected_rows["aaSeqCDR3"].isin(["ADDS"]).any(), "Outer merge failed because ADDS is only in the binding data and not in the other report. The merge failed"
    Model = TransformerBased()
    output_embed = Model.embedding_per_seq(sequences)
    assert output_embed.shape[0] == 7, "First dimension does not have 7 values for the 7 sequences"
    assert output_embed.shape[1] == 1024, "There are not 1024 dimensions in the second dimension"
    Model = TransformerBased("Rostlab/prot_bert")
    output_embed = Model.embedding_per_seq(sequences)
    assert output_embed.shape[0] == 7, "First dimension does not have 7 values for the 7 sequences"
    assert output_embed.shape[1] == 1024, "There are not 1024 dimensions in the second dimension"
    
    
        
    
    
    
    