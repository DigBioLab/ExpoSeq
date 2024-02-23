from src.ExpoSeq.plots.protein_embedding_umap import GetProteinProperty
import pandas as pd

def test_GetProteinProperty():
    sequencing_report_path = r"src/ExpoSeq/software_tests/test_files/test_show/sequencing_report.csv"
    sequencing_report = pd.read_csv(sequencing_report_path)
    sequences = sequencing_report["aaSeqCDR3"].head(10)
    Property = GetProteinProperty(sequences)
    Property.calc_attribute(attribute = "hydrophobicity")
    print(Property.sequence_property_interest)
    assert len(list(Property.sequence_property_interest.values()))
    Property = GetProteinProperty(sequences)
    Property.calc_attribute(attribute = "mass_charge_ratio")
    Property.calc_attribute(attribute = "isoelectric_point")
    Property.calc_attribute(attribute = "weight")