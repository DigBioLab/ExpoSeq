from src.ExpoSeq.augment_data.binding_data import *


def test_binding_data():
    binding_data, second_prompt = open_file_binding(r"src\ExpoSeq\software_tests\test_files\sequencing_report.csv")
    assert second_prompt == True, "Second prompt should be True"
    assert binding_data == None, "Binding data should be None"
    binding_data, second_prompt = open_file_binding(r"src\ExpoSeq\software_tests\test_files\binding_data.csv")
    assert binding_data.columns[0] == "Sequences", "First column should be Sequences"
    assert second_prompt == False, "Second prompt should be False"
    assert binding_data.shape[0] == 128, "There should be 130 rows in the binding data"
    