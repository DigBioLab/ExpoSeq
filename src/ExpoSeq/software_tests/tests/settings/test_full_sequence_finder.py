from src.ExpoSeq.settings.full_sequence_finder import FullSequence

def test_FullSequence():
    Finder = FullSequence(["CDR1", "FR2", "FR3", "CDR3", "FR4"])
    row_indexes = Finder.find_connecting_seq()
    assert type(row_indexes) == list
    assert row_indexes == ["FR3", "CDR3", "FR4"]
    Finder = FullSequence(["FR1", "CDR1", "FR3", "CDR3"])
    row_indexes = Finder.find_connecting_seq()
    assert row_indexes == ["FR3", "CDR3"]