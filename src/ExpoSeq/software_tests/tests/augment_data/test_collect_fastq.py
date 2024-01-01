from src.ExpoSeq.augment_data.collect_fastqs import CollectFastq
import os


def test_collectfastq():
    filepath = os.path.join(r"src/ExpoSeq/software_tests/test_files")
    filenames = CollectFastq(paired_end_sequencing=False).add_fastq_files(path_to_files=filepath)
    assert type(filenames)  == list, "filenames is not a list"
    assert type(filenames[0]) == str, "filenames[0] is not a str"
    for i in filenames:
        assert os.path.isfile(i) == True, "files are not added correctly"
    SingleEnd = CollectFastq(paired_end_sequencing=False)
    SingleEnd.get_files(path_to_forward=filepath)
    assert type(SingleEnd.paired) == list, "paired is not a list"
    assert type(SingleEnd.paired[0]) == list, "paired[0] is not a list"
    assert type(SingleEnd.paired[0][0]) == str, "paired[0][0] is not a str"
    assert os.path.isfile(SingleEnd.paired[0][0]) == True, "files are not added correctly"
    PairedEnd = CollectFastq(paired_end_sequencing=True)
    forward = os.path.abspath(os.path.join(filepath, r"paired/forward/"))
    backward = os.path.abspath(os.path.join(filepath, r"paired/backward/"))
    print(forward)
    PairedEnd.get_files(path_to_forward= forward, path_to_backward=backward)
    assert type(PairedEnd.paired) == list, "paired is not a list"
    assert type(PairedEnd.paired[0]) == list, "paired[0] is not a list"
    assert type(PairedEnd.paired[0][0]) == str, "paired[0][0] is not a str"
    assert os.path.isfile(PairedEnd.paired[0][0]) == True, "files are not added correctly"
    assert os.path.isfile(PairedEnd.paired[0][1]) == True, "files are not added correctly"


