from python_scripts.tidy_data.read_extract_data import read_alignment_report
import pandas as pd

def collect_boxplot_data(report_filename, local_pattern):
    with open(report_filename) as f:
        alignment_report = f.readlines()
    data_extractor = read_alignment_report(alignment_report)
    digit = data_extractor.digit
    total_sequencing_reads = data_extractor.find_word_align_rep(word = "Total sequencing reads",
                                                                regex = digit) # digit in class
    successfully_aligned_reads = data_extractor.find_word_align_rep(word = "Succesfully aligned reads",
                                                                    regex = digit)
    experiment = data_extractor.find_word_align_rep("Output file",
                                                    regex = local_pattern)
    data = {'Experiment': experiment,
            'Aligned_Reads': successfully_aligned_reads,
            'tot_sequenced_reads': total_sequencing_reads}
    alignment_frame = pd.DataFrame(data)
    return alignment_frame