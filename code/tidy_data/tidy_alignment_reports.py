from read_extract_data import read_intermediate_reports


with open(report_filename) as f:
    alignment_report = f.readlines()
data_extractor = read_intermediate_reports(filename, alignment_report, )
total_sequencing_reads = data_extractor.find_word_align_rep(word = "Total sequencing reads",
                                                            regex = digit)
total_sequencing_reads = int(total_sequencing_reads)
successfully_aligned_reads = data_extractor.find_word_align_rep(word = "Succesfully aligned reads",
                                                                regex = digit)
successfully_aligned_reads = int(successfully_aligned_reads)
experiment = data_extractor.find_word_align_rep("Output file",
                                                regex = local_pattern_more_digits)

data = {'Experiment': experiment,
        'Aligned_Reads': successfully_aligned_reads,
        'tot_sequenced_reads': total_sequencing_reads}
alignment_frame = df.DataFrame(data)