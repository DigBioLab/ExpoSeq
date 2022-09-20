from python_scripts.augment_data.structure_files import grouping_files
from python_scripts.augment_data.structure_files import pattern_generator
from python_scripts.augment_data.loop_collect_reports import collect_intermediate_files

class activate_common_variables:

    def __init__(self, filepattern, filenames,):
        local_pattern, self.local_pattern_more_digits = pattern_generator(filepattern)
        self.grouped_filenames = grouping_files(local_pattern,
                                                self.local_pattern_more_digits,
                                                filenames)

    def create_seq_report(self):
        sequencing_report_all = collect_intermediate_files(self.grouped_filenames,
                                                            self.local_pattern_more_digits)
        return sequencing_report_all
        