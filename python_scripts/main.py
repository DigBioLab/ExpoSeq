from python_scripts.augment_data.load_data import collectFiles
from python_scripts.augment_data.loop_collect_reports import collect_nocluster_files, load_alignment_reports
from python_scripts.tidy_data.trimming import trimming
from python_scripts.plots import usq_plot, plt_heatmap, barplot
from python_scripts.tidy_data import tidy_heatmap
from python_scripts.test_data.rename_labels import Library_2_to_panning
from python_scripts.plots.saveFig import saveFig
from python_scripts.augment_data.check_reports import check_completeness
def get_ready():
    filenames = collectFiles()
    sequencing_report_all, each_instance = collect_nocluster_files(filenames)
    all_alignment_reports = load_alignment_reports(filenames)
    sequencing_report_all, tidy_data = trimming(report = sequencing_report_all,
                                                divisible_by = 3,
                                                min_count = 1,
                                                new_fraction = "cloneFraction")
    all_alignment_reports = check_completeness(all_alignment_reports, sequencing_report_all)
    return sequencing_report_all, tidy_data, all_alignment_reports

tidy_data.map_to_column(object_to_map = Library_2_to_panning,
                        column_to_map = "Experiment",
                        new_column_name = "Panning Experiment")


unique_sequences, unique_experiments = tidy_heatmap.cleaning_data(protein = False,
                                                                  sequencing_report_input = sequencing_report_all,
                                                                  tidy_data = tidy_data)
matrix = plt_heatmap.plot_heatmap(unique_sequences, unique_experiments)
barplot.barplot(all_alignment_reports, sequencing_report_all)
usq_plot.plot_USQ(sequencing_report_all, local_pattern_more_digits)
saveFig()