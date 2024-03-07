from .collecting_all_arguments import ExpoSeqArgs, prep_args
from ..plots.matrices import jaccard_matrix, morosita_horn_matrix, relative_matrix, sorensen_matrix
import pandas as pd


def call_args():
    Args = ExpoSeqArgs()
    Args.add_save_csv()
    Args.add_sequencing_report()
    Args.add_samples(required = False, default = False)
    Args.add_region_plots()
    Args.add_matrix_type()
    return Args


args = call_args()
parser = prep_args(args)
sequencing_report = pd.read_csv(parser.sequencing_report)
if parser.matrix_type == "morosita_horn":
    matrix, unique_experiments = morosita_horn_matrix.PrepareData().cleaningPlot(parser.samples, sequencing_report, region_string = parser.region_plots, protein = True)
elif parser.matrix_type == "jaccard":
    matrix, unique_experiments =jaccard_matrix.PrepareData().cleaning_jaccard(sequencing_report, parser.region_plots, protein = True, specific_experiments = parser.samples)
elif parser.matrix_type == "sorensen":
    matrix, unique_experiments =sorensen_matrix.PrepareData().sorensen_matrix(sequencing_report, True, parser.region_plots, specific_experiments = parser.samples)
else:
    raise ValueError("Matrix type not recognized")

matrix.to_csv(parser.save_csv)