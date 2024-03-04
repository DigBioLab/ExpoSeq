from .collecting_all_arguments import ExpoSeqArgs, prep_args
from ..plots.rarefraction_curves import PrepareData  # import as package
import pandas as pd


def call_args():
    Args = ExpoSeqArgs()
    Args.add_save_csv()
    Args.add_sequencing_report()
    Args.add_samples(required=False)
    Args.add_region_plots()
    return Args


args = call_args()
parser = prep_args(args)


PrepData = PrepareData()
sequencing_report = pd.read_csv(parser.sequencing_report)
if (
    parser.samples == None or parser.samples[0] == ""
):  # enables one to choose per default all samples
    samples = sequencing_report["Experiment"].unique().tolist()
else:
    samples = parser.samples
results_plot = PrepData.tidy(
    sequencing_report,
    samples,
    parser.region_plots,
)

results_flat = results_plot.explode(["x_axis", "y_axis"])
results_flat.to_csv(parser.save_csv, index=False)
