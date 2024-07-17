from ExpoSeq.settings.collecting_all_arguments import ExpoSeqArgs, prep_args
from ExpoSeq.plots.logo_plot import PrepareData
import pandas as pd

def call_args():
    Args = ExpoSeqArgs()
    Args.add_save_csv()
    Args.add_sequencing_report()
    Args.add_chosen_seq_length()
    Args.add_samples()
    Args.add_region_plots()
    Args.add_method_logo()
    return Args


if __name__ == "__main__":
    args = call_args()
    parser = prep_args(args)
    sequencing_report = pd.read_csv(parser.sequencing_report)
    # plot prepare
    PrepData = PrepareData()
    local_report = sequencing_report.loc[sequencing_report["Experiment"].isin(parser.samples)]
    if parser.chosen_seq_length not in local_report["aaSeqCDR3"].astype(str).str.len().unique():
        print("The chosen sequence length is not in the available sequence lengths. Thus no processing will be conducted.")
        exit()
    else:
        aa_distribution = PrepData.cleaning(
            parser.samples,
            sequencing_report,
            parser.chosen_seq_length,
            parser.region_plots,
            parser.method_logo,
        )
        numbers_true = PrepData.get_labels(parser.region_plots, 
                                        parser.chosen_seq_length) # does not work for target sequence yet
        

        aa_distribution["labels"] = numbers_true
        color_schemes = ['chemistry', 'dmslogo_charge', 'dmslogo_funcgroup', 'skylign_protein', 'hydrophobicity', 'charge']
        cols = columns=aa_distribution.columns
        cols = cols.drop("labels")
        for scheme in color_schemes:
            color_list = PrepData.assign_color_to_table(aa_distribution, scheme)

            new_row_df = pd.DataFrame([color_list], columns = cols )
            # Append the new row to the DataFrame
            aa_distribution = pd.concat([aa_distribution, new_row_df], ignore_index=True)
        aa_distribution.to_csv(
            parser.save_csv, index = False
        )
