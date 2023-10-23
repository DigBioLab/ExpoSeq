
import pandas as pd


def cleaning(report, summed_clonefraction, max_num_reads):
    grouped = report.groupby('Experiment')

    # Step 2: Define a function to get the top rows adding up to 10%
    def get_top_10_percent(group):
        sorted_group = group.sort_values(by='cloneFraction', ascending=False)
        cumulative_sum = sorted_group['cloneFraction'].cumsum()
        return sorted_group[cumulative_sum <= summed_clonefraction][:max_num_reads]

    result_df = grouped.apply(get_top_10_percent)

    # Step 4: Reset index (if needed)
    result_df = result_df.reset_index(drop=True)
    report = result_df
    return report

