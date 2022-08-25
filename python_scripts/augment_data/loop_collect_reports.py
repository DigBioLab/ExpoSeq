import pandas as pd
import re
def collect_intermediate_files(grouped_filenames, local_pattern_more_digits):
    all_intermediate_files = pd.DataFrame()
    #max_files = len(grouped_filenames)
    #limit = int(input(f"how many do you want to plot? Max No. is {max_files}"))
    list_experiments = []
    for cluster in grouped_filenames:
        for file in cluster:
            local_intermediate_file = pd.read_table(file)
            if len(local_intermediate_file.columns) > 2:
                experiment_name = re.search(local_pattern_more_digits,
                                            file).group()
                local_intermediate_file.insert(0,
                                                "Experiment",
                                                experiment_name)

                if experiment_name not in list_experiments:
                    all_intermediate_files = pd.concat([all_intermediate_files, local_intermediate_file])
                    list_experiments.append(experiment_name)
                else:
                    pass
            else:
                pass
    return all_intermediate_files
