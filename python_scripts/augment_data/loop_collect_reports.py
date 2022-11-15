import pandas as pd
import re
from python_scripts.augment_data.structure_files import find_exp_name
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
                if experiment_name not in list_experiments:
                    local_intermediate_file.insert(0,
                                                   "Experiment",
                                                   experiment_name)
                    all_intermediate_files = pd.concat([all_intermediate_files, local_intermediate_file])
                    list_experiments.append(experiment_name)
                else:
                    pass
            else:
                pass
    return all_intermediate_files


def collect_nocluster_files(filenames): # with local_intermediate_report.shape[1] it only collect cdr3 files not uniquecdr. uniquecdr has only 7 columns while cdr3 has 9. ask chirs why
    all_intermediate_files = pd.DataFrame()
    each_instance = pd.DataFrame()
    list_experiments = []
    for file in filenames:
        try:
            local_intermediate_report = pd.read_table(file)
            experiment_name = find_exp_name(file)
            if local_intermediate_report.shape[1] > 8 and experiment_name not in list_experiments:
                local_intermediate_report.insert(0,
                                                 "Experiment",
                                                 experiment_name)
                instance = local_intermediate_report.iloc[:1]
                all_intermediate_files = pd.concat([all_intermediate_files, local_intermediate_report])
                each_instance = pd.concat([each_instance, instance])
                list_experiments.append(experiment_name)
            else:
                pass
        except:
            pass
    return all_intermediate_files, each_instance

