import shutil
import os

def remove_directories(directory_path):
    # Iterate over all items in the directory
    for item in os.listdir(directory_path):
        shutil.rmtree(item)

def original_settings():
    module_dir = os.path.abspath("expoSeq")
    experiments_dir = os.path.join(module_dir,
                                   "my_experiments")
    remove_directories(experiments_dir)
    common_vars_path = os.path.join(module_dir,
                                    "settings",
                                    "global_vars.txt")
    with open(common_vars_path, "r") as f:
        global_params = f.read()
    global_params["mixcr_path"] = ""
    global_params['last_experiment'] = ""





