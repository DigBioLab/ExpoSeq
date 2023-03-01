import os
def remove_files_and_directories(directory):
    for filename in os.listdir(directory):
        file_path = os.path.join(directory, filename)
        try:
            if os.path.isfile(file_path):
                os.remove(file_path)
            elif os.path.isdir(file_path):
                os.rmdir(file_path)
        except Exception as e:
            print(f"Failed to remove {file_path}: {e}")
def original_settings():
    module_dir = os.path.abspath("")
    experiments_dir = os.path.join(module_dir,
                                   "expoSeq",
                                   "my_experiments")
    remove_files_and_directories(experiments_dir)
    common_vars_path = os.path.join(module_dir,
                                    "settings",
                                    "global_vars.txt")
    with open(common_vars_path, "r") as f:
        global_params = f.read()
    global_params["mixcr_path"] = ""
    global_params['last_experiment'] = ""




