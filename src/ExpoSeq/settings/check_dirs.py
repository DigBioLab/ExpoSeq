import os


def check_dirs(module_dir, experiment):
    if not os.path.isdir(os.path.join(module_dir, "my_experiments")):
        os.mkdir(os.path.join(module_dir, "my_experiments"))
    if not os.path.isdir(os.path.join(module_dir, "my_experiments", experiment)):
        os.mkdir(os.path.join(module_dir, "my_experiments", experiment))      
    if not os.path.isdir(os.path.join(module_dir, "my_experiments", experiment, "alignment_reports")):
        os.mkdir(os.path.join(module_dir, "my_experiments", experiment, "alignment_reports"))
    if not os.path.isdir(os.path.join(module_dir, "temp")):
        os.mkdir(os.path.join(module_dir, "temp"))   
    if not os.path.isdir(os.path.join(module_dir, "my_experiments", experiment, "assembly_reports")):
        os.mkdir(os.path.join(module_dir, "my_experiments", experiment, "assembly_reports"))
    if not os.path.isdir(os.path.join(module_dir, "my_experiments", experiment, "clones_result")):
        os.mkdir(os.path.join(module_dir, "my_experiments", experiment, "clones_result"))
    if not os.path.isdir(os.path.join(module_dir,
                                     "my_experiments",
                                     experiment,
                                     "tables_mixcr")):
        os.mkdir(os.path.join(module_dir, "my_experiments", experiment, "tables_mixcr"))