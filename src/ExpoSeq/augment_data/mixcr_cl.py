from ..augment_data.create_mixcr_commands import CreateCommand
from ..augment_data.collect_fastqs import CollectFastq
import subprocess
import os
import pandas as pd
import pickle
import argparse
from ..settings.check_dirs import check_dirs
        
def mixcr_(path_to_mixcr, experiment_name,  path_to_forward, path_to_backward = None, method = "ampliseq-tcrb-plus-full-length", threads = 1, java_heap_size = None):
    print(path_to_mixcr)
    print(path_to_forward)
    print(path_to_backward)
    print(method)
    print(threads)
    print(java_heap_size)
    module_dir = os.path.abspath("")
    assert os.path.isfile(path_to_mixcr), "File path to mixcr.jar file is incorrect."
    assert path_to_forward != path_to_backward, "Directory to forward and backward files needs to be different"
    assert os.path.isdir(path_to_forward), "The given directory for the forward files does not exist."
    assert subprocess.run(["java", "-jar",path_to_mixcr,"--version"]), "The given path to mixcr is not correct. Or you have not installed mixcr correctly."
    assert type(experiment_name) == str, "Please enter a string as experiment name"
  #  assert os.path.isfile(module_dir, "my_experiments", experiment_name) == False, "The experiment name already exists, please enter another one."
    check_dirs(module_dir, experiment_name)
    if java_heap_size == None:
        java_heap_size = 1000
    if path_to_backward == None:
        paired_end_sequencing = False
    else:
        paired_end_sequencing = True

    Collect = CollectFastq(paired_end_sequencing)
    
    Collect.module_dir = module_dir

    Collect.get_files(cmd = True,
                      path_to_files = path_to_forward,
                      path_to_backward = path_to_backward)
    files = Collect.paired
    print(files)
    if paired_end_sequencing:
        print("Found pairs for the forward and reverse files")
        Collect.call_pairs()

            
    sequencing_report = pd.DataFrame([])
    for filename in files:
        Commands = CreateCommand(module_dir, path_to_mixcr, paired_end_sequencing, experiment_name, filename)
        print("go")
        subprocess.run(Commands.prepare_align(method))
        skip_sample = False
        
        align_quality, percent_quality = Commands.read_aligned_reads()
        if align_quality < 1000 and align_quality >10:
            print(f"Only {align_quality} reads could be aligned ({percent_quality}).")
        if align_quality < 10:
            print("Too less reads were aligned. Please choose the correct mixcr method.")
            skip_sample = True
        if skip_sample == False:
            subprocess.run(Commands.prepare_assembly())
            assembly_qual, realtive_assembly_quality = Commands.read_clones()
            if assembly_qual < 1000 and realtive_assembly_quality > 10:
                print(f"Only {assembly_qual} reads could be assembled ({realtive_assembly_quality}).")
            if assembly_qual < 10:
                print("Too less reads were assembled. Please choose the correct ")
                skip_sample = True
            if skip_sample == False:
                subprocess.run(Commands.prepare_clones())
                basename = Commands.basename

                for filename in os.listdir(os.path.join(module_dir, "temp")):
                    if filename.endswith(".tsv"):
                        clones_sample = pd.read_table(os.path.join(module_dir, "temp", filename))
                clones_sample["Experiment"] = basename
                sequencing_report = pd.concat([sequencing_report, clones_sample])

                for file in os.listdir(os.path.join(module_dir, "temp")):
                    os.remove(os.path.join(module_dir,
                                        "temp",
                                        file))
     #   except:
      #      print(f"File {filename} could not be processed. Please check if you have chosen the correct mixcr method or verify that the sequencing was successful.")

    sequencing_report.to_csv(os.path.join(module_dir, experiment_name, 'sequencing_report.csv'))
    if sequencing_report.shape[0] == 0:
        raise Exception(
            "The processing failed for all of your files. Most likely because you chose the wrong mixcr method.")

    unique_experiments = sequencing_report["Experiment"].unique()
    experiment_dic = {item: item for item in list(unique_experiments)}
    exp_names_dir = os.path.join(module_dir,
                                    "my_experiments",
                                    experiment_name,
                                    "experiment_names.pickle")
    with open(exp_names_dir, "wb") as f:
        pickle.dump(experiment_dic, f)
    print(
    "The output of your processing is a sequencing_report which contains the sequences and their sample names of all fastq files. Further, a folder with the alignment reports for each fastq file was created. You can have a look under ~/my_experiments/YOUR_EXPERIMENT_NAME")



parser = argparse.ArgumentParser(description='Processing')
parser.add_argument('path_to_mixcr', type=str, help="Path to the mixcr jar file")
parser.add_argument('experiment_name', type=str, help="Name of the experiment")
parser.add_argument('--path_to_forward',type = str, help='Directory to the fastq files with the forward reads')
parser.add_argument('--path_to_backward',type = str,default=None, help='Directory to the fastq files with the backward reads')
parser.add_argument('--threads', default=1, type=int, help='Number of threads to use')
parser.add_argument("--method", default="milab-human-tcr-dna-multiplex-cdr3", type=str, help="Method to use for the alignment")
parser.add_argument("--java_heap_size", default="1000", type=int, help="Memory allocation for the processing in MB")

args = parser.parse_args()
path_to_mixcr = args.path_to_mixcr
experiment_name = args.experiment_name
path_to_forward = args.path_to_forward
path_to_backward = args.path_to_backward
threads = args.threads
method = args.method
java_heap_size = args.java_heap_size

# Call the function
mixcr_(path_to_mixcr, experiment_name, path_to_forward, path_to_backward,  threads, method, java_heap_size)
#PlotManager(experiment = experiment_name)