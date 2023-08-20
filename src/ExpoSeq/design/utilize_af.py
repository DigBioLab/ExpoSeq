import subprocess
import os
from anno_sequence import SequenceManipulator

class Docker:
    def __init__(self):
        self.installed = False
        self.runs = False
    def is_docker_installed(self):
        try:
            subprocess.run(["docker", "--version"], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            self.installed = True

        except subprocess.CalledProcessError:
            self.installed = False

    def is_docker_daemon_running(self):
        try:
            subprocess.run(["docker", "info"], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            self.runs = True
        except subprocess.CalledProcessError:
            self.runs = False
            
    def install_message(self):
        if self.installed == False:
            print("Docker is not installed. Please install docker from: " + "https://docs.docker.com/get-docker/")
    
    def run_message(self):
        if self.runs == False:
            print("Docker is not running. Please check if you have created an account and that you are logged in.")

    def is_image_built(self, image_name ="alphafold"):
        cmd = ["docker", "images", "-q", image_name]
        result = subprocess.run(cmd, stdout=subprocess.PIPE, text=True)
        return bool(result.stdout.strip())



def run_alphafold(fasta_path, output_dir):
    cmd = [
        'docker', 'run',
        '-v', f"{fasta_path}:/data",
        '-v', f"{output_dir}:/output",
        'alphafold',
        '--fasta_paths=/data/my_sequence.fasta',
        '--max_template_date=YYYY-MM-DD',
        '--output_dir=/output'
    ]
    subprocess.run(cmd)



class DesignFasta():
    def __init__(self, template_sequence):
        self.reference = template_sequence
        self.new_query = template_sequence
        self.sequence_name = "my_sequence"
    
    def define_sequence_name(self):
        self.sequence_name = input("Please enter a name for your sequence: ")

    def manipulate_seq(self):
        app = SequenceManipulator(self.reference)
        app.run()
        self.new_query = app.sequence
        
    def write_fasta(self, module_dir):
        if not os.path.isdir("myfastasequences"):
            os.mkdir("myfastasequences")
        self.define_sequence_name()
        full_path = os.path.abspath(os.path.join(module_dir,"temp", self.sequence_name + ".fasta"))
        print(full_path)
        with open(full_path, "w") as f:
            f.write(">" + self.sequence_name + "\n")
            f.write(self.new_query.get())


def folding(template_seq, output_dir):
    module_dir = r"C:\Users\nilsh\my_projects\ExpoSeq"
    My_fasta = DesignFasta(template_seq)
    My_fasta.manipulate_seq()
    My_fasta.write_fasta(module_dir)
    filename = My_fasta.sequence_name
    My_docker = Docker()
    os.chdir(os.path.join(module_dir, "alphafold"))
    if My_docker.is_image_built():
        print(filename)
      #  print()
        subprocess.run(["python3",
                        r"C:\Users\nilsh\my_projects\ExpoSeq\alphafold\docker\run_docker.py",
                        "--fasta_paths=" + os.path.join(module_dir, "temp", filename + ".fasta"),
                        "--model_preset=monomer",
                        "--db_preset=reduced_dbs",
                        "--output_dir="+output_dir,
                        ])
     #   run_alphafold(os.path.join(module_dir, "temp", filename + ".fasta"), output_dir)
        print("Prediction_succesful")
    else:
        os.chdir(module_dir)
        My_docker.is_docker_installed()
        My_docker.install_message()
        My_docker.is_docker_daemon_running()
        My_docker.run_message()
        if not os.path.exists("alphafold"):
            subprocess.run(["git",
                    "clone",
                    "https://github.com/deepmind/alphafold.git"])
        module_dir = os.getcwd()
        if My_docker.runs == True and My_docker.installed == True:
            os.chdir(os.path.join(module_dir, "alphafold"))
            if not My_docker.is_image_built():
                subprocess.run(["docker",
                                "build",
                                "-f",
                                "docker/Dockerfile",
                                "-t",
                                "alphafold",
                                "."])  # "." indicates current directory as build context
            
            run_alphafold(os.path.join(module_dir, "temp", filename + ".fasta"), output_dir)
    os.chdir(module_dir)
        
folding("QMQLVQSGAEVKKPGSSVKVSCKASGGTFSSYAISWVRQAPGQGLEWMGGIIPIFGTANYAQKFQGRVTITADESTSTAYMELSSLRSEDTAVYYCARSSGYYGMDVWGQGTLVTVSS", r"C:\Users\nilsh\my_projects\ExpoSeq\my_experiments\test_true\alphafold_results")

