

import subprocess
from modeller import *
import os




class Use_NanoNet:
    def __init__(self,module_dir, fasta_path, modeller = False) -> None:
        self.git_available()
        self.check_tensorflow()
        self.prepare_NanoNet()
        self.base_commands += [fasta_path]
        avail_modeller = self.add_modeller()
        if modeller:
            command_modeller = ["-m True"]
            self.base_commands += command_modeller
        if os.path.isfile(os.path.join(module_dir, "NanoNet","NanoNet.py")):
            print("Please cite:")
            self.print_citation()
            subprocess.run(self.base_commands)
        else:
            print("Cannot access NanoNet repo.")
        
    def git_available(self):
        try:
            output = subprocess.check_output(['git', '--version'], stderr=subprocess.STDOUT, text=True)
            print(f"Git is available. Version: {output.strip()}")
            self.clone_repo()
        except subprocess.CalledProcessError as e:
            print("Git is not available.")
            return 
        
    def print_citation():
        dir = os.path.dirname(__file__)
        with open(os.path.join(dir, "citations", "nanonet.txt")) as f:
            text = f.read()

            # Print the text
            print(text)

            # Close the file
            f.close()
    
    @staticmethod
    def check_tensorflow():
        try:
            import tensorflow as tf
            print(f"TensorFlow version: {tf.__version__}")
        except ImportError:
            print("TensorFlow is not installed.")      
            subprocess.run(["pip install tensorflow"])

    def clone_repo():
        subprocess.run(["git", "clone", "https://github.com/dina-lab3D/NanoNet.git"])
    
    def prepare_NanoNet(self, fasta_file):
        self.base_commands = ["python", 
                        "NanoNet/NanoNet.py",
                        fasta_file,
                        ]
        
    def add_modeller(self, ):
        try:
            env = environ()
            aln = alignment(env)
            avail_modeller = True
        except NameError:
            print("Modeller is not installed. You can continue without.")
            avail_modeller = False
            
    
