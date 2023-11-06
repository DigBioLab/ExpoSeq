from .modelling.nanonet import Use_NanoNet
from .modelling.create_fasta import Fasta
from modelling import multi_sequence_align
from .settings.change_settings import Settings
from tkinter.filedialog import askopenfilename
import subprocess
import sys

class ModelProtein:
    
    def __init__(self) -> None:
        Sets = Settings.check_dirs()
        if self.git_available() == False:
            raise Exception("Please install git to be able to use this module.")
        
    def model_nanonet(modeller = False):
        """
        :param modeller: If set to true modeller will be used for side chain reconstruction.
        :result: Output is a pdb file which is the predicted 3D structure from the linear sequence which is the input.
        """
        fasta_path = askopenfilename()
        is_fasta = Fasta.check_fasta(fasta_path)
        if is_fasta:
            Use_NanoNet(fasta_path, modeller)
        
    def create_fasta(sequence, sequence_name, fasta_name):
        Fasta.make_fasta(sequence, sequence_name, fasta_name)
        
    def multi_sequence_alignment(fasta_path, out_path):
        multi_sequence_align.MSA(fasta_path, out_path)
        
    
    @staticmethod
    def git_available():
        try:
            output = subprocess.check_output(['git', '--version'], stderr=subprocess.STDOUT, text=True)
            print(f"Git is available. Version: {output.strip()}")
            git_avail = True
        except subprocess.CalledProcessError as e:
            print("Git is not available.")
            git_avail = False
            return git_avail
        