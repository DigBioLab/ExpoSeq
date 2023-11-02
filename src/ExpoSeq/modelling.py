
from .modelling.nanonet import Use_NanoNet
from .modelling.create_fasta import make_fasta
from modelling import multi_sequence_align
from .settings.change_settings import Settings
from tkinter.filedialog import askopenfilename



class ModelProtein:
    
    def __init__(self) -> None:
        Sets = Settings.check_dirs()

        
    def model_nanonet(modeller = False):
        """
        :param modeller: If set to true modeller will be used for side chain reconstruction.
        :result: Output is a pdb file which is the predicted 3D structure from the linear sequence which is the input.
        """
        fasta_path = askopenfilename()
        Use_NanoNet(fasta_path, modeller)
        
    def create_fasta(sequence, sequence_name, fasta_name):
        make_fasta(sequence, sequence_name, fasta_name)
        
    def multi_sequence_alignment(fasta_path, out_path):
        multi_sequence_align.MSA(fasta_path, out_path)
        
    
    