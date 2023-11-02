from Bio.Align.Applications import ClustalwCommandline
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from ExpoSeq.design.msa_design import ExtendedSequenceManipulator
import os
import pandas as pd
import tkinter as tk
from tkinter import filedialog

class MSA:
    def __init__(self, fasta_path, out_path) -> None:
        clustalw_path = self.get_clustalw_path()
        self.run_msa(clustalw_path, fasta_path, out_path)
    
    def get_clustalw_path(self, Settings):
        global_params = Settings.read_global_vars()
        clustalw_path = global_params["clustalw_path"]
        if clustalw_path == '':
            print("Please choose the path to the clustalw exe with the filechooser. If you don't have it, you can download it from http://www.clustal.org/download/current/")
            clustalw_path = filedialog.askopenfilename()
            Settings.change_global_vars("clustalw_path", clustalw_path)

        return clustalw_path
    
    def run_msa(clustalw_path, seqs_path, outfile):         
        clustalw_cline = ClustalwCommandline(clustalw_path, infile=seqs_path, outfile=outfile)
        stdout, stderr = clustalw_cline()