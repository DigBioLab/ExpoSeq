from Bio.Align.Applications import ClustalwCommandline
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from ExpoSeq.design.msa_design import ExtendedSequenceManipulator
import os
import pandas as pd



class MSA:
    def __init__(self, region_of_interest, SequencingReport,Settings, module_dir):
        self.region_of_interest= region_of_interest
        self.SeqReport = SequencingReport()
        self.Settings = Settings()
        self.seqs_path = os.path.join(module_dir, "temp", "sequences.fasta")
        self.outfile = os.path.join(module_dir, "temp", "aligned_clustal.fasta")
        
    
    def create_fasta(self, sequences):
        records = [SeqRecord(Seq(row["aaSeq_all"]), id=f"seq{i+1}", description='') for i, row in sequences.iterrows()]
        SeqIO.write(records, self.seqs_path, "fasta")
        
    
    def get_top_sequences(self, batch_size, samples):
        self.SeqReport.filter_longest_sequence()
        sequencing_report = self.SeqReport.sequencing_report
        sequencing_report = sequencing_report[sequencing_report['Experiment'].isin(samples)]
        sequencing_report = sequencing_report[~sequencing_report['aaSeq_all'].str.contains('\*')]
        sequencing_report = sequencing_report.groupby("Experiment").head(batch_size)
        sequences = sequencing_report["aaSeq_all"]
        sequences = pd.DataFrame(sequences.values, columns = ["aaSeq_all"])
        return sequences

    def get_clustalw_path(self):
        clustalw_path = self.Settings.global_params["clustalw_path"]
        try:
            import tkinter as tk
            from tkinter import filedialog
            if clustalw_path == '':
                print("Please enter the path to clustalw exe. If you don't have it, you can download it from http://www.clustal.org/download/current/")
                clustalw_path = filedialog.askopenfilename()
                self.Settings.change_global_vars("clustalw_path", clustalw_path)
        except:
            print("Tk not installed. Please install it to use this functionality and visualize the MSA")

    def run_MSA(self, batch_size, samples):
        path_clustalw = self.get_clustalw_path()
        if not path_clustalw == '' and os.path.isfile(path_clustalw):
            os.remove(self.outfile)
            sequences = self.get_top_sequences(batch_size, samples)
            self.create_fasta(sequences)
            clustalw_cline = ClustalwCommandline(path_clustalw, infile=self.seqs_path, outfile=self.outfile)
            stdout, stderr = clustalw_cline()
            example_seq = "YOUR_SEQUENCE_HERE"
            app = ExtendedSequenceManipulator(example_seq, self.outfile)
            app.run()
            os.remove(self.outfile)