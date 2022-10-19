from Bio import SeqIO
import gzip
from os.path import normpath
from tkinter.filedialog import askopenfilename

#file = gzip.open(r"C:\Users\nils\Desktop\Tropical Pharmaceutical Lab\Project\NG-28545_Library_1_F5_2.fastq")

#fq = SeqIO.parse(r"C:\Users\nils\Desktop\Tropical Pharmaceutical Lab\Project\Library_1_F1_2_UniqueCDR3_Exp_UniqueCDR3", "txt")
#fq = SeqIO.parse(r"C:\Users\nils\Desktop\Tropical Pharmaceutical Lab\Project\NG-28545_Library_1_F5_2.fastq", "fastq")

def read_fastq():
    sequences = []
    quality = []
    filename = askopenfilename()
    filename = normpath(filename)
    for seq_record in SeqIO.parse(filename, 'fastq'):
        local_sec = str(seq_record.seq)
        sequences.append(local_sec)
        quality.append(seq_record.letter_annotations['phred_quality'])
    return sequences, quality
#for read in fq:
 #   title_line = fq.records.gi_frame.f_locals["title_line"]
  #  seq_string = fq.records.gi_frame.f_locals["seq_string"]
   # read_id = fq.records.gi_frame.f_locals["id"]
    #read_name = fq.records.gi_frame.f_locals["name"]
    #read_descr = fq.records.gi_frame.f_locals["descr"]

    #read_quality = fq.records.gi_frame.f_locals["quality_string"]


    #read.letter_annotations is phred_quality but it is not the quality


