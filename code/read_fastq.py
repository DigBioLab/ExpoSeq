from Bio import SeqIO
from io import StringIO
import skbio # needs pip install pingouin
import izip
from skbio.parse.sequences import FastaIterator, FastqIterator
import itertools

file = gzip.open(r"C:\Users\nils\Desktop\Tropical Pharmaceutical Lab\Project\NG-28545_Library_1_F5_2.fastq")

fq = SeqIO.parse(r"C:\Users\nils\Desktop\Tropical Pharmaceutical Lab\Project\Library_1_F1_2_UniqueCDR3_Exp_UniqueCDR3", "txt")
fq = SeqIO.parse(r"C:\Users\nils\Desktop\Tropical Pharmaceutical Lab\Project\NG-28545_Library_1_F5_2.fastq", "fastq")

for read in fq:
    title_line = fq.records.gi_frame.f_locals["title_line"]
    seq_string = fq.records.gi_frame.f_locals["seq_string"]
    read_id = fq.records.gi_frame.f_locals["id"]
    read_name = fq.records.gi_frame.f_locals["name"]
    read_descr = fq.records.gi_frame.f_locals["descr"]

    read_quality = fq.records.gi_frame.f_locals["quality_string"]


    #read.letter_annotations is phred_quality but it is not the quality


