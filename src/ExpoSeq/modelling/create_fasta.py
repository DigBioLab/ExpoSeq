from Bio import SeqIO

def make_fasta(sequence, sequence_name, fasta_name):

    fasta_record = SeqIO.SeqRecord(
        seq=sequence,
        id=sequence_name,
        description="",
    )

    with open(f"{fasta_name}.fasta", "a") as f:
        SeqIO.write(fasta_record, f, "fasta")