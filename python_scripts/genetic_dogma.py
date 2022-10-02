from Bio.Seq import Seq

def genetic_dogma(dna):
    dna = Seq(dna)
    m_rna = dna.transcribe()
    protein = m_rna.translate()
    protein = str(protein)
    return protein


