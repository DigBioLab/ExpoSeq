from Bio import SeqIO
import re



class Fasta:
    def make_fasta(sequence, sequence_name, fasta_name):

        fasta_record = SeqIO.SeqRecord(
            seq=sequence,
            id=sequence_name,
            description="",
        )

        with open(f"{fasta_name}.fasta", "a") as f:
            SeqIO.write(fasta_record, f, "fasta")
        
    
        
    def check_fasta(self, fasta_file):
        """Checks if a FASTA file contains amino acid sequences.

        Args:
            fasta_file: The path to the FASTA file.

        Returns:
            True if the FASTA file contains amino acid sequences, False otherwise.
        """


        fasta_records = SeqIO.parse(fasta_file, "fasta")

        # Check the sequence characters
        for fasta_record in fasta_records:
            
            sequence = fasta_record.seq
            print(sequence)
            if not all(x in ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L",
                                "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]
                        for x in sequence):
                return False


        try:
            # Read the FASTA file
            fasta_records = SeqIO.parse(fasta_file, "fasta")

            # Check all of the headers
            for fasta_record in fasta_records:
                header = fasta_record.id
                print(header)
                if not self.is_valid_fasta_header(header):
                    return False

        except Exception as e:
                # If there is an exception, the FASTA file is not valid
            return False

        return True
    
    @staticmethod
    def is_valid_fasta_header(header):
        """Checks if a FASTA header is valid.

        Args:
            header: The FASTA header to check.

        Returns:
            True if the FASTA header is valid, False otherwise.
        """

        # Check if the header starts with a greater than sign
   #     if not header.startswith(">"):
   #         return False

        # Check if the header contains any whitespace characters
        if re.search(r"\s+", header) is not None:
            return False

        # Check if the header contains any invalid characters
        if re.search(r"[^\w\-]", header) is not None:
            return False

        return True
    
    
    
