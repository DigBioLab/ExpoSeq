import os
import subprocess
import pandas as pd
from ExpoSeq.settings.full_sequence_finder import FullSequence
from difflib import SequenceMatcher
import numpy as np

class ExtractFasta:
    def __init__(self, merged_fasta_file: str, save_dir):
        assert merged_fasta_file.endswith(".fasta") or merged_fasta_file.endswith(".fa")
        self.merge_fasta_file = merged_fasta_file
        self.save_dir = self.check_create_dir(save_dir)

    @staticmethod
    def check_create_dir(dir):
        save_dir = os.path.join(dir, "binder_import")
        if not os.path.exists(save_dir):
            os.makedirs(os.path.join(dir, "binder_import"))
        else:
            pass

        return save_dir

    def merge_fasta_files(self, dir):
        with open(self.merge_fasta_file, "w") as oh:
            first_line = True
            for f in os.listdir(dir):
                if f.endswith(".fasta") or f.endswith(".fa") or f.endswith(".fas"):
                    with open(os.path.join(dir, f)) as fh:
                        lines = fh.readlines()
                        for i, line in enumerate(lines):
                            if line.startswith(">"):
                                if first_line:
                                    oh.write(line)
                                    first_line = False
                                else:
                                    oh.write("\n" + line)
                            else:

                                oh.write(line)


    def align_mixcr(self, mixcr_path: str, preset_name="exom-full-length", species=""):
        """_summary_

        Args:
            mixcr_path (str): _description_
            preset_name (str, optional): _description_. Defaults to "milab-human-bcr-multiplex-full-length".
            species (str, optional): For available species check out:https://github.com/repseqio/library-imgt/releases . Defaults to "".
        """
        command = [
            "java",
            "-jar",
            mixcr_path,
            "analyze",
            preset_name,
            self.merge_fasta_file,
            os.path.join(self.save_dir, "result_mixcr_fasta"),
            "-f",
            "--no-json-reports",

        ]
        if species == "":
            pass
        else:
            command.append("--species")
            command.append(f"{species}")
        subprocess.run(command)

    @staticmethod
    def get_headers_and_sequences_from_fasta(fasta_file):
        headers = []
        sequences = []
        sequence = ""
        with open(fasta_file, "r") as file:
            for line in file:
                if line.startswith(">"):
                    headers.append(line.strip())
                    if sequence != "":
                        sequences.append(sequence)
                        sequence = ""
                else:
                    sequence += line.strip()
            sequences.append(sequence)  # for the last sequence
            
        df = pd.DataFrame({"Header": headers, "targetSequences": sequences})
        return df
    @staticmethod
    def sequence_in_other(row, df_long, supply_possible_cols = None):
        """_summary_

        Args:
            row (_type_): row from iteration from tsv from mixcr output.
            df_long (_type_): Sequences and headers from original fasta file

        Returns:
            _type_: _description_
        """
        if supply_possible_cols is not None:
            possible_cols = supply_possible_cols
        # Find the row in df_long where the 'nSeqtargetSequences' value is the same as in the current row
        for i, seq in enumerate(df_long["targetSequences"]):
            if row["nSeqCDR3"] in seq:
                matching_row = df_long.iloc[i, :]
                possible_header = matching_row.loc['Header']
                if supply_possible_cols is not None:
                    a = possible_header.split("_")
                    b = "_".join(a[:3])

                    modified_header = b[1:] + "_M"
                    if modified_header in possible_cols:
                        return modified_header
                    else:
                        pass
        # If no matching row was found, return NaN
        return np.nan


    @staticmethod   
    def translate_sequence(amino_acid_sequence):
        # Define the codon table
        codon_table = {
            "F": "TTT",
            "L": "TTA",
            "I": "ATT",
            "M": "ATG",
            "V": "GTT",
            "S": "TCT",
            "P": "CCT",
            "T": "ACT",
            "A": "GCT",
            "Y": "TAT",
            "H": "CAT",
            "Q": "CAA",
            "N": "AAT",
            "K": "AAA",
            "D": "GAT",
            "E": "GAA",
            "C": "TGT",
            "W": "TGG",
            "R": "CGT",
            "G": "GGT",
            "*": "TAA",
            "_": "",
            
        }

        # Translate the amino acid sequence to a nucleotide sequence
        nucleotide_sequence = "".join([codon_table[aa] for aa in amino_acid_sequence])

        return nucleotide_sequence

        
    def merge_and_read_seq(
        self,
        tsv_filename = None, 
        naming_col = None
    ):
        all_seqs = self.get_headers_and_sequences_from_fasta(self.merge_fasta_file)
        if tsv_filename is None:
            tsv_filename = os.path.join(
                self.save_dir, "result_mixcr_fasta" + ".clones_IGH" + ".tsv"
            )
        else:
            tsv_filename = tsv_filename
        tsv = pd.read_table(tsv_filename)

        #    tsv_merge = tsv_merge[~tsv_merge.apply(lambda row: row.astype(str).str.contains('region not covered').any(), axis=1)]
        aaSeq_columns = [col for col in tsv.columns if col.startswith("aaSeq")]
        aaSeq_columns = [
            col.replace("aaSeq", "") for col in tsv.columns if col.startswith("aaSeq")
        ]
        row_indexes = FullSequence(avail_regions=aaSeq_columns).find_connecting_seq()
        avail_regions = ["".join("aaSeq" + region) for region in row_indexes]
        tsv["merged_regions"] = tsv[avail_regions].apply(
            lambda row: "".join(row.values.astype(str)), axis=1
        )
        tsv = tsv[~tsv['merged_regions'].str.contains('region_not_covered')]
        tsv["nSeqtargetSequences"] = tsv["merged_regions"].apply(self.translate_sequence)
        tsv['Header'] = tsv.apply(lambda row: self.sequence_in_other(row, df_long=all_seqs, supply_possible_cols=naming_col), axis=1)
        tsv = tsv[["Header", "merged_regions"]]
        tsv = tsv[["Header", "merged_regions"]].dropna(subset=["Header"])
        tsv = tsv.rename(columns={"merged_regions": "aaSeqtargetSequences"})
        cols = list(tsv.columns)
        cols.insert(0, cols.pop(cols.index("aaSeqtargetSequences")))
        print(tsv)
        tsv = tsv.loc[:, cols]
        
        
        tsv.to_csv(tsv_filename.replace(".tsv", "_merged.csv"))

