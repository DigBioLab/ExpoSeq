import os
import subprocess
import pandas as pd
from .full_sequence_finder import FullSequence

class ExtractFasta:
    def __init__(self, merged_fasta_file:str, save_dir) :
        assert merged_fasta_file.endswith('.fasta') or merged_fasta_file.endswith('.fa')
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
        with open(self.merge_fasta_file, 'w') as oh:
            first_line = True
            for f in os.listdir(dir):
                if f.endswith('.fasta') or f.endswith('.fa') or f.endswith('.fas'):
                    with open(os.path.join(dir, f)) as fh:
                        lines = fh.readlines()
                        for i, line in enumerate(lines):
                            if line.startswith('>'):
                                if first_line:
                                    oh.write(line)
                                    first_line = False
                                else:
                                    oh.write('\n' + line)
                            else:
                                for j in range(0, len(line), 80):
                                    oh.write(line[j:j+80])
                                    if j + 80 < len(line) or i != len(lines) - 1:
                                        oh.write('\n')
                
    def align_mixcr(self, mixcr_path:str, preset_name = "exom-full-length", species = ""):
        """_summary_

        Args:
            mixcr_path (str): _description_
            preset_name (str, optional): _description_. Defaults to "milab-human-bcr-multiplex-full-length".
            species (str, optional): For available species check out:https://github.com/repseqio/library-imgt/releases . Defaults to "".
        """
        command = ["java",
                   "-jar", 
                   mixcr_path, 
                   "analyze",
                   preset_name,
                   "--keep-non-CDR3-alignments",
                   self.merge_fasta_file,
                   os.path.join(self.save_dir, "result_mixcr_fasta"),
                   "-f",
                   "--no-json-reports"
                #   f"{str(os.path.join(self.save_dir,'sample_table.tsv'))}"
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
        sequence = ''
        with open(fasta_file, 'r') as file:
            for line in file:
                if line.startswith('>'):
                    headers.append(line.strip())
                    if sequence != '':
                        sequences.append(sequence)
                        sequence = ''
                else:
                    sequence += line.strip()
            sequences.append(sequence)  # for the last sequence
        df = pd.DataFrame({'Header': headers, 'targetSequences': sequences})
        return df
    @staticmethod
    def sequence_in_other(row, df_long):
        mask = df_long['targetSequences'].str.contains(row['nSeqCDR3'])
        if mask.any():
            return df_long.loc[mask, 'Header'].iloc[0]
        else:
            return None

                                   
    def merge_and_read_seq(self, ):
        all_seqs = self.get_headers_and_sequences_from_fasta(self.merge_fasta_file)
        tsv_filename = os.path.join(self.save_dir, "result_mixcr_fasta" + ".clones_IGH" + ".tsv")
        tsv = pd.read_table(tsv_filename)
        tsv['Header'] = tsv.apply(self.sequence_in_other, axis=1, df_long=all_seqs)  
    #    tsv_merge = tsv_merge[~tsv_merge.apply(lambda row: row.astype(str).str.contains('region not covered').any(), axis=1)]
        aaSeq_columns = [col for col in tsv.columns if col.startswith('aaSeq')]
        aaSeq_columns = [col.replace('aaSeq', '') for col in tsv.columns if col.startswith('aaSeq')]
        row_indexes = FullSequence(avail_regions=aaSeq_columns).find_connecting_seq()
        avail_regions = ["".join("aaSeq" + region) for region in row_indexes]
        tsv['merged_regions'] = tsv[avail_regions].apply(lambda row: ''.join(row.values.astype(str)), axis=1)
        tsv = tsv[["Header", "merged_regions"]]
        tsv = tsv[["Header", "merged_regions"]].dropna(subset=['Header'])
        tsv = tsv.rename(columns={'merged_regions': 'aaSeqtargetSequences'})
        cols = list(tsv.columns)
        cols.insert(0, cols.pop(cols.index('aaSeqtargetSequences')))
        tsv = tsv.loc[:, cols]
        tsv.to_csv(tsv_filename.replace(".tsv", "_merged.csv"))
        

Extract = ExtractFasta("my_test_merged.fasta", save_dir = r"C:\Users\nilsh\my_projects\ExpoSeq\my_experiments\max_both")
#Extract.merge_fasta_files(r"C:\Users\nilsh\OneDrive\Desktop\master_thesis\binding_data\RE_ binding data from Discovery of broadly neutralizing nanobodies using designed consensus antigens\Nanobody sequences\Sequences")
Extract.align_mixcr(r"C:\Users\nilsh\OneDrive\Desktop\DTU\NGS_pipeline\mixcr-4.2.0\mixcr.jar", species= "lamaGlama")
Extract.merge_and_read_seq()