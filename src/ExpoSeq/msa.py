

class RUNMSA:
    def get_clustalw_path(self):
        global_params = self.Settings.read_global_vars()
        clustalw_path = global_params["clustalw_path"]
        try:
            import tkinter as tk
            from tkinter import filedialog
            if clustalw_path == '':
                print("Please choose the path to the clustalw exe with the filechooser. If you don't have it, you can download it from http://www.clustal.org/download/current/")
                clustalw_path = filedialog.askopenfilename()
                self.Settings.change_global_vars("clustalw_path", clustalw_path)
        except:
            print("Tk not installed. Please install it to use this functionality and visualize the MSA")
        return clustalw_path
    
    def run_MSA(self, batch_size, samples):
        path_clustalw = self.get_clustalw_path()
        if not path_clustalw == '' and os.path.isfile(path_clustalw):
            try:
                os.remove(self.outfile)
            except:
                pass
            sequences = self.get_top_sequences(batch_size, samples)
            self.create_fasta(sequences)
            clustalw_cline = ClustalwCommandline(path_clustalw, infile=self.seqs_path, outfile=self.outfile)
            stdout, stderr = clustalw_cline()
            os.remove(self.outfile)