from argparse import ArgumentParser
import os


class TestArgs:
    def __init__(self, args):
        self.args = args
        
    def check_tsv_dir(self):
        assert os.path.isdir(self.args.tsv_dir) == True, "Path to directory with tsv files does not exist"

    def check_region(self):
        assert self.args.region in  ["targetSequences", "aaSeqFR1", "aaSeqCDR1", "aaSeqFR2", "aaSeqCDR2","aaSeqFR3", "aaSeqCDR3", "aaSeqFR4"]

    def check_length_threshold(self):
        assert self.args.length_threshold >=0, "Length threshold cannot be negative"
        
    def check_min_read_count(self):
        assert self.args.min_read_count >= 0, "Minimal read count cannot be negative"
        
    def check_remove_gaps(self):
        assert type(self.args.remove_gaps) == bool, "--remove_gaps must be boolean"
        
    def check_remove_errors(self):
        assert type(self.args.remove_errors) == bool, "--remove_errors must be boolean"
        
    def check_save_csv(self):
        dir_csv = os.path.dirname(self.args.save_csv) 
        assert os.path.isdir(dir_csv), f"The directory for {os.path.basename(self.args.csv_csv)} does not exist"
        
class ExpoSeqArgs:
    def __init__(self, **kwargs) -> None:
        self.parser = ArgumentParser(description = "Argument parser for plot function",
                                     allow_abbrev = True,
                                     **kwargs)
        self.chosen_tests = []
    def save_csv(self, flag = "--save_csv"):
        self.parser.add_argument(flag,
                                 help = "This must be parsed as path to tell where the output should be saved.",
                                 required = True)
        self.chosen_tests.append("check_save_csv")
        
    def sequencing_report(self, shortcut = "-r", flag = "--sequencing_report"):
        self.parser.add_argument(shortcut, 
                                 flag,
                                 help = "Insert the sequencing report (merged pd.dataframe of all samples) here.",
                                 required = True)
    
    def samples(self, shortcut = "-s", flag = "--samples"):
        self.parser.add_argument(shortcut,
                                 flag, 
                                 help = "Add a list of sample names here.",
                                 required = True,
                                 type = list)
        
    def region_of_interest(self, flag = "--region", default = "aaSeqCDR3"):
        self.parser.add_argument(flag, 
                                 help = "This must be given to indicate which region should be analysed.",
                                 required=True,
                                 type = str)
        self.chosen_tests.append("check_region")
        
    def pca_components(self, flag = "--pca_components", default_value = 50):
        self.parser.add_argument(flag, 
                                 help = "Add an integer for the number of pca_components you want to reduce the embedding to",
                                 type = int,
                                 default = default_value)
        
    def perplexity(self, flag = "--perplexity", default_value = 20 ):
        self.parser.add_argument(flag,
                                 help ="Add an integer for the perplexity value the t-SNE algorithm needs to arrange your points in a certain cluster size.",
                                 type = int,
                                 default = default_value)
    def iterations_tsne(self, flag = "--iterations_tsne", default_value = 1500):
        self.parser.add_argument(flag, 
                                 help = "Add an integer for the number of iterations the t-SNE algorithm should take to arrange the coordiantes for the points.",
                                 default = default_value,
                                 type=int)
    def batch_size(self, flag = "--batch_size", default_value = 600):
        self.parser.add_argument(flag, 
                                 help = "Add an integer for the batch size you want to take per sample which equals to the number of sequences with the top clone counts",
                                 type = int,
                                 default = default_value)
        
    def model_type(self, flag = "--model_type", default_value = 'Rostlab/prot_t5_xl_half_uniref50-enc'):
        self.parser.add_argument(flag,
                                 help = "Add a str as model type from the available models to embed your sequences in a high dimensional space.", 
                                 type = str,
                                 default= default_value)
    def strands(self, flag = "--show_strands", default_value = True):
        self.parser.add_argument(flag, 
                                 help = "Enter a boolean to show a batch of sequences in your plots as identifier. False does not show any.",
                                 type = bool,
                                 default = default_value)
        
    def tsv_dir(self, flag = "--tsv_dir"):
        self.parser.add_argument(flag,
                                 help = "Add a path to a directory with tsv files here to generate the sequencing report.",
                                 type = str,
                                 required= True)
        self.chosen_tests.append("check_tsv_dir")
    
    def length_threshold(self, flag = "--length_threshold", default_value = 5):
        self.parser.add_argument(flag, 
                                 help = "Enter the threshold for the minimal sequence length in the sequencing report here.",
                                 default = default_value,
                                 type = int
                                 )
        self.chosen_tests.append("check_length_threshold")
        
    def min_read_count_threshold(self, flag = "--min_read_count", default_value = 3):
        self.parser.add_argument(flag,
                                 help = "Add a value for the minimum required read count certain sequences should have to remain in the sequencing report.",
                                 default = default_value,
                                 type = int)
        self.chosen_tests.append("check_min_read_count")
        
    def remove_gaps(self, flag = "--remove_gaps", default_value = False):
        self.parser.add_argument(flag,
                                 help = "Here you can enter whether you would like to remove gaps in your sequencing data. That means that your target sequence should be divisible by 3.",
                                 default=default_value,
                                 type = bool)
        self.chosen_tests.append("check_remove_gaps")
        
    def remove_errors(self, flag = "--remove_errors", default_value = True):
        self.parser.add_argument(flag,
                                 help = "If set to True this will remove all sequences (rows) in the sequencing report which have a *.",
                                 default=default_value,
                                 type = bool)
        self.chosen_tests.append("check_remove_errors")
        
    
        