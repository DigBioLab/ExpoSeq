from argparse import ArgumentParser
import os
import pandas as pd

class TestArgs:
    def __init__(self, args):
        self.args = args
        if hasattr(self.args, 'sequencing_report'):
            self.report = self.check_sequencing_report()
        if hasattr(self.args, "region_plots"):
            self.region_of_interest = self.check_region_plots()
        
    def check_sequencing_report(self):
        assert self.args.sequencing_report.endswith(".csv"), "Your filename does not end with csv"
        assert os.path.isfile(self.args.sequencing_report), "Could not find the path to the sequencing report file."
        report = pd.read_csv(self.args.sequencing_report)
        assert isinstance(report, pd.DataFrame), "Could not read sequencing report. Please verify that your input file is a csv file."
        return report 
    
    def check_tsv_dir(self):
        assert os.path.isdir(self.args.tsv_dir) == True, "Path to directory with tsv files does not exist"

    def check_region(self):
        assert self.args.region in  ["targetSequences", "FR1", "CDR1", "FR2", "CDR2", "FR3", "CDR3", "FR4"]
            
            
    def check_length_threshold(self):
        assert self.args.length_threshold >=0, "Length threshold cannot be negative"
        
    def check_min_read_count(self):
        assert self.args.min_read_count >= 0, "Minimal read count cannot be negative"
        
    def check_remove_gaps(self):
        assert type(self.args.remove_gaps) == bool, "--remove_gaps must be boolean"
        
    def check_remove_errors(self):
        assert type(self.args.remove_errors) == bool, "--remove_errors must be boolean"
        
    def check_save_csv(self):
        assert self.args.save_csv.endswith(".csv") == True, "Your output filename must end with .csv"
        dir_csv = os.path.dirname(self.args.save_csv) 
        assert os.path.isdir(dir_csv), f"The directory for {os.path.basename(self.args.save_csv)} does not exist"
    
    def check_samples(self):
        assert len(self.args.samples) >= 1, "Please insert a sample name"
        for sample in self.args.samples:
            assert sample in self.report["Experiment"].unique(), f"Could not find {sample} in sequnecing_report"
        
    def check_pca_components(self):
        if self.args.batch_size:
            assert self.args.pca_components < self.args.batch_size, "Your number of pca components must be larger than your batch size"
        assert type(self.args.pca_components) == int
        assert self.args.pca_components > 2, "You must insert more than two pca components"
        
    def check_perplexity(self):
        assert self.args.perplexity >= 2, "Please enter more than 2 for perplexity"
        if self.args.batch_size:
            assert self.args.perplexity < self.args.batch_size, "Your batch size must be larger than the perplexity value"
            
    def check_iterations_tsne(self):
        assert self.args.iterations_tsne > 251, "Please enter more than 251 for the number of the iterations in tsne"
        
    def check_model_type(self):
        assert self.args.model_type in ["facebook/esm2_t6_8M_UR50D", "Rostlab/ProstT5_fp16", "Rostlab/prot_t5_xl_uniref50", "Rostlab/prot_t5_base_mt_uniref50", "Rostlab/prot_bert_bfd_membrane", "Rostlab/prot_t5_xxl_uniref50",
                       "Rostlab/ProstT5", "Rostlab/prot_t5_xl_half_uniref50-enc", "Rostlab/prot_bert_bfd_ss3", "Rostlab/prot_bert_bfd_localization",
                       "Rostlab/prot_t5_xl_bfd", "Rostlab/prot_bert", "Rostlab/prot_xlnet", "Rostlab/prot_bert_bfd", "Rostlab/prot_t5_xxl_bfd"], "please enter a valid model"
        
    def check_batch_size(self):
        assert type(self.args.batch_size) == int
                
    def check_binding_data(self):
        binding_data = self.args.binding_data
        if binding_data != None:
            assert binding_data.endswith("csv")
            df_bind = pd.read_csv(binding_data) 
            bind_cols = df_bind.columns.to_list()     
            assert self.region_of_interest in bind_cols,"Your region you want to analyse must exist in your binding data"
            assert len(bind_cols) >= 2,"Besides your sequences you must also have values."
            assert df_bind.shape[0] >= 1, "you must have at least one value"
                    
            
            
    def check_region_plots(self):
        region_of_interest = self.args.region_plots
        regions = ["aaSeqtargetSequences", "aaSeqFR1", "aaSeqCDR1", "aaSeqFR2", "aaSeqCDR2", "aaSeqFR3", "aaSeqCDR3", "aaSeqFR4"]
        assert region_of_interest in regions, f"Please enter a valid region of interest from {regions}"
        assert region_of_interest in self.report.columns.to_list(), "The region you would like to analyse does not exist in the current sequencing report. Please change that in the highest level."
        return region_of_interest
    
    def check_antigen_names(self):
        if self.args.antigen_names != None:
            assert self.args.binding_data != None, "You cannot parse a value here without having binding data."
            df_bind = pd.read_csv(self.args.binding_data)
            for antigen in self.args.antigen_names:
                assert antigen in df_bind.columns.to_list(), f"{antigen} is not in your binding data"
            
    def check_method_diversity(self):
        assert self.args.method in ["Shannon", "InverseSimpson"], "please enter a one of the following values: [Shanon, InverseSimpson] "
                
    
    
        
class ExpoSeqArgs:
    def __init__(self, **kwargs) -> None:
        self.parser = ArgumentParser(description = "Argument parser for plot function",
                                     allow_abbrev = True,
                                     **kwargs)
        self.chosen_tests = []
    def add_save_csv(self, flag = "--save_csv"):
        self.parser.add_argument(flag,
                                 help = "This must be parsed as path to tell where the output should be saved.",
                                 required = True)
        self.chosen_tests.append("check_save_csv")
        
    def add_sequencing_report(self, shortcut = "-r", flag = "--sequencing_report"):
        self.parser.add_argument(shortcut, 
                                 flag,
                                 help = "Insert the path to the sequencing report (merged csv file of all samples) here.",
                                 required = True)
        self.chosen_tests.append("check_sequencing_report")
    
    def add_samples(self, shortcut = "-s", flag = "--samples", required = True):
        self.parser.add_argument(shortcut,
                                 flag, 
                                 help = "Add a list of sample names here.",
                                 required = required,
                                 type = str,
                                 nargs = "+")
        if required:
            self.chosen_tests.append("check_samples")
        
    def add_region_of_interest(self, flag = "--region", default = "CDR3"):
        self.parser.add_argument(flag, 
                                 help = "This must be given to indicate which region should be analysed.",
                                 required=True,
                                 type = str)
        self.chosen_tests.append("check_region")
        
    def add_region_plots(self, flag = "--region_plots", default = "aaSeqCDR3"):
        self.parser.add_argument(flag,
                                 help = "This parses the region of interest starting with aaSeq to the plot.",
                                 type = str,
                                 required = True)
        self.chosen_tests.append("check_region_plots")
        
    def add_pca_components(self, flag = "--pca_components", default_value = 50):
        self.parser.add_argument(flag, 
                                 help = "Add an integer for the number of pca_components you want to reduce the embedding to",
                                 type = int,
                                 default = default_value)
        self.chosen_tests.append("check_pca_components")
        
    def add_perplexity(self, flag = "--perplexity", default_value = 20 ):
        self.parser.add_argument(flag,
                                 help ="Add an integer for the perplexity value the t-SNE algorithm needs to arrange your points in a certain cluster size.",
                                 type = int,
                                 default = default_value)
        self.chosen_tests.append("check_perplexity")
        
    def add_iterations_tsne(self, flag = "--iterations_tsne", default_value = 1500):
        self.parser.add_argument(flag, 
                                 help = "Add an integer for the number of iterations the t-SNE algorithm should take to arrange the coordiantes for the points.",
                                 default = default_value,
                                 type=int)
        self.chosen_tests.append("check_iterations_tsne")
        
    def add_batch_size(self, flag = "--batch_size", default_value = 600):
        self.parser.add_argument(flag, 
                                 help = "Add an integer for the batch size you want to take per sample which equals to the number of sequences with the top clone counts",
                                 type = int,
                                 default = default_value)
        self.chosen_tests.append("check_batch_size")
        
    def add_model_type(self, flag = "--model_type", default_value = 'Rostlab/prot_t5_xl_half_uniref50-enc'):
        self.parser.add_argument(flag,
                                 help = "Add a str as model type from the available models to embed your sequences in a high dimensional space.", 
                                 type = str,
                                 default= default_value)
        self.chosen_tests.append("check_model_type")
        
    def add_strands(self, flag = "--show_strands", default_value = True):
        self.parser.add_argument(flag, 
                                 help = "Enter a boolean to show a batch of sequences in your plots as identifier. False does not show any.",
                                 type = bool,
                                 default = default_value)
        
    def add_tsv_dir(self, flag = "--tsv_dir"):
        self.parser.add_argument(flag,
                                 help = "Add a path to a directory with tsv files here to generate the sequencing report.",
                                 type = str,
                                 required= True)
        self.chosen_tests.append("check_tsv_dir")
    
    def add_length_threshold(self, flag = "--length_threshold", default_value = 5):
        self.parser.add_argument(flag, 
                                 help = "Enter the threshold for the minimal sequence length in the sequencing report here.",
                                 default = default_value,
                                 type = int
                                 )
        self.chosen_tests.append("check_length_threshold")
        
    def add_min_read_count_threshold(self, flag = "--min_read_count", default_value = 3):
        self.parser.add_argument(flag,
                                 help = "Add a value for the minimum required read count certain sequences should have to remain in the sequencing report.",
                                 default = default_value,
                                 type = int)
        self.chosen_tests.append("check_min_read_count")
        
    def add_remove_gaps(self, flag = "--remove_gaps", default_value = False):
        self.parser.add_argument(flag,
                                 help = "Here you can enter whether you would like to remove gaps in your sequencing data. That means that your target sequence should be divisible by 3.",
                                 default=default_value,
                                 type = bool)
        self.chosen_tests.append("check_remove_gaps")
        
    def add_remove_errors(self, flag = "--remove_errors", default_value = True):
        self.parser.add_argument(flag,
                                 help = "If set to True this will remove all sequences (rows) in the sequencing report which have a *.",
                                 default=default_value,
                                 type = bool)
        self.chosen_tests.append("check_remove_errors")
        
    
    def add_binding_data(self, flag = "--binding_data", default_value = None):
        self.parser.add_argument(flag,
                                 help = "Add the path to the csv file with the binding data here, and it will be included in the plot.",
                                 type = str,
                                 default = default_value)
        self.chosen_tests.append('check_binding_data')
        
    def add_antigen_names(self, flag = "--antigen_names", default_value = None):
        self.parser.add_argument(flag, 
                                 help = "Your antigen names must be the column names of your binding data. You can insert multiple and they must be a string. You cannot parse names without having binding data.",
                                 nargs = "+",
                                 type = str,
                                 default = default_value)
        self.chosen_tests.append("check_antigen_names")
        
    def method_diversity(self, flag = "--method"):
        self.parser.add_argument(flag, 
                    help = "Choose the method to calculate the diversity in your samples. Default is Shannon for the ShannonIndex",
                    type = str,
                    default = "Shannon"
                    )
        self.chosen_tests.append("check_method_diversity")
        
def prep_args(args):
    parser = args.parser.parse_args()
    tests = args.chosen_tests
    TestInput = TestArgs(parser)
    for test in tests: # tests basically input values for chosen arguments. they should be all tested in the pipeline but this just prevents wrong inputs on a lower level.
        method = getattr(TestInput, test)
        method()
    return parser