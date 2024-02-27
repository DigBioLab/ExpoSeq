from argparse import ArgumentParser

class ExpoSeqArgs:
    def __init__(self, kwargs) -> None:
        self.parser = ArgumentParser(description = "Argument parser for plot function",
                                     allow_abbrev = True,
                                     **kwargs)
    def save_csv(self, flag = "--save_csv"):
        self.parser.add_argument(flag,
                                 help = "This must be parsed as path to tell where the output should be saved.",
                                 required = True)
        
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
    
        
        
    
        