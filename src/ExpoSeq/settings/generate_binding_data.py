import pandas as pd
import numpy as np 
from argparse import ArgumentParser

parser = ArgumentParser(
    description="Argument parser generating binding data", allow_abbrev=True 
)
parser.add_argument(
                    "--sequencing_report",
                    help = "Path to the sequencing report",
                    type = str,
                    required = True)
parser.add_argument("--save_csv",
                    help = "Path to save the csv",
                    type = str,
                    required = True)
parser.add_argument("--samples",
                    help = "Samples to include in the plot",
                    nargs = "+",
                    required = True)
parser.add_argument("--no_sequences",
                    help = "Number of sequences to generate",
                    type = int,
                    default = 10)

args = parser.parse_args()
sequencing_report = args.sequencing_report
sequencing_report = pd.read_csv(sequencing_report)
report = sequencing_report[sequencing_report['Experiment'].isin(args.samples)]
grouped = report.groupby('Experiment')
# Initialize an empty list to store the top sequences for each sample
top_sequences = []
top_n = args.no_sequences
binding_data = pd.DataFrame(columns=[ 'aaSeqCDR3', 'Antigen 1'])
# Iterate over each group (sample) and select the top 10 sequences
for sample, group_data in grouped:
    top_sequences_per_sample = group_data['aaSeqCDR3'].value_counts().nlargest(top_n).index.tolist()
    values = np.random.randint(low=1000, high=1000000, size=len(top_sequences_per_sample))  
    sample_values = pd.DataFrame({'aaSeqCDR3': top_sequences_per_sample, 'Antigen 1': values})
    top_sequences.extend(top_sequences_per_sample)
    binding_data = pd.concat([binding_data, sample_values])

binding_data.to_csv(args.save_csv, index=False)