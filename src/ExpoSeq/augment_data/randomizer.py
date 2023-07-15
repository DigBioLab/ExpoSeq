import numpy as np
import random
import string
import os
import pandas as pd
from Bio.Seq import Seq
import pickle

import math
def create_fractions(num_fractions, top_fraction):
    total_sum = 1
    decrement_factor = math.pow(top_fraction, 1/(num_fractions-1))

    # Create a list of fractions that decrease geometrically
    fractions = [top_fraction * math.pow(decrement_factor, i) for i in range(num_fractions)]
    fractions_sum = sum(fractions)

    # Normalize the fractions so they sum to one
    fractions = [fraction/fractions_sum for fraction in fractions]
    return fractions


def create_sequencing_report(num_experiments = 30, mean_length = 48, stddev_length = 36,panning_rounds = 3, threshold_seq_length = 5 ):
    # Define the stop codons
    stop_codons = ["TAA", "TAG", "TGA"]
    num_experiments = round(num_experiments/panning_rounds)
    sequencing_report = pd.DataFrame([])
    for experiment in range(num_experiments):
        # Define the set of possible nucleotides
        nucleotides = "ATCG"
        sequences = []
        num_sequences = random.randint(1000, 20000)
        # Generate multiple random DNA sequences
        for i in range(num_sequences):
            length = int(np.random.normal(mean_length, stddev_length))
            length = length - (length % 3)  # Round down to the nearest multiple of 3
            while True:
                sequence = ''.join(random.choice(nucleotides) for _ in range(length))
                if not any(sequence[i:i + 3] in stop_codons for i in range(0, len(sequence), 3)):
                    break


            dna_seq = sequence
            sequences.append(dna_seq)
            dna_seq = Seq(dna_seq)
            aa_seq = dna_seq.translate()
        sequences = [string for string in sequences if len(string) >= threshold_seq_length]
        aminoacid_seqs = [str(Seq(string).translate()) for string in sequences]
        num_sequences = len(sequences)
        sample_name = ["sample " + str(experiment)] * num_sequences
        random_float = random.uniform(0, 1) * 0.1
        fractions = create_fractions(num_sequences, random_float)
        random_read = random.randint(1000000, 10000000)
        readsCount = np.array(fractions) * random_read
        readsCount = list(np.round(readsCount).astype(np.uint32))
        length_cdr3 = [len(string) for string in sequences]
        data = {"Experiment": sample_name,
                    "aaSeqCDR3": aminoacid_seqs,
                    "nSeqCDR3": sequences,
                    "clonesFraction": fractions,
                    "readCount": readsCount,
                    "lengthOfCDR3": length_cdr3}
        intermediate_frame = pd.DataFrame(data)
        sequencing_report = pd.concat([sequencing_report, intermediate_frame])
        for panning_round in range(panning_rounds-1):
            random_float = random.uniform(0.6, 1)
            fraction_to_sample = 0.5
            num_values_to_sample = int(fraction_to_sample * len(sequences))
            sequences_new = random.sample(sequences, num_values_to_sample)
            aminoacid_seqs_new = [str(Seq(string).translate()) for string in sequences_new]
            num_sequences = len(sequences_new)
            sample_name = ["sample " + str(experiment) + " round " + str(panning_round + 2)] * num_sequences
            random_float = random.uniform(0, 1) * 0.1
            fractions_new = create_fractions(num_sequences, random_float)
            random_read = random.randint(1000000, 10000000)
            readsCount_new = np.array(fractions_new) * random_read
            readsCount_new = list(np.round(readsCount_new).astype(np.uint32))
            length_cdr3_new = [len(string) for string in sequences_new]
            data = {"Experiment": sample_name,
                    "aaSeqCDR3": aminoacid_seqs_new,
                    "nSeqCDR3": sequences_new,
                    "clonesFraction": fractions_new,
                    "readCount": readsCount_new,
                    "lengthOfCDR3": length_cdr3_new}
            intermediate_frame = pd.DataFrame(data)
            sequencing_report = pd.concat([sequencing_report, intermediate_frame])
    return sequencing_report



def create_binding_report(sequencing_report, num_antigen):
    threshold_aa_seq = 5
    sampled_df = sequencing_report[~sequencing_report.duplicated(subset='aaSeqCDR3', keep=False)]
    sampled_df = sampled_df.groupby("Experiment").sample(frac = 0.01)
    sampled_df = sampled_df[sampled_df['aaSeqCDR3'].apply(len) >= threshold_aa_seq]
    aaSeq = sampled_df["aaSeqCDR3"].to_list()
    fractions = sampled_df.clonesFraction
    binding_values = []
    for i in fractions:
        random_binding = random.randint(1000000, 10000000000)
        binding_value = i * random_binding
        binding_values.append(binding_value)
    num_splits = num_antigen
    # generate random indices where the splits should occur
    split_indices = np.sort(np.random.choice(range(1, len(binding_values)), num_splits, replace=False))
    # split the list into random parts
    split_lists = [(binding_values[i:j], aaSeq[i:j]) for i, j in zip([0]+split_indices, split_indices+[len(binding_values)])]
    for i in range(len(split_lists)):
        part_binding = split_lists[i][0]
        part_aaSeq = split_lists[i][1]
        data = {"aaSeqCDR3": part_aaSeq, "Antibody " + str(i): part_binding}
        intermediate_binding = pd.DataFrame(data)
        if i == 0:
            binding_data = intermediate_binding
        else:
            binding_data = pd.merge(binding_data, intermediate_binding,how = "left",on = "aaSeqCDR3")
    
    return binding_data


def create_test_dataset(module_dir):
    sequencing_report = create_sequencing_report()
    binding_data = create_binding_report(sequencing_report, 5)
    binding_data_filepath = os.path.join(module_dir, "test_data", "test_files", "binding_data.csv")
    assert os.path.isfile(binding_data_filepath), "binding data filepath does not exist"
    binding_data.to_csv(binding_data_filepath)
    sequencing_report_filepath = os.path.join(module_dir, "test_data", "test_files", "sequencing_report.csv")
    assert os.path.isfile(sequencing_report_filepath), "sequencing report filepath does not exist"
    sequencing_report.to_csv(sequencing_report_filepath)
    unique_experiments = sequencing_report["Experiment"].unique()
    experiment_dic = {item: item for item in list(unique_experiments)}
    exp_names_dir = os.path.join(module_dir, "test_data", "test_files", "experiment_names.pickle")
    assert os.path.isfile(exp_names_dir), "experiment names filepath does not exist"
    with open(exp_names_dir, "wb") as f:
        pickle.dump(experiment_dic, f)
        
