from gensim.models import Word2Vec
from gensim.models.keyedvectors import KeyedVectors
from Bio import SeqIO
import numpy as np

from Bio.Align import substitution_matrices
# Load the BLOSUM62 matrix from Biopython
blosum = np.array(substitution_matrices.load("BLOSUM62"), dtype=int)
blosum = blosum[:20, :20]
# Define a function to generate training examples

fraction = 0.01

# group your data by sample
grouped = sequencing_report.groupby(['Experiment', 'aaSeqCDR3'])

# loop over each sample and sample the desired fraction of sequences
samples = []
for name, group in grouped:

    sampled_group = sequencing_report.groupby("Experiment")["aaSeqCDR3"].head(100).sample(frac=1)
    # append the sampled sequences to the list of samples
    samples.append(sampled_group)

# concatenate the samples into a new DataFrame
sampled_data = pd.concat(samples)



def generate_training_examples(seq, window=5):
    examples = []
    collect_all = []
    for i, target in enumerate(seq):
        for j in range(max(0, i- window), min(len(seq), i + window + 1)):
            if i != j:
                context = seq[j]
                collect_all.append((target, context))
    sampling = random.sample(collect_all, 10)
    examples.extend(sampling)


    return examples


# Define a function to generate embeddings
def generate_embeddings(sequences_full, embedding_size=100, window=6, iter=5):
    # Read the protein sequences from a FASTA file

    aminoacids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    # Generate the training examples
    training_examples = []
    for seq in sequences_full:
        seq_vec = []
        for aa in seq:
            if aa in aminoacids:
                seq_vec.append(list(blosum[aminoacids.index(aa)]))
        training_examples.extend(generate_training_examples(seq_vec, window=window))

    batch = []
    for tup in training_examples:
        for word in tup:
            batch.append(word)
    # Train the Word2Vec model
    model = Word2Vec(vector_size=embedding_size, window=window, min_count=1, sg=1)
    model.build_vocab(batch)
    model.train(batch, total_examples=model.corpus_count, epochs=model.epochs)

    # Save the embeddings to a file
    model.wv.save_word2vec_format('protvec.txt')

    # Load the embeddings into a KeyedVectors object
    embeddings = KeyedVectors.load_word2vec_format('protvec.txt', binary=False)

    return embeddings