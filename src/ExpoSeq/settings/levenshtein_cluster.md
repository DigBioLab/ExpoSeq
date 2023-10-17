### Clustering based on levenshtein distance 


The Levenshtein distance is a metric that measures the similarity between two strings. It is calculated by counting the minimum number of edits (insertions, deletions, and substitutions) required to transform one string into the other.

For example, the Levenshtein distance between the strings "CAT" and "DOG" is 2, because we need to insert the letter "D" and substitute the letter "C" with the letter "G" to transform "CAT" into "DOG".

The Levenshtein distance can be used to analyze peptide sequences in a variety of ways. For example, it can be used to:

1. Identify similar CDR3 regions. This can be useful for identifying antibodies with similar antigen-binding specificities.
2. Detect mutations in CDR3 regions. This can be useful for identifying mutations that may affect the affinity or specificity of an antibody.
3. Cluster CDR3 regions into groups. This can be useful for identifying groups of CDR3 regions with similar antigen-binding specificities.


The connected components in the plot (shown by a a network of lines) are groups of CDR3 sequences that are similar to each other, based on a threshold for a levenshtein distance (default = 2). The size of each circle in the plot represents the number of clones for the corresponding CDR3 sequence relative to the clone counts of other sequences.

Overall, the plot can be used to get a general overview of the diversity and similarity of the CDR3 sequences in the library.

Furthermore, a table is given which has the sorted sequences based on the clone count. That means the top 10 sequences represent the sequences with the highest clone fraction in that sample. This table, coming from a report can be used to identify specificity for a certain antigen based on a density of these sequences in certain clusters.

Disadvantages of the levenshtein distance

- Only captures global similarity of sequences and does not focues on specific regions
- Does not take any chemical properties or structural properties into account