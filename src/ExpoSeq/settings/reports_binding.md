### Sequence embedding and T-SNE - REPORT

Additionally to the plot, a report was created which contains the results of the t-SNE embedding. This report contains the following columns:
- tsne1: Contains the coordinates for the first dimension (x-axis) of the plot
- tsne2: Contains the coordinates for the second dimension (y-axis) of the plot
- experiments_factorized: Contains a factorized form for the experiment column. You can ignore this column
- experiments_string: Contains the sample name from the NGS data, used for this plot.
- binding: Contains the binding value which was taken from the binding data with the corresponding antigen. NGS sequences have the value 0
- sequences: Amino Acid sequence for the given coordinates and binding value
- sequence_id: Can be used to identify the position of the sequence in the right plot with the numbers

The plot can be used to identify sequences which could have a potential high binding based on the values of neighboring sequences. Further, cluster may be identified which have a high local binding and sequence attributes may be derived based on this.

**Levenshtein Distance based Clustering - REPORT**

The Report for this plots can be used to localize sequence in the network plots. The report contains the following columns:
- Clusters_[YOUR_SAMPLENAME]_No: Contains the cluster number for the given sequence. This is useful if you want to see the sequences which are localized in the same cluster and thus have a high sequence similarity.
- Sequences_[YOUR_SAMPLE_NAME]: The sequences for the corresponding cluster number

This plot can be particular useful if you look for potential high binders only based on point mutations.