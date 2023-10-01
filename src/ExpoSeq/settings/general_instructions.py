def print_instructions():
    instructions = '''

-----------------------
General Instructions
-----------------------
All plots can be called with 'plot' followed by the name of the plot, separated by a dot.
To analyze a specific sample or investigate possible options, use 'help()' and 'plot' followed by the name of the plot.

-----------------------
Plot Commands
-----------------------
- plot.alignment_quality(): 
    Creates a figure showing the Overall Reads per sample and the number of aligned reads.

- plot.aa_distribution(): 
    Returns a plot showing the amino acid distribution in the given sequence range.

- plot.rarefraction_curves(): 
    Shows the rarefraction curves for the given samples.

- plot.logoPlot_single(): 
    A logo Plot showing the composition of amino acids per position.

- plot.logoPlot_multiple(): 
    A logo Plot showing the composition of amino acids per position for multiple samples.

- plot.lengthDistribution_single(): 
    Shows the length distribution of the given sample.

- plot.relative_abundance_multi(): 
    Shows the relative abundance of the given samples.

- plot.lengthDistribution_multi(): 
    Shows the length distribution of the given samples.

- plot.rel_seq_abundance(): 
    Shows a bar plot of the frequencies of the most abundant sequences with optional Levenshtein distance.

- plot.basic_cluster(): 
    Clusters a batch of sequences based on a threshold for the Levenshtein distance.

- plot.cluster_one_AG(): 
    Clusters sequences based on Levenshtein distance and shows binding data against a specific antigen.

- plot.tsne_cluster_AG(): 
    Embeds sequences in a vector space, reduces dimensions, and clusters them. Also plots binding data.

- plot.embedding_tsne(): 
    Transforms sequences into a vector space and uses PCA and t-SNE for dimensionality reduction.

- plot.morosita_horn(): 
    Returns a matrix of the identity between your samples based on the Morosita Horn Index.

- plot.jaccard(): 
    Returns a matrix of the identity between your samples based on the Jaccard Index.

- plot.sorensen(): 
    Returns a matrix of the identity between your samples based on the Sorensen Index.

- plot.relative(): 
    Returns a matrix of the identity between your samples based on the Relative Index.

- plot.levenshtein_dendrogram(): 
    Returns a dendrogram of the sequences based on the Levenshtein distance.

- plot.dendro_bind(): 
    Shows a dendrogram based on Levenshtein distance and a bar plot with binding data.

- plot.MSA_design(): 
    Creates a multiple sequence alignment and offers an interface for sequence design.

-----------------------
Non-Plot Commands
-----------------------
- plot.add_binding_data(): 
    Allows you to add more binding data to the pipeline.

- plot.change_preferred_antigen(): 
    Sets the antigen you'd like to analyze in the plots.

- plot.change_preferred_sample(): 
    Sets the sample you'd like to analyze in the plots.

- plot.change_preferred_region(): 
    Sets the region you'd like to analyze in the plots.

- plot.change_filter(): 
    Adjusts the thresholds for filtering sequences in your analysis.

- plot.discard_samples(): 
    Allows you to remove samples from the analysis.

- plot.merge_bind_seq_report(): 
    Merges binding data with NGS data.

- plot.print_antigens(): 
    Shows the available antigen names in your data.

- plot.print_samples(): 
    Shows the available samples in your data.

- plot.save(): 
    Allows you to save a plot.

- plot.change_experiment_names(): 
    Allows you to change the names of the samples.
'''
    print(instructions)