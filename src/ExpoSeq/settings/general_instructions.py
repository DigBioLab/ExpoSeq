def print_instructions():
    instructions = '''

-----------------------
General Instructions
-----------------------
All plots can be called with 'plot' followed by the name of the plot, separated by a dot.
To analyze a specific sample or investigate possible options, use 'help()' and 'plot' followed by the name of the plot.

-----------------------
Important Variables
-----------------------

- plot.sequencing_report 
    Contains the processed and trimmed NGS data.
    
- plot.binding_data
    Contains the uploaded binding data.
    
- plot.Report.origin_seq_report
    Contains the raw NGS data from mixcr which is not trimmed.
    
- plot.ControlFigure.fig
    Is the matplotlib figure object
    
- plot.ControlFigure.ax
    Is the matplotlib axis object
    
- plot.alignment_report
    A dataframe which contains the data from all alignment reports.
    
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

- plot.lengthDistribution_multi(): 
    Shows the length distribution of the given samples.

- plot.rel_seq_abundance(): 
    Shows the clone fractions of a given sample in a tree map plot

- plot.ls_distance_binding():
    Advanced and new version of clustser_one_AG. You can include multiple samples and antigens here. 
    You will cluster sequences based on Levenshtein distance and combine them with your binding data.

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
    
- plot.cluster_binding_data():
    Can cluster multiple samples and sequences with binding data from multiple antigens. Integrates sequence embedders from huggingface.

- plot.connect_samples():
    Connects all samples in a network plot based on the Levenshtein distance.
    
- plot.sample_diversity()
    Shows the diversity of your samples based on the Inverse Simpson Index or the Shannon Index.
    
- plot.length_distribution_all()
    Shows you the length distribution for all of your samples separately in one plot as Violin or Boxplot.
-----------------------
Non-Plot Commands
-----------------------
- plot.add_binding_data(): 
    Allows you to add more binding data to the pipeline.

- plot.change_preferred_antigen(): 
    Sets the antigen you'd like to analyze in the default plots.

- plot.change_preferred_sample(): 
    Sets the sample you'd like to analyze in the default plots.

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
    
- plot.chat()
    Allows you to use artifical intelligence to investigate your data without coding. 
   
- plot.show()
    Sometimes the plot does not show up. You can force the pipeline to do so by typing this command. 

'''
    print(instructions)