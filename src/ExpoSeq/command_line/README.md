# Highest layer: Sequencing report

Before you can start anything you need to prepare the sequencing report as input for all plots and functionalities in ExpoSeq. The output of the corresponding command line script is a csv file which contains all the data from your samples in only one table. If you clone the repository you can call the script with:

```bash
python -m src.ExpoSeq.command_line.cl_sequencing_report --tsv_dir "DIR_TO_TSV_FILES" --save_csv "OUTPUT_FILENAME" --region "CDR3"

```

- --tsv_dir: Directory to the tsv files with the tables from mixcr export tables
- -- save_csv: directory + filename ending with .csv for the output table
- --region: Name of the region you want to analyze. Available regions are: FR1, CDR1, FR2, CDR2, FR3, CDR3, FR4, targetSequences
- --length_threshold: Add a value for the minimum required read count certain sequences should have to remain in the sequencing report. Default is 3
- --remove_gaps: Here you can enter whether you would like to remove gaps in your sequencing data. That means that your target sequence should be divisible by 3. Default is True
- --remove_errors: If set to True this will remove all sequences (rows) in the sequencing report which have a *.

**NOTE:** If you decide to choose targetSequences as region, it is recommended to set remove_gaps to False, since all sequences with only one gap or more will be removed otherwise, which can be quiet substantial.

## Highest layer (optional): Binding Report



## Plots - Lower layers

### Embedding with t-SNE

If you want to cluster your data without any binding data you can control the output as shown in the following:

**WITHOUT BINDING DATA**

```bash
python -m src.ExpoSeq.command_line.cl_protein_embedding --sequencing_report "PATH_TO_CSV_FILE" --region "aaSeqCDR3" --save_csv "OUTPUT_FILENAME" --samples "SAMPLE_NAME1" "SAMPLE_NAME2" --batch_size 100 --pca_components 50 --perplexity 25 --iterations_tsne 252 
```

- --sequencing_report: Path to the report you generated in the highest layer
- --region string of the region you would like to analyse. It has to start with aaSeq. You cannot analyse other regions unless you change it in the highest layer under region
- --save_csv: directory + filename ending with .csv for the output table
- --batch_size: Number of sequences you would like to choose per sample for clustering.
- --pca_components: Default is 50. Has to be applied for better accuracy of t-SNE. You can indirectly change the described variance with this.
- --perplexity: Default is 20. It roughly determines the number of nearest neighbors that are considered in the embedding. A higher perplexity value results in a more global structure in the low-dimensional embedding, while a lower perplexity value emphasizes local structure. The optimal perplexity value for a given dataset depends on the dataset's intrinsic dimensionality, and it is usually determined by trial and error
- --iterations_tsne: Default is 1500. number of times that the algorithm will repeat the optimization process for reducing the cost function. The optimization process aims to minimize the difference between the high-dimensional and low-dimensional representations of the data. More iterations result in a more optimized low-dimensional representation, but also increases the computational cost.

### Embedding with t-SNE - with binding data

**WITH BINDING DATA**

```bash
python -m src.ExpoSeq.command_line.cl_protein_embedding --sequencing_report "PATH_TO_CSV_FILE" --region "aaSeqCDR3" --save_csv "OUTPUT_FILENAME" --samples "SAMPLE_NAME1" "SAMPLE_NAME2" --batch_size 100 --pca_components 50 --perplexity 25 --iterations_tsne 252 --binding_data "PATH_TO_CSV_FILE" --antigen_names "ANTIGEN1" "ANTIGEN2"
```


### Diversity Plot

```bash
python -m src.ExpoSeq.plots.command_line.cl_diversity_plot --sequencing_report "PATH_TO_CSV_FILE" --region "aaSeqCDR3" --save_csv "OUTPUT_FILENAME"
```


### Rarefraction curves

```bash
python -m src.ExpoSeq.command_line.cl_rarefraction_curves -r "PATH_TO_CSV_FILE" --region_plots "aaSeqCDR3" --save_csv "OUTPUT_FILENAME"
```

