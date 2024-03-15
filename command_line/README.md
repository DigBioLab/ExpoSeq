
# First step: Clone repo

```bash
git clone https://github.com/nilshof01/ExpoSeq
```

You can find all necessary scripts with:

```bash
cd src/ExpoSeq/command_line
```

However, the following guides will be all executed from the root dir of the repository.


# Highest layer: Sequencing report

Before you can start anything you need to prepare the sequencing report as input for all plots and functionalities in ExpoSeq. The output of the corresponding command line script is a csv file which contains all the data from your samples in only one table. If you clone the repository you can call the script with:

```bash
python command_line/cl_sequencing_report.py --tsv_dir "DIR_TO_TSV_FILES" --save_csv "C:\Users\nilsh\my_projects\ExpoSeq\tmp_test\test_report.csv" --region "CDR3"

```

- --tsv_dir: Directory to the tsv files with the tables from mixcr export tables
- -- save_csv: directory + filename ending with .csv for the output table
- --region: Name of the region you want to analyze. Available regions are: FR1, CDR1, FR2, CDR2, FR3, CDR3, FR4, targetSequences
- --length_threshold: Add a value for the minimum required read count certain sequences should have to remain in the sequencing report. Default is 3
- --remove_gaps: Here you can enter whether you would like to remove gaps in your sequencing data. That means that your target sequence should be divisible by 3. Default is True
- --remove_errors: If set to True this will remove all sequences (rows) in the sequencing report which have a *.

**NOTE:** If you decide to choose targetSequences as region, it is recommended to set remove_gaps to False, since all sequences with only one gap or more will be removed otherwise, which can be quiet substantial.

## Highest layer (optional): Binding Report

For testing purposes: You can generate binding data for your sequences with the following command line script:

```bash
python command_line/generate_binding_data.py --sequencing_report "C:\Users\nilsh\my_projects\ExpoSeq\tmp_test\test_report.csv" --save_csv "C:\Users\nilsh\my_projects\ExpoSeq\tmp_test\binding_data_test.csv" --samples "GeneMind_TRABkit_DNA77_300ng_repl1_L01_R1_001" "GeneMind_TRABkit_DNA80_300ng_repl1_L01_R1_001"
```


## Plots - Lower layers

### Embedding with t-SNE: Scatterplot

If you want to cluster your data without any binding data you can control the output as shown in the following:

**WITHOUT BINDING DATA**

```bash
python command_line/cl_protein_embedding.py -r "C:\Users\nilsh\my_projects\ExpoSeq\tmp_test\test_report.csv" --save_csv "C:\Users\nilsh\my_projects\ExpoSeq\tmp_test\embedding_umap_test.csv" --region_plots "aaSeqCDR3" --samples "GeneMind_TRABkit_DNA77_300ng_repl1_L01_R1_001" "GeneMind_TRABkit_DNA80_300ng_repl1_L01_R1_001" --batch_size 100 --pca_components 50 --perplexity 25 --iterations_tsne 252 --model_type "Rostlab/prot_t5_xl_half_uniref50-enc" --embedding_vector_path "temp/current_array.npz" --eps 0.1 --min_pts 2    
```

- --sequencing_report: Path to the report you generated in the highest layer
- --region string of the region you would like to analyse. It has to start with aaSeq. You cannot analyse other regions unless you change it in the highest layer under region
- --save_csv: directory + filename ending with .csv for the output table
- --batch_size: Number of sequences you would like to choose per sample for clustering.
- --pca_components: Default is 50. Has to be applied for better accuracy of t-SNE. You can indirectly change the described variance with this.
- --perplexity: Default is 20. It roughly determines the number of nearest neighbors that are considered in the embedding. A higher perplexity value results in a more global structure in the low-dimensional embedding, while a lower perplexity value emphasizes local structure. The optimal perplexity value for a given dataset depends on the dataset's intrinsic dimensionality, and it is usually determined by trial and error
- --iterations_tsne: Default is 1500. number of times that the algorithm will repeat the optimization process for reducing the cost function. The optimization process aims to minimize the difference between the high-dimensional and low-dimensional representations of the data. More iterations result in a more optimized low-dimensional representation, but also increases the computational cost.

### Embedding with t-SNE - with binding data: Scatterplot

**WITH BINDING DATA**

```bash
python command_line/cl_protein_embedding.py -r "C:\Users\nilsh\my_projects\ExpoSeq\tmp_test\test_report.csv" --save_csv "C:\Users\nilsh\my_projects\ExpoSeq\tmp_test\embedding_umap_test.csv" --region_plots "aaSeqCDR3" --samples "GeneMind_TRABkit_DNA77_300ng_repl1_L01_R1_001" "GeneMind_TRABkit_DNA80_300ng_repl1_L01_R1_001" --batch_size 100 --pca_components 50 --perplexity 25 --iterations_tsne 252 --model_type "Rostlab/prot_t5_xl_half_uniref50-enc" --embedding_vector_path "temp/current_array.npz" --eps 0.1 --min_pts 2 --binding_data "C:\Users\nilsh\my_projects\ExpoSeq\tmp_test\binding_data_test.csv" --antigen_names "Antigen 1"   
```

### Embedding with UMAP: Scatterplot

**WITHOUT BINDING DATA: Clustering multiple samples and find sequence clusters**


```bash
python command_line/cl_protein_embedding_umap.py -r "C:\Users\nilsh\my_projects\ExpoSeq\tmp_test\test_report.csv" --save_csv "C:\Users\nilsh\my_projects\ExpoSeq\tmp_test\embedding_umap_test.csv" --region_plots "aaSeqCDR3" --samples "GeneMind_TRABkit_DNA77_300ng_repl1_L01_R1_001" "GeneMind_TRABkit_DNA80_300ng_repl1_L01_R1_001" --batch_size 100 --pca_components 50 --n_neighbors 25 --min_dist 0.1 --metric "euclidean" --eps 0.5 --min_pts 2 --point_size 300 --model_type "Rostlab/prot_t5_xl_half_uniref50-enc" --embedding_vector_path "temp/current_array.npz" --n_jobs 1
```

- --save_csv: directory + filename ending with .csv for the output table
- --region_plots: string of the region you would like to analyse. It has to start with aaSeq. You cannot analyse other regions unless you change it in the highest layer under region
- --samples: Is required. You can enter as many samples as you want. The samples have to be in the report you are using.
- --batch_size: Number of sequences you would like to choose per sample for clustering.
- --pca_components: Default is 50. Has to be applied for better accuracy of UMAP. You can indirectly change the described variance with this.
- --n_neighbors: Default is 15. The size of local neighborhood (in terms of number of neighboring sample points) used for manifold approximation. Larger values result in more global views of the manifold, while smaller values result in more local data being preserved. In general, the value should be smaller than the number of samples.
- --min_dist: Default is 0.1. The effective minimum distance between embedded points. Smaller values will result in a more clustered/clumped embedding where nearby points on the manifold are drawn closer together, while larger values will result on a more even dispersal of points. The value should be set with respect to the scale of the data: in general smaller values will result in a more detailed embedding.
- --metric: Default is "euclidean". The metric to use to compute distances in high dimensional space. If X is a sparse matrix, it is recommended to set this to "cosine". If X is a dense matrix, it is recommended to set this to "euclidean".
- --eps and --min_pts: These are the parameters for the DBSCAN algorithm. The default values are 0.5 and 2. 
- --point_size: Default is 300. Adjusts the points in the plot relative to their clone fraction
- --model_type: Default is "Rostlab/prot_t5_xl_half_uniref50-enc". The model you would like to use for the embedding. 
- --embedding_vector_path: This is quiet useful since it can reduce the time for computing the plot drastically. The idea is that if you only change parameters for umap, pca or DBSCAN but leave the same parameters for the data you want to use the part with the sequence embedding is skipped and the report from the previous run with the same settings is used instead.
- ---n_jobs: Default is 1. The number of jobs to use for the computation. If -1, you will all of your available processes.

**WITHOUT BINDING DATA: Clustering multiple samples and supervise clustering based on certain sequence attribute**

```bash
python command_line/cl_protein_embedding_umap.py -r "C:\Users\nilsh\my_projects\ExpoSeq\tmp_test\test_report.csv" --save_csv "C:\Users\nilsh\my_projects\ExpoSeq\tmp_test\embedding_umap_test.csv" --region_plots "aaSeqCDR3" --samples "GeneMind_TRABkit_DNA77_300ng_repl1_L01_R1_001" "GeneMind_TRABkit_DNA80_300ng_repl1_L01_R1_001" --batch_size 100 --pca_components 50 --n_neighbors 25 --min_dist 0.1 --metric "euclidean" --eps 0.5 --min_pts 2 --point_size 300 --model_type "Rostlab/prot_t5_xl_half_uniref50-enc" --embedding_vector_path "temp/current_array.npz" --n_jobs 1 --characteristic "length"
```
Additional to the plot above you can add this parameter:

- --characteristic: Here you add the characteristic you want to use for the supervised clustering. You can choose one of the following:"isoelecrtric_point","aliphatic_index""hydrophobicity","weight","mass_charge_ratio","length","binding",

**WITH BINDING DATA: Clustering mutliple samples with functional assay data**

```bash
python -m command_line/cl_protein_embedding_umap.py -r "C:\Users\nilsh\my_projects\ExpoSeq\tmp_test\test_report.csv" --save_csv "C:\Users\nilsh\my_projects\ExpoSeq\tmp_test\embedding_umap_test.csv" --region_plots "aaSeqCDR3" --samples "GeneMind_TRABkit_DNA77_300ng_repl1_L01_R1_001" "GeneMind_TRABkit_DNA80_300ng_repl1_L01_R1_001" --batch_size 100 --pca_components 50 --n_neighbors 25 --min_dist 0.1 --metric "euclidean" --eps 0.5 --min_pts 2 --point_size 300 --model_type
 "Rostlab/prot_t5_xl_half_uniref50-enc" --embedding_vector_path "temp/current_array.npz" --n_jobs 1 --binding_data "C:\Users\nilsh\my_projects\ExpoSeq\tmp_test\binding_data_test.csv" --antigen_names "Antigen 1"
```

- --binding_data: Path to the binding data you want to use for the embedding. The binding data has to be in a csv file with the following columns: "aaSeqCDR3", "ANTIGEN1", "ANTIGEN2". The binding value has to be a float value.
- --antigen_names: Here you can add the names of the antigens you want to use for the embedding. You can enter as many antigens as you want. The names have to be in the binding data you are using.


### Diversity Plot: Barplot

```bash
python command_line/cl_diversity_plot.py --sequencing_report "PATH_TO_CSV_FILE" --region "aaSeqCDR3" --save_csv "OUTPUT_FILENAME"
```


### Rarefraction curves: line Plot

```bash
python command_line/cl_rarefraction_curves.py -r "PATH_TO_CSV_FILE" --region_plots "aaSeqCDR3" --save_csv "OUTPUT_FILENAME" --samples "SAMPLE_NAME1" "SAMPLE_NAME2"
```

- --samples: Is not reuqired, so if you do not enter the flag it will take all available samples.

### Data for logo plots: LOGOPLOT or stacked bar plot

```bash
python command_line/cd_logoplot.py -r "C:\Users\nilsh\my_projects\ExpoSeq\tmp_test\test_report.csv" --save_csv "C:\Users\nilsh\my_projects\ExpoSeq\tmp_test\logoplot_test.csv" --region_plots "aaSeqCDR3" --samples "GeneMind_TRABkit_DNA77_300ng_repl1_L01_R1_001" --method "proportional"
```

- --save_csv: directory + filename ending with .csv for the output table
- --region_plots: string of the region you would like to analyse. It has to start with aaSeq. You cannot analyse other regions unless you change it in the highest layer under region
- --samples: Is required. You can enter as many samples as you want. The samples have to be in the report you are using.
- --method: Default is "proportional". You can also choose "bits" which will calculate the output based on bits.

--> YOU can plot this also as stacked barplot

### Length Distribution data : Barplot

```bash
python command_line/cl_length_distribution.py -r "C:\Users\nilsh\my_projects\ExpoSeq\tmp_test\test_report.csv" --save_csv "C:\Users\nilsh\my_projects\ExpoSeq\tmp_test\length_distribution_test.csv" --region_plots "aaSeqCDR3" --single_sample "GeneMind_TRABkit_DNA77_300ng_repl1_L01_R1_001"
```

- --save_csv: directory + filename ending with .csv for the output table
- --region_plots: string of the region you would like to analyse. It has to start with aaSeq. You cannot analyse other regions unless you change it in the highest layer under region
- --single_sample: Is required. You can enter only one sample. The sample has to be in the report you are using.


## Levenshtein Clustering 2D

### Levenshtein Network

```bash
python -m command_line.cl_levenshtein_clustering --save_csv "C:\Users\nilsh\my_projects\ExpoSeq\tmp_test\levenshtein_2d.csv" --region_plots "aaSeqCDR3" --samples "GeneMind_TRABkit_DNA77_300ng_repl1_L01_R1_001" -r "C:\Users\nilsh\my_projects\ExpoSeq\tmp_test\test_report.csv" --batch_size 100 --binding_data "C:\Users\nilsh\my_projects\ExpoSeq\tmp_test\binding_data_test.csv" --antigen_names "Antigen 1
```
Note: You can also run that without binding data by leaving out the corresponding flags (antigen_names, binding_data)

### Levenshtein Histogram

```bash
python -m command_line.cl_lvst_histogram --save_csv "C:\Users\nilsh\my_projects\ExpoSeq\tmp_test\levenshtein_dendro.csv" --region_plots "aaSeqCDR3" --samples "GeneMind_TRABkit_DNA77_300ng_repl1_L01_R1_001" -r "C:\Users\nilsh\my_projects\ExpoSeq\tmp_test\test_report.csv" --batch_size 200
```

