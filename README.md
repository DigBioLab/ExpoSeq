# Welcome to ExpoSeq

ExpoSeq is a powerful pipeline for processing and analyzing FASTQ files from sequencing phage Display panning samples. It utilizes [MiXCR](https://docs.milaboratories.com/mixcr/getting-started/installation/) to align and assemble the data which you can subsequently analyze in multiple plots. The pipeline focuses on analysing the identity between samples but also applies various clustering techniques to analyse the relation between the sequences. Besides, you can add binding data to relate the clusters to affinity.  ![overview](pictures_gen/expoSeq_overview.png)

## Installation

Open a virtual environment and type ```pip install ExpoSeq```. Ensure that you have python > 3.11 installed.

To get started, please download and follow the instructions for MiXCR under the following link: https://docs.milaboratories.com/mixcr/getting-started/installation/ 
You can also only use the test version of ExpoSeq without installing it.

## Importing the Plotting Tool

To access the plotting tool, you will need to import it into your console by running the following command:
<br>
```python
from ExpoSeq.pipeline import PlotManager
```

## Using the PlotManager

The PlotManager is the main interface for creating various plots using your FASTQ data. You can create an instance of the PlotManager by running the following command:
<br>

```python
plot = PlotManager()
```
After that you will be automatically guided through the data processing and preparation. As soon as this has been finished you the pipeline will automatically create an analysis of your data and will store the plots in ~/my_experiments/YOUR_EXPERIMENT_NAME/plots .
<br>
If you want to create some plots by yourself, please take a look at the [Jupyter script](ExpoSeq_handsOn.ipynb).
<br>

In the following you can obtain an insight in the worklow of the pipeline after the initial call. There, the blue boxes indicate your input, gray are optional inputs while black and red are processing steps and output, respectively.
<br>
![relative_path_to_image](pictures_gen/workflow_ExpoSeq.png)
<br>
If you just want to test the pipeline and see its functions you can call:
<br>

```python
plot = PlotManager(test_version = True)
```
<br>
If you would like to have details about the inputs and functions of the PlotManager call:
<br>

```python
help(plot)
``` 

You can also call for specific plots, for instance:
<br>

```python
help(plot.jaccard)
```

## Upload binding data 

If you have conducted DELFIA or other techniques to receive binding data for certain sequences (usually sanger sequenced), you can upload these in a certain format and use these for clustering to potentially find other suitable sequences with high binding.  You need to import the data as csv file where the first column starts in the first row with the header: aaSeqCDR3 which are the sequences. It is very important to keep the header at this position. In the second column you can put the binding data for your epitope which you can name in the first row however you prefer. You can have a look in [this csv file](src/ExpoSeq/test_data/test_files/binding_data.csv) to see the general structure of the file. Moreover, you can download it and import it in Excel. Therefore, open Excel and choose under "Data" in the Excel header "From Text/CSV". Then make sure to delete the first column which contains the row number. After that you can delete the random data in that excel sheet and add your own. Finally you can export the data as a csv and import it with the pipeline either in the initial uploading process which will be prompted or with the command 
<br>

```python
plot.add_binding_data()
```
<br>

**Note**: If you decide to add more binding data to your analysis you can just use the same command and choose the new file with the filechooser and it will be added to the existing data.

## Data processing on a Cluster

First pull the folder with the scripts for the processing to your working directory

```bash
wget https://github.com/nilshof01/ExpoSeq/tree/final_master/bash_processing
```

I have prepared and example jobscript for working on an LSF cluster. You can have a look under ~/bash_processing/example_LSF_cluster.sh . 
To run your script interactively you can call:

```bash
python ~/bash_processing/mixcr_cl.py $PATH_TO_MIXCR $YOUR_EXPERIMENT_NAME $PATH_TO_FORWARD_FILES 
```

**NOTE**: You need to have mixcr installed in your working directory to be able to start the processing.
To use multithreading and increase the RAM allocation have a look at the following arguments you can input:

- `--path_to_mixcr`: Is the filepath to the mixcr.jar file.
- `--experiment_name`: A string which is the name of your experiment.
- `--path_to_forward`: The directory of the fastq files with the forward reads. 
- `--path_to_backward`: The directory of the fastq files with the backward reads.
- `--threads`: The number of threads you would like to use for the processing.
- `--method`: The mixcr method to align and assemble the reads you would like to use. Default is milab-human-tcr-dna-multiplex-cdr3
- `--java_heap_size`: Memory for proecssing in MB.

**NOTE**: If you only want to process forward reads, then you do not need to add the path to the directory with the backward reads. Further, if you would like to analyze paired end sequencing data, please make sure that forward and backward fastq files are in separate folders.

After the processing has been finished you can use the directory ~/my_experiments/YOUR_EXPERIMENT_NAME for the plotmanager. As soon as you call it you can choose to upload the directory.

## References
[1] Dmitriy A. Bolotin, Stanislav Poslavsky, Igor Mitrophanov, Mikhail Shugay, Ilgar Z. Mamedov, Ekaterina V. Putintseva, and Dmitriy M. Chudakov. "MiXCR: software for comprehensive adaptive immunity profiling." Nature methods 12, no. 5 (2015): 380-381.


[2] Dmitriy A. Bolotin, Stanislav Poslavsky, Alexey N. Davydov, Felix E. Frenkel, Lorenzo Fanchi, Olga I. Zolotareva, Saskia Hemmers, Ekaterina V. Putintseva, Anna S. Obraztsova, Mikhail Shugay, Ravshan I. Ataullakhanov, Alexander Y. Rudensky, Ton N. Schumacher & Dmitriy M. Chudakov. "Antigen receptor repertoire profiling from RNA-seq data." Nature Biotechnology 35, 908–911 (2017)

[3] (1, 2) Tareen A, Kinney JB (2019) Logomaker: beautiful sequence logos in Python. Bioinformatics btz921. bioRxiv doi:10.1101/635029.

[4] M.A. Larkin and others, Clustal W and Clustal X version 2.0, Bioinformatics, Volume 23, Issue 21, November 2007, Pages 2947–2948, https://doi.org/10.1093/bioinformatics/btm404





