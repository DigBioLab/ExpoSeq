# Welcome to ExpoSeq

ExpoSeq is a powerful pipeline for processing and analyzing FASTQ files from high-throughput sequencing of samples from display techniques such as phage display. The pipeline was mainly designed for analyzing phage display campaigns, but RNA-Seq data for analysing T-cell receptor data can also be given as input. It utilizes [MiXCR](https://docs.milaboratories.com/mixcr/getting-started/installation/) to align and assemble the data which you can subsequently analyze in multiple plots. The pipeline offers a huge flexibility for various sequencing strategies and targets which are defined by the presets which are offered by MiXCR and are implemented in the pipeline. You can take a look at the available presets [here](src/ExpoSeq/settings/available_presets.py). The pipeline focuses on analysing the quality of your sequencing and the sequence relations between your samples. Further, you can upload binding data to the pipeline to analyse these with various clustering techniques.

## Installation

Make sure you have installed [Python](https://www.python.org/downloads/) on your system. After that you can install ExpoSeq in the terminal with
```
pip install ExpoSeq
```
 Ensure that you have python > 3.9 installed.

To get started, please download and follow the instructions for MiXCR at their [official documentation](https://docs.milaboratories.com/mixcr/getting-started/installation/ )
You can also only use the test version of ExpoSeq without installing it.

## Importing the Plotting Tool

To access the plotting tool, you will need to import it into your console by running the following command:
<br>
```python
from ExpoSeq.pipeline import PlotManager
```
The PlotManager is the main interface for creating various plots using your FASTQ data. You can create an instance of the PlotManager by running the following command:
<br>

```python
plot = PlotManager()
```

After that you will be automatically guided through the data processing and preparation. As soon as this has been finished the pipeline will automatically create an analysis of your data and will store the plots in
<br>
```
📦my_experiments
 ┣ 📂YOUR_EXPERIMENT_NAME
 ┃ ┣ 📂plots
 ┃ ┃ ┗ 📂aaSeqCDR3
```
<br>

If you would like to know more about ExpoSeq's generated data structure, please have a look [here](bash_processing/utilize_all_data.md).


### Generate automatically a report of your whole analysis (optional)

The pipeline is able to generate a whole report automatically based on the results generated by the pipeline. You can have a look at a generated example [here](additional_material/test_show.html) by downloading and opening the html file. Before you can use this funcitonality you need to install [Quarto](https://quarto.org/docs/get-started/) on your system. It was used with [VS Code](https://code.visualstudio.com/). After you have installed quarto you can run: 
<br>
```python
plot.create_report()
```
The pipeline will automatically check if Quarto exists on your device. You can use the [quarto file](additional_material/test_show.qmd) generated by the pipeline to add your findings and individualize the report by utilizing the strength of markdown. If you have made changes to the file you can render it with quarto by running:
<br>
```
quarto render YOUR/PATH/TO/FILE.qmd
```

**NOTE**: Do not change the default structure of the files in the folder with the plots which is generated by the pipeline. Otherwise, this functionality might not work. Further the pipeline should contain the binding data if you have uploaded it and created the default plots with it.


## Individualize Plots by using the PlotManager (optional)

<br>

If you want to create some plots by yourself, please take a look at the [jupyter script](additional_material/ExpoSeq_handsOn.ipynb) .
<br>

In the following you can obtain an insight in the worklow of the pipeline after the initial call. There, the blue boxes indicate your input, gray are optional inputs while black and red are processing steps and output, respectively.
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

## Upload binding data (optional)

If you have conducted DELFIA or other techniques to receive binding data for certain sequences (usually sanger sequenced), you can upload these in a certain format and use these for clustering to potentially find other suitable sequences with high binding. The table has to have the following format and can be created in excel.
| aaSeqCDR3|Antigen 1|Antigen 2|Antigen 3|
|:----|:----|:----|:----|
|AIEAAAC|10000|30294|0|
|AEMNW|1000|0|0|
|PEICEES|0|1929|100000|

  
  You can upload the table as csv or xlsx file but make sure that the first column's name is *aaSeqCDR3* and *is in row 1*. Besides the given example, you can have a look at [an example file](src/ExpoSeq/software_tests/test_files/binding_data.csv). If you decide to work with that file make sure to delete the first column which contains the row number.

If you have prepared your data you can upload it with:
<br>

```python
plot.add_binding_data()
```
<br>

**Note**: If you decide to add more binding data to your analysis you can just use the same command and choose the new file with the filechooser and it will be added to the existing data. This can be also useful if you cannot manage to merge multiple antigens on the sequences in excel. Then you can just upload for each antigen separately the binding data.

## Data processing on a Cluster (optional)

First pull the folder with the scripts for the processing to your working directory

```bash
git clone https://github.com/nilshof01/ExpoSeq
```

I have prepared an example jobscript for working on a LSF cluster. You can have a look under
```bash
cd ExpoSeq/bash_processing
nano example_LSF_cluster.sh
```
To run your script interactively you can call:

```bash
python ~/ExpoSeq/bash_processing/mixcr_cl.py $PATH_TO_MIXCR $YOUR_EXPERIMENT_NAME $PATH_TO_FORWARD_FILES 
```

**NOTE**: You need to have installed MiXCR in your working directory to be able to start the processing.
To use multithreading and increase the RAM allocation have a look at the following parameter you can define:

- `--path_to_mixcr`: Is the filepath to the mixcr.jar file.
- `--experiment_name`: A string which is the name of your experiment.
- `--path_to_forward`: The directory of the fastq files with the forward reads. 
- `--path_to_backward`: The directory of the fastq files with the backward reads.
- `--threads`: The number of threads you would like to use for the processing.
- `--method`: The MiXCR method to align and assemble the reads you would like to use. Default is milab-human-tcr-dna-multiplex-cdr3
- `--java_heap_size`: Memory for proecssing in MB.
- `--add_arguments`: A string with further command line arguments, which are separated by a space. For instance, if you want to add the argument species, you would write `add_arguments "--species hsa"`. 
- `--mixcr_chain`: The chain you would like to extract from the data. 

**NOTE**: If you only want to process forward reads, then you do not need to add the path to the directory with the backward reads. Further, if you would like to analyze paired end sequencing data, please make sure that forward and backward fastq files are in separate folders.

After the processing has been finished you can import the folder with the processed files for the plotmanager. You can find the corresponding folder under
<br>
```
📦my_experiments
 ┣ 📂YOUR_EXPERIMENT_NAME
```
<br>
As soon as you call it you need to press 2 to upload the directory with the files.

## Talk with your data (optional)

I implemented [pandasai](https://docs.pandas-ai.com/en/latest/) in the pipeline to give you the option to investigate your data and even create plots without any further knowledge in programming with python. You can call the following command to start the chat with pandasai
<br>

```python
plot.chat()
```
The pipeline will prompt you to enter an API. This is necessary, since you need to connect the engine to a large language model. You can obtain the API key from OpenAI. 
<br>




## References
[1] Dmitriy A. Bolotin, Stanislav Poslavsky, Igor Mitrophanov, Mikhail Shugay, Ilgar Z. Mamedov, Ekaterina V. Putintseva, and Dmitriy M. Chudakov. "MiXCR: software for comprehensive adaptive immunity profiling." Nature methods 12, no. 5 (2015): 380-381.

[2] Dmitriy A. Bolotin, Stanislav Poslavsky, Alexey N. Davydov, Felix E. Frenkel, Lorenzo Fanchi, Olga I. Zolotareva, Saskia Hemmers, Ekaterina V. Putintseva, Anna S. Obraztsova, Mikhail Shugay, Ravshan I. Ataullakhanov, Alexander Y. Rudensky, Ton N. Schumacher & Dmitriy M. Chudakov. "Antigen receptor repertoire profiling from RNA-seq data." Nature Biotechnology 35, 908–911 (2017)

[3] (1, 2) Tareen A, Kinney JB (2019) Logomaker: beautiful sequence logos in Python. Bioinformatics btz921. bioRxiv doi:10.1101/635029.

[4] M.A. Larkin and others, Clustal W and Clustal X version 2.0, Bioinformatics, Volume 23, Issue 21, November 2007, Pages 2947–2948, https://doi.org/10.1093/bioinformatics/btm404

[5] Aric A. Hagberg, Daniel A. Schult and Pieter J. Swart, "Exploring network structure, dynamics, and function using NetworkX", in Proceedings of the 7th Python in Science Conference (SciPy2008), Gäel Varoquaux, Travis Vaught, and Jarrod Millman (Eds), (Pasadena, CA USA), pp. 11–15, Aug 2008


## Collaboration and Contact

Please feel encouraged to discuss on our discussions page on the repository. Further you can contact Nils Hofmann (nilhof01@gmail.com) if you have any questions or requests for future implementations. 

