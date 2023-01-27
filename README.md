# Welcome to ExpoSeq

ExpoSeq is a powerful pipeline for processing and analyzing FASTQ files. It utilizes [MiXCR](https://docs.milaboratories.com/mixcr/getting-started/installation/) to align and assemble from high-throughput sequencing data.

## Installation

To get started, please download and follow the instructions for MiXCR under the following link: https://docs.milaboratories.com/mixcr/getting-started/installation/

Once MiXCR is installed and set up, you can install the ExpoSeq repository by running the following command:
```pip install ExpoSeq```

## Importing the Plotting Tool

To access the plotting tool, you will need to import it into your console by running the following command:
```from expoSeq import PlotManager```
## Using the PlotManager

The PlotManager is the main interface for creating various plots using your FASTQ data. You can create an instance of the PlotManager by running the following command:
```plotting = PlotManager()```
To use the PlotManager to create plots, you will need to upload your FASTQ data to the pipeline. Please refer to the documentation for instructions on how to do this.

Once your data is uploaded, you can use the PlotManager to create a variety of plots, such as a Morosita-Horn plot. Here is an example of how to create this type of plot:
```plotting.morosita_horn()```


## For Collaborators
open a new project and pull the repository (link in browser) to your IDE. Then, install the libraries with pip install requirements.txt
