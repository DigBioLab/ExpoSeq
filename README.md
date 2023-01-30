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
To use the PlotManager to create plots, you will need to upload your FASTQ data to the pipeline. This will automatically happen as soon as you start call the PlotManager. If you just want to test the pipeline and see its function you can call: ```plotting = PlotManager(test_version = True)```

Once your data is uploaded, you can use the PlotManager to create a variety of plots, such as a Morosita-Horn plot. Here is an example of how to create this type of plot:
```plotting.morosita_horn()```
If you want to change the style of the plot you can use the class. If you have called it plotting you can do for instance the following: ```plotting.style.title_xaxis("your_title")``` 
If you want to implement further plot change you can also refer to the matplotlib.pyplot library and change it in the same way as following:
```import matplotlib.pyplot as plt```
plt.xlabel("your_title")



## For Collaborators
open a new project and pull the repository (link in browser) to your IDE. Then, install the libraries with pip install -r requirements.txt
