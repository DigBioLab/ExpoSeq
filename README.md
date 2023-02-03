# Welcome to ExpoSeq

ExpoSeq is a powerful pipeline for processing and analyzing FASTQ files from sequencing phage Display panning samples. By using the pipeline you can receive your plots in only 2 steps. The pipeline focuses on analysing the identity between samples but also applies various clustering techniques to analyse the relation between the sequences. Furthermore you can add binding data to relate the clusters to affinity. It utilizes [MiXCR](https://docs.milaboratories.com/mixcr/getting-started/installation/) to align and assemble from high-throughput sequencing data.

## Installation

open a new project and pull the repository (link in browser) to your IDE. Then, install the libraries with ```pip install -r requirements.txt```

To get started, please download and follow the instructions for MiXCR under the following link: https://docs.milaboratories.com/mixcr/getting-started/installation/ 
You can also only use the test version without installing it.

## Importing the Plotting Tool

To access the plotting tool, you will need to import it into your console by running the following command:
```from expoSeq import PlotManager```
## Using the PlotManager

The PlotManager is the main interface for creating various plots using your FASTQ data. You can create an instance of the PlotManager by running the following command:
```plot = PlotManager()```
To use the PlotManager to create plots, you will need to upload your FASTQ data to the pipeline. This will automatically happen as soon as you start call the PlotManager. In the following you can obtain an insight in the worklow which follows after you have called the class. There the blue boxes indicate your input, gray are optional inputs while black and red are processing steps and output, respectively.
![relative_path_to_image](workflow_ExpoSeq.png)
If you just want to test the pipeline and see its functions you can call: ```plotting = PlotManager(test_version = True)```

Once you have called the test version or have finished the data processing, you can use the PlotManager to create a variety of plots, such as a Morosita-Horn plot. Here is an example of how to create this type of plot:
```plot.morosita_horn()```
If you want to change the style of the plot you can use the class. If you have called it plotting you can do for instance the following: ```plotting.style.title_xaxis("your_title")``` 
If you want to implement further plot change you can also refer to the matplotlib.pyplot library and change it in the same way as following:
```import matplotlib.pyplot as plt```
```plt.xlabel("your_title")```
If you need further help you can also ask the chatpot which is trained on gpt3. Therefore call: ```plotting.askMe()``` . Yet it is not trained on the package but you can use it to customize your graphs. 
If you would like to have details about the input and function of the plot functions call for instance: ```help(plotting.usqPlot)



## For Collaborators



