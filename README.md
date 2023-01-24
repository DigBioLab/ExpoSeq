# ExpoSeq
Installation
1. For processing your fastq files correctly in the pipeline download and follow the instructions for mixcr under the following link.
  https://docs.milaboratories.com/mixcr/getting-started/installation/
2. install the repository with 
  /code pip install ExpoSeq
3. import the plotting tool in your console with:
  from expoSeq import PlotManager
4. call the plotmanager with for instance:
    plotting = PlotManager()
5. use plotting to make various plot for instance:
   plotting.morosita_horn()
   for that you have to follow some steps in the beginning where you upload your data
6. 