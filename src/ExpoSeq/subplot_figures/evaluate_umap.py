from src.ExpoSeq.pipeline import PlotManager
from time import time
import os


class IterativeAnalysis:
    def __init__(self, experiment):
        """
         Initialize the plot manager. This is called by __init__ and should not be called directly. It is the responsibility of the user to call this method if they want to do something other than setup the experiment before it is passed to the plot manager.
         
         @param experiment - The experiment to use for plotting. If None a default experiment will be used
        """
        self.check_experiment(experiment)
        self.plot = PlotManager(experiment = experiment, allow_binding_data=False, test_version = True, no_automation=True )
        self.allowed_funcs = {"cluster_binding_data_umap": self.plot.cluster_binding_data_umap,
         "umap_clustering_characteristic": self.plot.umap_clustering_characteristic,
         "umap_sample_cluster": self.plot.umap_sample_cluster,
         }
        self.create_main_dir()
    
    @staticmethod
    def check_experiment(experiment):
        """
         Check if the experiment exists. This is a helper function to be used by other functions in this module
         
         @param experiment - Name of the experiment to
        """
        # If the directory to experiment does not exist raise an exception
        if not os.path.isdir(os.path.join("my_experiments", experiment)):
            raise Exception(f"The directory to {experiment} does not exist")
            
    
    @staticmethod
    def create_main_dir():
        """
         Create the directory to store plots. This is called by plot_to_pdf and plot_to
        """
        # Create the generated subplots directory if it doesn t exist.
        if os.path.isdir("GeneratedSubplots") == False:
            os.mkdir("GeneratedSubplots")
        
    def iterate_and_visualize(self,iteration_name, func, attribute_name1:str, attribute_name2:str, attributes1:list, attributes2:list, **kwargs):
        """
         Iterates over the attributes and visits them. This is a wrapper for the pymel. visualize method.
         
         @param iteration_name - Name of the iteration to use for the visualization
         @param func - Function to be called for each attribute
         @param attribute_name1 - Name of the first attribute ( str )
         @param attribute_name2 - Name of the second attribute ( str )
         @param attributes1 - List of attribute names that should be plotted in the first iteration ( list ). If this is a list it is assumed that the first attribute is a list of strings.
         @param attributes2 - List of attribute names that should be plotted in the second iteration ( list
        """

        dir = os.path.join("GeneratedSubplots", f"{attribute_name1}_{attribute_name2}", iteration_name)     
        all_captures = []
        # Create a directory if it doesn t exist.
        if not os.path.isdir(dir):
            os.makedirs(dir)   
        n = 0
        # This function is called by the plot.
        for attr1 in attributes1:
            # This function is called by the plot.
            for  att2 in attributes2:
                params ={attribute_name1: attr1, attribute_name2: att2, **kwargs} 
                filename = os.path.join(dir, f"{n}_.png")
                # if filename is not a file or if filename is not a file then the function will be called with the same parameters as the filename.
                if not os.path.isfile(filename):
                    self.allowed_funcs[func](**params)
                  #  print(f"Could not run {func} with {params}")
                    self.plot.add_to_subplot(figure_name = f"{attribute_name1}_{attribute_name2}", capture = f"{attr1}_{att2}", dir = dir)
               #     self.plot.ControlFigure.fig.savefig(filename, dpi = 300, format = "png")
                else: 
                    pass
                all_captures.append(f"{attr1}_{att2}")
                n += 1
        self.plot.show_subplot(figure_name=f"{attribute_name1}_{attribute_name2}", dir = dir, captures=all_captures)
    
    def performance(self, func, name:str, attribute_value:list):
        """
         Measure the performance of a function. This is a wrapper around : meth : ` allowed_funcs ` to allow custom functions to be used in the plot.
         
         @param func - The name of the function to measure. It must be one of the allowed functions in self. allowed_funcs.
         @param name - The name of the function that is being measured.
         @param attribute_value - The list of attribute values that are being measured.
         
         @return A dictionary of time in seconds for each attribute value in the list. The key is the attribute value
        """
        self.plot.is_test = True
        times = {}
        # This method is used to calculate the time of each attribute value.
        for attribute in attribute_value:
            start_time = time.time()
            method = self.allowed_funcs[func]
            method(name = attribute)
            end_time = time.time()
            times[attribute_value] = end_time - start_time
        return times

import pandas as pd
Iteration = IterativeAnalysis("max_both")
binding_data = pd.read_csv(r"C:\Users\nilsh\OneDrive\Desktop\DTU\NGS_pipeline\data\Binding_data\Chris_main_df.csv")
#pla2 = binding_data.columns.tolist().index("Ecarpholin_S_DDEL_5ug/mL")
Iteration.plot.binding_data = binding_data
Iteration.plot.preferred_antigen = ["Ecarpholin_S_DDEL_5ug/mL"]
Iteration.iterate_and_visualize("binding_ecarpholin",
                                "cluster_binding_data_umap", 
                                attribute_name1="n_neighbors",
                                attribute_name2="min_dist",
                                attributes1=[10,20,30, 40],
                                attributes2 = [0.01,0.1, 0.4],
                                samples = ["Analysis1_R1", "Analysis2_R1", "Analysis3_R1", "Analysis4_R1"],
                                batch_size = 200,
                                show_strands = True,
                                antigens = ["Ecarpholin_S_DDEL_5ug/mL"],
                                )