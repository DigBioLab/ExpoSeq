from src.ExpoSeq.pipeline import PlotManager
from time import time
import os


class IterativeAnalysis:
    def __init__(self, experiment):
        self.check_experiment(experiment)
        self.plot = PlotManager(experiment = experiment, allow_binding_data=False, test_version = True, no_automation=True )
        self.allowed_funcs = {"cluster_binding_data_umap": self.plot.cluster_binding_data_umap,
         "umap_clustering_characteristic": self.plot.umap_clustering_characteristic,
         "umap_sample_cluster": self.plot.umap_sample_cluster,
         }
        self.create_main_dir()
    
    @staticmethod
    def check_experiment(experiment):
        if not os.path.isdir(os.path.join("my_experiments", experiment)):
            raise Exception(f"The directory to {experiment} does not exist")
            
    
    @staticmethod
    def create_main_dir():
        if os.path.isdir("GeneratedSubplots") == False:
            os.mkdir("GeneratedSubplots")
        
    def iterate_and_visualize(self,iteration_name, func, attribute_name1:str, attribute_name2:str, attributes1:list, attributes2:list, **kwargs):
        """Generates all pngs for the the correspondi

        Args:
            func (_type_): _description_
            attribute_name1 (str): _description_
            attribute_name2 (str): _description_
            attributes1 (list): _description_
            attributes2 (list): _description_
        """
        dir = os.path.join("GeneratedSubplots", f"{attribute_name1}_{attribute_name2}", iteration_name)     
        all_captures = []
        if not os.path.isdir(dir):
            os.makedirs(dir)   
        n = 0
        for attr1 in attributes1:
            for  att2 in attributes2:
                params ={attribute_name1: attr1, attribute_name2: att2, **kwargs} 
                filename = os.path.join(dir, f"{n}_.png")
                if not os.path.isfile(filename):
                    try:
                        self.allowed_funcs[func](**params)
                    except:
                        print(f"Could not run {func} with {params}")
                    self.plot.add_to_subplot(figure_name = f"{attribute_name1}_{attribute_name2}", capture = f"{attr1}_{att2}", dir = dir)
               #     self.plot.ControlFigure.fig.savefig(filename, dpi = 300, format = "png")
                else: 
                    pass
                all_captures.append(f"{attr1}_{att2}")
                n += 1
        self.plot.show_subplot(figure_name=f"{attribute_name1}_{attribute_name2}", dir = dir, captures=all_captures)
    
    def performance(self, func, name:str, attribute_value:list):
        self.plot.is_test = True
        times = {}
        for attribute in attribute_value:
            start_time = time.time()
            method = self.allowed_funcs[func]
            method(name = attribute)
            end_time = time.time()
            times[attribute_value] = end_time - start_time
        return times

Iteration = IterativeAnalysis("max_both")
Iteration.iterate_and_visualize("test_characteristic",
                                "umap_sample_cluster", 
                                attribute_name1="n_neighbors",
                                attribute_name2="min_dist",
                                attributes1=[10,20,30],
                                attributes2 = [0.01,0.1, 0.4],
                                batch_size = 100,
                                show_strands = True,
                                characteristic = "length"
                                )