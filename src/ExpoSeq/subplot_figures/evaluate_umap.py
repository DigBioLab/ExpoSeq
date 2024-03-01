from src.ExpoSeq.pipeline import PlotManager
from time import time



class IterativeAnalysis:
    def __init__(self):
        self.plot = PlotManager(experiment = experiment)
        {"cluster_binding_data_umap": self.plot.cluster_binding_data_umap,
         "umap_clustering_characteristic": self.plot.umap_clustering_characteristic,
         "umap_sample_cluster": self.plot.umap_sample_cluster,
         }
        
    def iterate_and_visualize(self, func):
        getattr(func)()
        self.plot.add_as_subplot()
    
    
    def performance(self, func, name:str, attribute_value:list):
        times = []
        for attribute in attribute_value:
            start_time = time.time()
            method = getattr(func)
            method(name = attribute)
            end_time = time.time()
            times.append((end_time - start_time))
        return times

plot = PlotManager()
neighbors = [10, 20, 30, 40, 50]
min_distances = [0.001, 0.01, 0.1, 0.5]
for n_neighbor, min_dist in zip(neighbors, min_distances):  
    plot.umap_sample_cluster(samples = ["Analysis1_R1", "Analysis2_R1"], n_neighbors=n_neighbor, min_dist=min_dist)
    plot.add_to_subplot(capture = f"{n_neighbor}_{min_dist}")