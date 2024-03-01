from typing import Any
from src.ExpoSeq.pipeline import PlotManager




class IterativeAnalysis:
    def __init__():
        self.plot = PlotManager(experiment = experiment)
        {"cluster_binding_data_umap": plot.cluster_binding_data_umap,
         "umap_clustering_characteristic": plot.umap_clustering_characteristic,
         "umap_sample_cluster": plot.umap_sample_cluster,
         }

plot = PlotManager()
neighbors = [10, 20, 30, 40, 50]
min_distances = [0.001, 0.01, 0.1, 0.5]
for n_neighbor, min_dist in zip(neighbors, min_distances):  
    plot.umap_sample_cluster(samples = ["Analysis1_R1", "Analysis2_R1"], n_neighbors=n_neighbor, min_dist=min_dist)
    plot.add_to_subplot(capture = f"{n_neighbor}_{min_dist}")