from .collecting_all_arguments import ExpoSeqArgs,  prep_args
from ..plots.protein_embedding_umap import PrepareData # import as package
import pandas as pd

def call_args():
    Args = ExpoSeqArgs()
    Args.add_save_csv()
    Args.add_sequencing_report()
    Args.add_samples()
    Args.add_region_plots()
    Args.add_pca_components()
    Args.add_n_neighbors()
    Args.add_min_dist()
    Args.add_batch_size()
    Args.add_model_type()
    Args.add_strands()
    Args.add_binding_data()
    Args.add_antigen_names()
    Args.add_metric()
    Args.add_characteristic()
    Args.add_eps_dbscan()
    Args.add_min_pts_dbscan()
    Args.add_point_size()
    Args.add_n_jobs()
    
    return Args


args = call_args()
parser = prep_args(args)
# plot prepare
PrepData = PrepareData()
sequencing_report = pd.read_csv(parser.sequencing_report)
PrepData.tidy(sequencing_report, parser.samples, parser.region_plots, batch_size = parser.batch_size, 
              pca_components=parser.pca_components, perplexity=parser.n_neighbors, mind_dist=parser.min_dist,
              metric = parser.metric, characteristic = parser.characteristic, eps_dbscan = parser.eps, min_pts_dbscan = parser.min_pts, 
              n_jobs = parser.n_jobs, add_clone_size = parser.point_size,
              model_choice=parser.model_type, binding_data = parser.binding_data, antigens = parser.antigen_names)


PrepData.umap_results.to_csv(parser.save_csv)
