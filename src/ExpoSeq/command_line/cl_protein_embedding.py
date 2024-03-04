from .collecting_all_arguments import ExpoSeqArgs, prep_args
from ..plots.protein_embedding import PrepareData  # import as package
import pandas as pd


def call_args():
    Args = ExpoSeqArgs()
    Args.add_save_csv()
    Args.add_sequencing_report()
    Args.add_samples()
    Args.add_region_plots()
    Args.add_pca_components()
    Args.add_perplexity()
    Args.add_iterations_tsne()
    Args.add_batch_size()
    Args.add_eps_dbscan()
    Args.add_min_pts_dbscan()
    Args.add_model_type()
    Args.add_strands()
    Args.add_binding_data()
    Args.add_antigen_names()
    Args.add_point_size()
    Args.add_embedding_vector_path()
    return Args


args = call_args()
parser = prep_args(args)
# plot prepare
PrepData = PrepareData()
sequencing_report = pd.read_csv(parser.sequencing_report)
if parser.binding_data is not None:
    binding_data = pd.read_csv(parser.binding_data)
else:
    binding_data = None
PrepData.tidy(
    sequencing_report,
    parser.samples,
    parser.region_plots,
    batch_size=parser.batch_size,
    pca_components=parser.pca_components,
    perplexity=parser.perplexity,
    iterations_tsne=parser.iterations_tsne,
    model_choice=parser.model_type,
    binding_data=binding_data,
    antigens=parser.antigen_names,
    eps_dbscan=parser.eps,
    min_pts_dbscan=parser.min_pts,
    X_path=parser.embedding_vector_path,
)


PrepData.tsne_results.to_csv(parser.save_csv)
