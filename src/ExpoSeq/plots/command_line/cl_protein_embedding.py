from .collecting_all_arguments import ExpoSeqArgs
from ExpoSeq.plots.protein_embedding import PrepareData # import as package

def call_args():
    Args = ExpoSeqArgs()
    Args.save_csv()
    Args.sequencing_report()
    Args.samples()
    Args.region_of_interest()
    Args.pca_components()
    Args.perplexity()
    Args.iterations_tsne()
    Args.batch_size()
    Args.model_type()
    Args.strands()
    return Args

Args = call_args()
Args.parser.parse_args()
PrepData = PrepareData()
PrepData.tidy(Args.sequencing_report, Args.samples, Args.region_of_interest, batch_size = Args.batch_size, 
              pca_components=Args.pca_components, perplexity=Args.perplexity, iterations_tsne=Args.iterations_tsne,
              model_choice=Args.model_type)
PrepData.make_csv(Args.save_csv)
