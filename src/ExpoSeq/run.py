from ExpoSeq.full_pipe import run_pipeline
from ExpoSeq.pipeline import PlotManager

def run_automatic(test = False, plotmanager = PlotManager()):
    if not test:
        run_pipeline(plot = plotmanager, test = test)
    else:
        run_pipeline(plot = plotmanager, test = True)