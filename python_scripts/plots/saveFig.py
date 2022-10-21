from python_scripts.plots.plot_params.open_txtfiles import openParams
from matplotlib.pyplot import savefig

def saveFig():
    filename = input("how do you want to call the plot?")
    params = openParams("savePlot_params.txt")
    filepath = params['fname']
    format = params['format']
    path = filepath + '\\' + filename + '.' + format
    del params['fname']
    savefig(path, **params)