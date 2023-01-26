from os.path import dirname, abspath
from ast import literal_eval

def openParams(params_file):
    ROOT = dirname(abspath('ExpoSeq'))
    with open(ROOT + r"\python_scripts\plots\plot_params" + "\\" + params_file) as f:
        params = f.read()
    params = literal_eval(params)
    return params