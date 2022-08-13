import pandas as pd
from tkinter import filedialog
import os

def read_raw_data():
    filename = filedialog.askopenfilename()
    filename = os.path.abspath(filename)
    pandas_raw_table = pd.read_table(filename)
    return pandas_raw_table