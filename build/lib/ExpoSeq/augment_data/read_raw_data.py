import pandas as pd

import os
try:
    import tkinter as tk
    from tkinter import filedialog
except:
    pass


def read_raw_data():
    try:
        filename = filedialog.askopenfilename()
    except:
        while True:
            filename = input("copy and paste the path to your file here")
            if os.path.isfile(filename) == True:
                break
            else:
                print("Please enter a valid filepath")
    filename = os.path.abspath(filename)
    pandas_raw_table = pd.read_table(filename)
    return pandas_raw_table