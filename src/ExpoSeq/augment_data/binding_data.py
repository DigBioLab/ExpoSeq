import pandas as pd
import os

try:
    import tkinter as tk
    from tkinter import filedialog
except:
    pass

def collect_binding_data(binding_data = None):
    if binding_data == None:
        binding_data = pd.DataFrame([])
    else:
        pass
    while True:
        # prompt the user to add a file
        print("add your excel sheet with the binding data with the file chooser")

        try:
            binding_file = filedialog.askopenfilename()
        except:
            while True:
                binding_file = input("copy and paste the path to your binding report")
                if os.path.isfile(os.path.abspath(binding_file)):
                    break
                else:
                    print("Please enter a valid filepath. ")
        if binding_file.endswith(".xlsx"):
            binding_new = pd.read_excel(binding_file)
        elif binding_file.endswith(".csv"):
            binding_new = pd.read_csv(binding_file)
        binding_data = pd.concat([binding_data, binding_new])
        response = input("Do you want to continue adding files? (Y/n) ")
        if response.lower() == "n":
            break
    return binding_data