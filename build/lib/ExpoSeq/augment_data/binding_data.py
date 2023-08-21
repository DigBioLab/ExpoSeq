import pandas as pd
import os

try:
    import tkinter as tk
    from tkinter import filedialog
except:
    pass

def collect_binding_data(binding_data = None):
    if binding_data is None:
        binding_data = pd.DataFrame([])
    else:
        pass
    while True:
        # prompt the user to add a file
        print("You can either add an excel sheet or a csv file which contains the binding data. Note: the first column must contain the CDR3 sequences and its column name has to be aaSeqCDR3.")

        try:
            binding_file = filedialog.askopenfilename()
        except:
            while True:
                binding_file = input("copy and paste the path to your binding report")
                if os.path.isfile(os.path.abspath(binding_file)):
                    break
                else:
                    print("Please enter a valid filepath. ")
        while True:
            if binding_file.endswith(".xlsx") or binding_file.endswith(".csv"):
                if binding_file.endswith(".xlsx"):
                    binding_new = pd.read_excel(binding_file)
                elif binding_file.endswith(".csv"):
                    binding_new = pd.read_csv(binding_file)
                if binding_new.columns.to_list()[0] == "aaSeqCDR3":
                    break
                else:
                    print("Please change the header of the first column in your csv file to aaSeqCDR3")
            else:
                print("Please enter a valid filepath to a csv or xlsx file")


        binding_data = pd.concat([binding_data, binding_new])
        response = input("Do you want to continue adding files? (Y/n) ")
        if response.lower() == "n":
            break
        print("The first five rows of your binding data look like this:")
        print(binding_data.head(5))
    return binding_data