import pandas as pd
from tkinter import filedialog

def collect_binding_data():
    binding_data = pd.DataFrame([])
    while True:
        # prompt the user to add a file
        print("add your excel sheet with the binding data with the file chooser")
        binding_file = filedialog.askopenfilename()
        binding_data = pd.concat([binding_data, binding_file])
        response = input("Do you want to continue adding files? (Y/n) ")
        if response.lower() == "n":
            break
    return binding_data