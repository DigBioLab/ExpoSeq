import pandas as pd
import os
from tabulate import tabulate

class ChooseMethod:
    def __init__(self) -> None:
        preset_path = os.path.join(os.getcwd(), "settings", "preset_list.csv")
        self.open_presets(preset_path)
        self.print_table()
        

        
        
    def open_presets(self, preset_path):
        
        if os.path.isfile(preset_path):
            self.preset_table = pd.read_csv(preset_path, sep=";")
        else:
            FileNotFoundError
            
    def print_table(self):
        print(tabulate(self.preset_table, headers='keys', tablefmt='psql'))      
    
    def ask_user(self):
        avail_preset = self.preset_table["Preset"].tolist()
        avail_preset.append("")
        while True:
            use_method = input(
                "Per default you will align your data with the following method: milab-human-tcr-dna-multiplex-cdr3.\nPress enter if you want to continue.\nOtherwise type in the method of your choice.\nIt has to be the exact same string which is given on the Mixcr documentation.")    
            if use_method in avail_preset:
                if use_method == "":
                    use_method = "milab-human-tcr-dna-multiplex-cdr3"
                else:
                    pass
                break
            else:
                print("The method you entered is not in the table. Please enter a valid method.")
        return use_method
    
    def check_args(self, chosen_method):
        args_table = self.preset_table[self.preset_table["Preset"] == chosen_method]
        args_table = args_table.iloc[:, 2:]
        add_arguments = {}
        method_vals = args_table.values.tolist()[0]
        gen_columns = args_table.columns.tolist()
        for vals in method_vals: 
            if vals != str(0) and vals != int(0):
                index = method_vals.index(vals)
                add_arguments[gen_columns[index]] = vals
            else:
                pass
        return add_arguments
    
    @staticmethod
    def change_args(add_arguments):
        if len(add_arguments) >= 1:
            for arg, value in add_arguments.items():
                usr_input = input(f"Do you want to enter another value for {arg}.\nThe default value is {value}. Press enter to continue, but if you want to change the input enter the corresponding value.")
                if usr_input == "":
                    continue
                else:
                    add_arguments[arg] = usr_input
        else:
            pass
        return add_arguments
    
    

class Test:
    def __init__(self) -> None:
        Test_Choose = ChooseMethod()
        assert Test_Choose.check_args("milab-human-tcr-dna-multiplex-cdr3") == {}, "output must be an empty dict"
        assert Test_Choose.check_args("ont-rna-seq-vdj-full-length") == {"species": "hsa"}, "output must contain the hsa value with species key"
        assert Test_Choose.change_args({}) == {}, "No prompting should appear"
        Test_Choose.open_presets("") == FileNotFoundError, "An error should appear"
    
