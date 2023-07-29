from .pipeline import PlotManager
import os
import pandas as pd
import shutil
import pkg_resources
from ast import literal_eval

def run_pipeline(antigens = None, usq_multiple = "n", pass_interrupt = True, test = False):
    pkg_path = pkg_resources.resource_filename("ExpoSeq", "")
    save_settings_path = os.path.join(pkg_path,
                 "settings",
                 "save_settings.txt")
    with open(save_settings_path, "r") as f:
        save_settings = f.read()
    save_settings = literal_eval(save_settings)
    original_path = save_settings["fname"]
    if antigens == True or antigens == False or antigens == None:
        pass
    else:
        raise ValueError("Enter the antigens in a list and seperate them by comma")
    if usq_multiple in ["Y", "y", "n", "N"]:
        pass
    else:
        raise ValueError("Invalid value. Please enter Y or n")
    if pass_interrupt in [True, False]:
        pass
    else: 
        raise ValueError("Invalid value. Please enter True or False")

    
    if not os.path.isdir("results"):
        os.mkdir("results")
    else:
        pass
    result_dir = os.path.abspath("results")

    plot = PlotManager(test_version = test)
    folder_name = plot.experiment
    experiment_name = folder_name
    if not test or not os.path.isdir(os.path.join(result_dir, "Test")):
        os.mkdir(os.path.join(result_dir, folder_name))
        os.mkdir(os.path.join(result_dir, folder_name, "Quality"))
    assert os.path.isdir(os.path.join(result_dir,
                                      folder_name,
                                      "Quality")), "The directory could not be created"
    save_path = os.path.abspath(os.path.join(result_dir,
                                             folder_name,
                                             "Quality"))
    plot.settings_saver.change_save_path(path = save_path), "The save path could not be changed"

    plot.lengthDistribution_multi()
    if test:
        print("length Distribution successfully created")
    plot.save(name = "length_Distribution" + folder_name)
    plot.logoPlot_multi()
    if test:
        print("Logo Plot successfully created")
    plot.save(name = "logo Plots " + folder_name)
    usq_multiple = usq_multiple
    unique_experiments = list(plot.unique_experiments.values())
    if usq_multiple.lower() in ["Y", "y"]:
        more_plots = True
        while more_plots == True:
            plot.print_samples()
            sample_names = input("Choose the sample from the list above which you want to include. Please seperate your input with commas")
            
            if not type(sample_names) == list:
                sample_names = sample_names.split(",")
                sample_names = [value.replace("'", "").replace('"', '') for value in sample_names]
                sample_names = [value.replace(" ", "") for value in sample_names]
                sample_names = [value for value in sample_names if value]
            else:
                pass
            for value in sample_names:
                if value not in unique_experiments and pass_interrupt == False:
                    interrupt = input("The system could not find " + value + " in your sample name. If you want to try again, press Y. Otherwise the analysis will be conducted without this value")
                    while True:
                        if interrupt.lower() in ["Y", "y"]:
                            break
                        else:
                            sample_names.remove(value)
                
            plot.usqPlot(samples = sample_names)
            while True:
                usq_name = input("How do you want to call your USQ Plot?")
                if not os.path.isfile(os.path.join(save_path, usq_name)):
                    plot.save(usq_name)
                    break
                else:
                    print("The file already exists. Please try another value.")
            while True:
                more_plots = input("Do you want to create more USQ Plots with different samples? (Y/n)")
                if more_plots.lower() in ["Y", "y", "n", "N"]:
                    break
                else:
                    print("Please enter a valid value")
            if more_plots.lower() in ["Y", "y"]:
                more_plots = True
                unique_experiments = [value for value in unique_experiments if value not in sample_names]
            else:
                more_plots = False
                break
    else:
        for experiment in unique_experiments:
            try:
                plot.usqPlot(samples = [experiment])

                plot.save("Usq Plot " + experiment)
            except:
                print("A problem appeard for " + experiment + ". The USQ plot could not be created")
    unique_experiments = list(plot.unique_experiments.values())
    plot.morosita_horn()
    if test:
        print("Morosita Horn successfully created")
    plot.save("morosita_horn_identity")
    plot.jaccard()
    if test:
        print("Jaccard successfully created")
    plot.save("jaccard_identity")
    plot.sorensen()
    if test:
        print("Sorensen Dice successfully created")
    plot.save("Sorensen_Dice_identity")
    plot.relative()
    if test:
        print("Relative Identity successfully created")
    plot.save("Relative_Identity")

    os.mkdir(os.path.join(result_dir,
                          folder_name,
                          "Sequence Analysis"))
    assert os.path.isdir(os.path.join(result_dir,folder_name, "Sequence Analysis")), "The directory : " + os.path.join(result_dir,folder_name, "Sequence Analysis") + " could not be created"
    for experiment in unique_experiments:
        try:
            os.mkdir(os.path.join(result_dir,folder_name, "Sequence Analysis", experiment))
            save_path = os.path.abspath(os.path.join("results", folder_name, "Sequence Analysis", experiment))
            plot.settings_saver.change_save_path(path = save_path)
            plot.logoPlot_single(experiment)
            print("Logo Plot for " + experiment + " created")
            plot.save("Logo Plot for " + experiment)
            plot.rel_seq_abundance(samples = [experiment])
            print("Relative Sequence Abundance for " + experiment + " created")
            plot.save("Sequence Abundance for " + experiment)
     #      plot.basic_cluster(experiment)
      #      plot.save("Clustering based on Levenshtein Distance")
            try:
                plot.aa_distribution(experiment, region = [3, 8])
                plot.save("AA Distribution " + "in region 3 - 8")
                print("AA Distribution for " + experiment + " created")
            except:
                print("sequence logo plot for " + experiment + "in region: "+ "3, 8" + " could not be created")
            try:
                plot.aa_distribution(sample = experiment, region = [9, 14])
                
                plot.save("AA Distribution "+ "in region 9 - 14")
            except: 
                print("sequence logo plot for " + experiment + "in region: "+ "9, 14" + " could not be created")
            try:
                plot.aa_distribution(sample = experiment, region = [15, 20])
                plot.save("AA Distribution " + "in region 15 - 20")
            except:
                print("sequence logo plot for " + experiment + "in region: "+ "15, 20" + " could not be created")

        except: pass

    if isinstance(plot.binding_data, pd.DataFrame) == True:
        print("You have binding data for the following Antigens")
        plot.print_antigens()
        all_antigens = plot.binding_data.columns.to_list()[1:-1]
        if test != True and antigens == None:
            count_AG = input("If you want to create a TSNE Embedding for all Antigens press 1. If you want to select certain Antigens, press 2.")
            while True:
                if count_AG in ["1", "2"]:
                    break
                else:
                    print("Please enter a valid value")
        else:
            count_AG = "1"
        if count_AG == str(1):
            antigens = all_antigens
        else:
           
            if not type(antigens) == list:
                antigens = antigens.split(",")
                antigens = [value.replace("'", "").replace('"', '') for value in antigens]
                antigens = [value.replace(" ", "") for value in antigens]
                assert type(antigens) == list, "The antigens have to be in a list. The program could not solve this error"
            else:
                pass
            for value in antigens:
                if value not in all_antigens:
                    while True:
                        interrupt = input("The system could not find " + value + " in your sample name. If you want to try again, press Y. Otherwise press n and the analysis will be conducted without this value")
                        if interrupt in ["Y", "y", "n", "N"]:
                            break
                        else:
                            print("Please enter a valid value")
                    if interrupt.lower() in ["Y", "y"] and not test:
                        break
                    else:
                        antigens.remove(value)
                            
        for antigen in antigens:
            os.mkdir(os.path.join(result_dir, folder_name, "Binding Analysis"))
            os.mkdir(os.path.join(result_dir, folder_name, "Binding Analysis", antigen))
            save_path = os.path.abspath(os.path.join("results", folder_name, "Binding Analysis", antigen))
            for experiment in unique_experiments:
                plot.tsne_cluster_AG(sample = experiment, toxins = [antigen], toxin_names = False)
                if test:
                    break

    else:
        print("You have not given any binding data. So this part of the analysis has not been carried out.")
    if test:
        shutil.rmtree(os.path.join(result_dir, folder_name))

    plot.settings_saver.change_save_path(path = original_path)
    
    
    
