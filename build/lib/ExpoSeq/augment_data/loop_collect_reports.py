import pandas as pd
import csv

def load_mixed_files(filenames, testing, value = None):
    all_alignment_reports = pd.DataFrame([])
    sequencing_report = pd.DataFrame([])
    if not testing:
        column_choice = input("Please check if your report contains a column which contains the name of the samples. If the name of the column is Experiment you can skip this part by pressing 1. If it has another column name you can enter it here. If your files do not contain the column you have to enter the name of the experiment for each file individually. Therefore, press 3 to continue.")
        while True:
            if column_choice in ["1", '2', '3']:
                break
            else:
                print("Please enter a valid value")
                column_choice = input()
    else:
        column_choice = value
    for i in filenames:
        if "AlignmentReport" in i:
            report = pd.read_table(i)
            assert report.shape[1] == 1, "The alignment report should only one column. Please check your file."
            splitted_report = report.iloc[:, 0].str.split(":",
                                                          expand=True)
            splitted_report = splitted_report.iloc[:, 0:2].T
            transposed_report = splitted_report.rename(columns=splitted_report.iloc[0]).drop(splitted_report.index[0])
            all_alignment_reports = pd.concat([all_alignment_reports, transposed_report])
        else:
            try:
                with open(i, newline='') as f:
                    reader = csv.reader(f)
                    # The csv.reader object can read the file line by line
                    # If the file is not a CSV file, this will raise a csv.Error exception
                    # Otherwise, it will read the file successfully
                    for row in reader:
                        pass
                is_csv = True
            except:
                is_csv = False
            try:
                if is_csv == True:

                    local_intermediate_report = pd.read_csv(i)
                else:
                    local_intermediate_report = pd.read_table(i)

                while True:
                    if column_choice == "1":
                        break
                    elif column_choice == "3":
                        sample_name = input("Please enter the name of the sample for file " + i)
                        local_intermediate_report["Experiment"] = sample_name
                        break
                    else:
                        if column_choice in local_intermediate_report.columns:
                            local_intermediate_report = local_intermediate_report.rename(columns={column_choice: 'Experiment'})
                            break
                        else:
                            if not testing:
                                while True:
                                    column_choice = input("Apparently, you have entered the wrong column. Please enter the correct column name. If the name of the column is Experiment you can skip this part by pressing 1.")
                                    if column_choice == "1" or column_choice in local_intermediate_report.columns:
                                        break
                                    else:
                                        print("Please enter a valid value")
                            else: 
                                pass

                sequencing_report = pd.concat([sequencing_report, local_intermediate_report])
            except:
                print(i + " could not be read.")
    try:
        del all_alignment_reports["index"]
    except:
        pass
    return sequencing_report, all_alignment_reports





def load_alignment_reports(filenames):
    all_alignment_reports = pd.DataFrame()
    for i in filenames:
        if "AlignmentReport" in i:
            report = pd.read_table(i)
            splitted_report = report.iloc[:, 0].str.split(":",
                                                          expand = True)
            splitted_report = splitted_report.iloc[:, 0:2].T
            transposed_report = splitted_report.rename(columns=splitted_report.iloc[0]).drop(splitted_report.index[0])
            all_alignment_reports = pd.concat([all_alignment_reports, transposed_report])
    all_alignment_reports = all_alignment_reports.reset_index()
    try:
        del all_alignment_reports["index"]
    except:
        pass
    return all_alignment_reports