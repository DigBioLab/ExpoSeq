from src.ExpoSeq.pipeline import PlotManager
import os

def test_create_quarto():
    plot = PlotManager(experiment = "my_test_experiment", 
                       module_dir = "src/ExpoSeq/software_tests/test_files",
                       allow_binding_data=False)
    plot.create_report()
    html_file = os.path.join(plot.plot_path, plot.experiment + ".html")
    qmd_file = os.path.join(plot.plot_path, plot.experiment + ".qmd")
    assert os.path.isfile(html_file)
    assert os.path.isfile(qmd_file)
    with open(qmd_file, 'r', encoding='utf-8') as file:
        for line in file.readlines():
            assert "Cluster Antigens" not in line, "Cluster Antigens can't be a headline since no binding data was uploaded"
            
    plot = PlotManager(experiment = "my_test_experiment", 
                       module_dir = "src/ExpoSeq/software_tests/test_files",
                       allow_binding_data=r"C:\Users\nilsh\my_projects\ExpoSeq\src\ExpoSeq\software_tests\test_files\binding_data.csv")
    plot.create_report()
    html_file = os.path.join(plot.plot_path, plot.experiment + ".html")
    qmd_file = os.path.join(plot.plot_path, plot.experiment + ".qmd")
    assert os.path.isfile(html_file)
    assert os.path.isfile(qmd_file)
    with open(qmd_file, 'r', encoding='utf-8') as file:
        for line in file.readlines():
            if "Cluster Antigens" in line:
                binding_test = True
                break
            else:
                binding_test = False
                pass
    assert binding_test == True, "Cluster Antigens should be a headline since binding data was uploaded"


    