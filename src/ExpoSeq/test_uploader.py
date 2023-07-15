#from augment_data.uploader import upload
#from .augment_data.uploader import upload
from ExpoSeq.augment_data.uploader import upload
import os
import shutil
from ExpoSeq.full_pipe import run_pipeline
import pkg_resources

def run_initial_tests():
    pkg_dir = pkg_resources.resource_filename("ExpoSeq", "")
    module_dir = os.path.abspath("")
    if not os.path.isdir(os.path.join(module_dir, "my_experiments")):
        os.mkdir(os.path.join(module_dir, "my_experiments"))
    if not os.path.isdir(os.path.join(module_dir, "temp")):
        os.mkdir(os.path.join(module_dir, "temp"))
    if not os.path.isdir(os.path.join(module_dir, "my_experiments", "test_directory")):
        new_dir = os.path.join(module_dir, "my_experiments", "test_directory")
        os.mkdir(new_dir)
        source_seq_report = os.path.join(pkg_dir, "test_data","test_files", "sequencing_report.csv")
     #   source_align_report = os.path.join(pkg_dir, "test_data", "all_alignment_reports.pickle")
        source_experiments = os.path.join(pkg_dir, "test_data","test_files", "experiment_names.pickle")
        assert shutil.copy2(source_seq_report, new_dir), "could not copy sequencing report"
     #   shutil.copy2(source_align_report, new_dir)
        shutil.copy2(source_experiments, new_dir)
    else:
        pass
    
    #upload(testing=True, continue_analysis = "n", upload_type = "2", choose_exp = '1', paired_end_test = 'n', experiment_column = '1')
    assert upload(testing=True, continue_analysis = "n", upload_type = "1", choose_exp = '1', paired_end_test = 'n', experiment_column = '1'), "test 1 with upload_type 1,choose_exp = 1 and negative analysis continuation failed"
    print("test 1 passed")
    assert upload(testing=True, continue_analysis = "n", upload_type = "1", choose_exp = '2', paired_end_test = 'n', experiment_column = '1'), "test 2 with upload_type 1,choose_exp = 2 and negative analysis continuation failed"
    print("test 2 passed")
    assert upload(testing=True, continue_analysis = "n", upload_type = "3", choose_exp = '1', paired_end_test = 'n', experiment_column = '1'), "test 3 with upload_type 3,choose_exp = 1 and negative analysis continuation failed"	
    print("test 3 passed")
    assert upload(testing=True, continue_analysis = "Y", upload_type = "2", choose_exp = '1', paired_end_test = 'n', experiment_column = '1'), "test 4 with upload_type 2,choose_exp = 1 and positive analysis continuation failed"
    print("test 4 passed")
    
    print(os.path.abspath('run.py'))
    run_pipeline(test = True)
    print("all tests passed")

    
run_initial_tests()