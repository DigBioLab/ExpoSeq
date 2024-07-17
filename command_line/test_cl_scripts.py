import subprocess
import os
import unittest
import pandas as pd

class TestCLScripts(unittest.TestCase):
    def setUp(self):
        self.sequencing_report_path = os.path.join(os.path.dirname(__file__),"data", "example_report.csv")
        print(self.sequencing_report_path)
        self.scripts_dir = os.path.dirname(__file__)
        self.temp_dir = os.path.join(self.scripts_dir, "temp_dir")
        self.sequencing_report = self.prep_seq_report()
        self.seq_report_fr3 = self.prep_seq_report(region = "FR3")
        self.region_plots = "aaSeqCDR3"
        self.sample_name = "Sample1"
        self.second_sample_name = "Sample2"
        self.alternative_region = "aaSeqFR3"
        self.binding_data = self.sample_binding_data()
        
    def prep_seq_report(self, region = "CDR3",):
        sequencing_report = os.path.join(self.temp_dir, f"sequencing_report_{region}.csv")
        result = subprocess.run(["python", os.path.join(self.scripts_dir, "cl_sequencing_report.py"), 
                                 "--tsv_dir", self.sequencing_report_path,
                                 "--region", region, "--save_csv", sequencing_report])
        self.assertEqual(result.returncode, 0)
        return sequencing_report
        
    @staticmethod
    def remove_file(file_path):
        if os.path.exists(file_path):
            os.remove(file_path)
    
    def sample_binding_data(self):
        """Returns sampled binding data. So it is not real. Contains binding data for one virtual antigen called Antigen 1
        """
        result = subprocess.run(["python", os.path.join(self.scripts_dir, "generate_binding_data.py"), "--sequencing_report", self.sequencing_report,
                                 "--save_csv", os.path.join(self.temp_dir, "binding_data.csv"), "--samples", self.sample_name], capture_output=True, text=True)
        return os.path.join(self.temp_dir, "binding_data.csv")
    
    def test_cd_cf_distribution(self):
        csv_save = os.path.join(self.scripts_dir, "temp_dir", "cd_cf_distribution.csv")
        result = subprocess.run(['python', os.path.join(self.scripts_dir, "cd_cf_distribution.py"), "--save_csv", 
                                csv_save, "-r", self.sequencing_report, "--single_sample", self.sample_name, ], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0)
        self.remove_file(csv_save)

    def test_cd_cf_length_distribution(self):
        csv_save = os.path.join(self.scripts_dir, "temp_dir", "cf_length_distribution.csv")
        result = subprocess.run(["python", os.path.join(self.scripts_dir, "cl_cf_length_distribution.py"), "--save_csv", csv_save, "--sequencing_report", 
                                self.sequencing_report, "--single_sample", self.sample_name, "--region_plots", self.region_plots], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0)
        self.remove_file(csv_save)
        
    def test_test_diversity(self):
        csv_save = os.path.join(self.scripts_dir, "temp_dir", "diversity.csv")
        result = subprocess.run(["python", os.path.join(self.scripts_dir, "cl_diversity_plot.py"), "--save_csv", csv_save, "--sequencing_report", self.sequencing_report, "--region_plots", self.region_plots], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0)
        self.remove_file(csv_save)
        seq_report = self.prep_seq_report(region = "FR3")
        result = subprocess.run(["python", os.path.join(self.scripts_dir, "cl_diversity_plot.py"), "--save_csv", csv_save, "--sequencing_report", self.seq_report_fr3, "--region_plots", self.alternative_region], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, "Alternative region failed")
        self.remove_file(csv_save)
        
    def test_cl_length_distribution(self):
        csv_save = os.path.join(self.scripts_dir, "temp_dir", "cl_length_distribution.csv")
        result = subprocess.run(["python", os.path.join(self.scripts_dir, "cl_length_distribution.py"),
                                 "--save_csv", csv_save, "--sequencing_report", self.sequencing_report, "--single_sample", self.sample_name, "--region_plots", self.region_plots], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0)
        self.remove_file(csv_save)
        
        result = subprocess.run(["python", os.path.join(self.scripts_dir, "cl_length_distribution.py"),
                                 "--save_csv", csv_save, "--sequencing_report", self.seq_report_fr3, "--single_sample", self.sample_name,
                                 "--region_plots", self.alternative_region], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, "Alternative region failed")
        self.remove_file(csv_save)
        
    def test_cl_levenshtein_clustering(self):
        csv_save = os.path.join(self.scripts_dir, "temp_dir", "cl_levenshtein_clustering.csv")
        result = subprocess.run(["python", os.path.join(self.scripts_dir, "cl_levenshtein_clustering.py"), 
                                 "--save_csv", csv_save, "--sequencing_report", self.sequencing_report, 
                                 "-s", f"{self.sample_name}", "--batch_size", "100", 
                                 "--region_plots", self.region_plots, "--binding_data", self.binding_data, 
                                 "--antigen_names", "Antigen 1"], capture_output=True, text=True)
        
        self.assertEqual(result.returncode, 0)
        self.remove_file(csv_save)
        
    def test_logoplot(self):
        csv_save = os.path.join(self.scripts_dir, "temp_dir", "logoplot.csv")
        result = subprocess.run(["python", os.path.join(self.scripts_dir, "cl_logoplot.py"), "--save_csv", csv_save, "--sequencing_report", self.sequencing_report, 
                                 "--samples", self.sample_name, self.second_sample_name, "--chosen_seq_length", "15", "--region_plots", self.region_plots, "--method_logo", "bits"], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0)
        self.remove_file(csv_save)
        
    def test_cl_lvst_histogram(self):
        csv_save = os.path.join(self.scripts_dir, "temp_dir", "cl_lvst_histogram.csv")
        result = subprocess.run(["python", os.path.join(self.scripts_dir, "cl_lvst_histogram.py"), "--save_csv", csv_save, "--sequencing_report", self.sequencing_report, 
                                 "-s", self.sample_name, self.second_sample_name, "--region_plots", self.region_plots, "--batch_size", "100"], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0)
        self.remove_file(csv_save)
        result = subprocess.run(["python", os.path.join(self.scripts_dir, "cl_lvst_histogram.py"), "--save_csv", csv_save, "--sequencing_report", self.seq_report_fr3, 
                                 "-s", self.sample_name, self.second_sample_name, "--region_plots", self.alternative_region, "--batch_size", "100"], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, "Alternative region failed")
        self.remove_file(csv_save)
        
    def test_cl_matrix_identity(self):
        csv_save = os.path.join(self.scripts_dir, "temp_dir", "cl_matrix_identity.csv")
        result = subprocess.run(["python", os.path.join(self.scripts_dir, "cl_matrix_identity.py"), "--save_csv", csv_save,
                                 "--sequencing_report", self.sequencing_report, 
                                 "--region_plots", self.region_plots, "--matrix_type", "morosita_horn"], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0)
        self.remove_file(csv_save)
        result = subprocess.run(["python", os.path.join(self.scripts_dir, "cl_matrix_identity.py"), "--save_csv", csv_save,
                                 "--sequencing_report", self.sequencing_report, 
                                 "--region_plots", self.region_plots, "--matrix_type", "sorensen"], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0)
        self.remove_file(csv_save)
        result = subprocess.run(["python", os.path.join(self.scripts_dir, "cl_matrix_identity.py"), "--save_csv", csv_save,
                                 "--sequencing_report", self.sequencing_report, 
                                 "--samples", self.sample_name, self.second_sample_name,
                                 "--region_plots", self.region_plots, "--matrix_type", "jaccard"], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0)
        self.remove_file(csv_save)
        result = subprocess.run(["python", os.path.join(self.scripts_dir, "cl_matrix_identity.py"), "--save_csv", csv_save,
                                 "--sequencing_report", self.seq_report_fr3, 
                                 "--samples", self.sample_name, self.second_sample_name,
                                 "--region_plots", self.alternative_region, "--matrix_type", "jaccard"], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, "Alternative region failed")
        self.remove_file(csv_save)
        
    def test_cl_rarefraction_curves(self):
        csv_save = os.path.join(self.scripts_dir, "temp_dir", "cl_rarefraction_curves.csv")
        result = subprocess.run(["python", os.path.join(self.scripts_dir, "cl_rarefraction_curves.py"), "--save_csv", csv_save, "--sequencing_report", self.sequencing_report, 
                                 "--samples", self.sample_name, self.second_sample_name, "--region_plots", self.region_plots], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0)
        self.remove_file(csv_save)
        result = subprocess.run(["python", os.path.join(self.scripts_dir, "cl_rarefraction_curves.py"), "--save_csv", csv_save, "--sequencing_report", self.sequencing_report, 
                                 "--region_plots", self.region_plots], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0)
        self.remove_file(csv_save)
        result = subprocess.run(["python", os.path.join(self.scripts_dir, "cl_rarefraction_curves.py"), "--save_csv", csv_save, "--sequencing_report", self.seq_report_fr3, 
                                 "--region_plots", self.alternative_region], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, "Alternative region failed")
        self.remove_file(csv_save)
        
    def test_cl_protein_embedding_umap(self):
        csv_save = os.path.join(self.scripts_dir, "temp_dir", "cl_protein_embedding_umap.csv")
        embedding_vector = os.path.join(self.temp_dir, "embedding_vector.npz")
        result = subprocess.run(["python", os.path.join(self.scripts_dir, "cl_protein_embedding_umap.py"),
                                 "--save_csv", csv_save, "--sequencing_report", self.sequencing_report, 
                                 "--samples", self.sample_name, self.second_sample_name, "--region_plots", self.region_plots,
                                 "--embedding_vector_path", embedding_vector, "--batch_size", "100"], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr + result.stdout)
        self.remove_file(csv_save)
        result = subprocess.run(["python", os.path.join(self.scripts_dir, "cl_protein_embedding_umap.py"),
                                 "--save_csv", csv_save, "--sequencing_report", self.sequencing_report, 
                                 "--samples", self.sample_name, self.second_sample_name, "--region_plots", self.region_plots,
                                 "--embedding_vector_path", embedding_vector, "--binding_data", self.binding_data, "--antigen_names", "Antigen 1", "--batch_size", "100"], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0)
        self.remove_file(csv_save)
        result = subprocess.run(["python", os.path.join(self.scripts_dir, "cl_protein_embedding_umap.py"),
                                 "--save_csv", csv_save, "--sequencing_report", self.seq_report_fr3, 
                                 "--samples", self.sample_name, self.second_sample_name, "--region_plots", self.alternative_region,
                                 "--embedding_vector_path", embedding_vector, "--binding_data", self.binding_data, "--antigen_names", "Antigen 1", "--batch_size", "100"], capture_output=True, text=True)
        self.assertEqual(result.returncode, 1) # because binding data does not have the col for the corresponding sequence
        self.remove_file(csv_save)
        result = subprocess.run(["python", os.path.join(self.scripts_dir, "cl_protein_embedding_umap.py"),
                                 "--save_csv", csv_save, "--sequencing_report", self.seq_report_fr3, 
                                 "--samples", self.sample_name, self.second_sample_name, "--region_plots", self.alternative_region,
                                 "--embedding_vector_path", embedding_vector], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, "Alternative region failed")
        self.remove_file(csv_save)
        
    def test_cl_protein_embedding_tsne(self):
        csv_save = os.path.join(self.scripts_dir, "temp_dir", "cl_protein_embedding_tsne.csv")
        embedding_vector = os.path.join(self.temp_dir, "embedding_vector.npz")

        result = subprocess.run(["python", os.path.join(self.scripts_dir, "cl_protein_embedding_tsne.py"), "--save_csv", csv_save, "--sequencing_report", self.sequencing_report,
                                 "--samples", self.sample_name, self.second_sample_name, "--region_plots", self.region_plots, "--batch_size", "300",
                                 "--embedding_vector_path", embedding_vector], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr + result.stdout)
        self.remove_file(csv_save)
        result = subprocess.run(["python", os.path.join(self.scripts_dir, "cl_protein_embedding_tsne.py"), "--save_csv", csv_save, "--sequencing_report", self.sequencing_report,
                                 "--samples", self.sample_name, self.second_sample_name, "--region_plots", self.region_plots,"--batch_size", "300",
                                 "--embedding_vector_path", embedding_vector, "--model_type", "nilsho01/LittleNano"], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr + result.stdout)
        self.remove_file(csv_save)
        
        

        
    
    