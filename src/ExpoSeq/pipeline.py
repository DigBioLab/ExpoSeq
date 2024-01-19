from .plots import barplot, cluster_embedding, embedding_with_binding, hist_lvst_dist, length_distribution, \
    logo_plot, plt_heatmap, protein_embedding, protein_network_embedding, rarefraction_curves, relative_sequence_abundance, stacked_aa_distribution, levenshtein_clustering, sample_cluster, \
        clone_fraction, diversity_plot
import matplotlib.pyplot as plt
from .augment_data.binding_data import collect_binding_data
from .augment_data.uploader import upload
import pandas as pd
import pickle
import os
from .settings import change_settings, reports, plot_styler
import subprocess
from .augment_data.uploader import create_alignment_report
from .settings.general_instructions import print_instructions
from .tidy_data.heatmaps.read_matrix import read_matrix
from .augment_data.mixcr_nils import check_mixcr
from .settings.aumotative_report import AumotativeReport
from .settings.figure import MyFigure, save_matrix
from .settings.markdown_builder import create_quarto


class PlotManager:
    def __init__(self,experiment = None, test_version=False,  length_threshold=6,
                 min_read_count=3, no_automation = False, module_dir = None, allow_binding_data = True, remove_gaps = True):
        """
        :param experiment: the name of the experiment you want to analyse
        :param test_version: Some settings are different if the pipeline is launched in test mode
        :param length_threshold: Filter all amino acid sequences smaller than this value
        :param min_read_count: Filter all amino acid sequences with less reads than this value
        :param no_automation: If True, the pipeline will not launch the automatic plotting tool
        :param module_dir: If you want to use the pipeline in a different directory than the main directory, you can specify it here
        :param allow_binding_data: Default is True. If set to False, the pipeline will not ask for binding data
        :param remove_gaps: Default is True. If set to False, the pipeline will not remove sequences with gaps
        """
        self.is_test = test_version
        self.Settings = change_settings.Settings()
        self.Settings.check_dirs()
        self.global_params = self.Settings.read_global_vars()
        self.remove_gaps = remove_gaps
        if module_dir:
            self.module_dir = module_dir
        else:    
            self.module_dir = self.Settings.module_dir

        if experiment == None:
            self.sequencing_report, self.alignment_report, self.experiment = upload()
        else:
            self.sequencing_report  = pd.read_csv(os.path.join(self.module_dir,
                                                                "my_experiments",
                                                                experiment,
                                                                "sequencing_report.csv"))
            self.alignment_report = create_alignment_report(self.module_dir, experiment)
            self.experiment = experiment
        self.region_string = self.global_params["region_of_interest"]
        if self.region_string == '':
            self.region_string = "CDR3"
        binding_report = reports.BindingReport(self.module_dir,
                                                self.experiment)
        if allow_binding_data != False:
            if os.path.isfile(allow_binding_data):
                self.binding_data = pd.read_csv(allow_binding_data)
            else:
                self.binding_data = binding_report.ask_binding_data()
        else:
            self.binding_data = None
    # check general settings and load them
        self.region_of_interest = "aaSeq" + self.region_string
        self.plot_path, self.mixcr_plots_path, self.experiment_path, self.report_path = self.Settings.check_dirs_automation(self.experiment,
                                                                                                                            self.region_of_interest)
        self.font_settings = self.Settings.read_font_settings()
        self.legend_settings = self.Settings.read_legend_settings()
        self.colorbar_settings = self.Settings.read_colorbar_settings()
        #experiment_path = self.Settings.get_experiment_path(self.experiment)
        
        self.Report = reports.SequencingReport(self.sequencing_report)
        self.Report.check_sample_name(module_dir=self.module_dir, experiment_name = self.experiment)
        self.experiments_list  = self.sequencing_report["Experiment"].unique().tolist()
        self.Report.prepare_seq_report(self.region_string,
                                       length_threshold=length_threshold,
                                       min_read_count=min_read_count,
                                       remove_gaps = self.remove_gaps)
       # self.Report.map_exp_names(self.unique_experiments)
        self.avail_regions = self.Report.get_fragment()
        self.sequencing_report = self.Report.sequencing_report

        self.ControlFigure = MyFigure(test_version)
        self.ControlFigure.set_backend() 
        self.style = plot_styler.PlotStyle(self.ControlFigure.ax, self.ControlFigure.plot_type)
        # self.settings_saver = change_save_settings.Change_save_settings()
        self.preferred_sample = self.experiments_list[0]
        self.java_heap_size = 1000
        self.change_preferred_antigen()
        self.Automation = None
        if self.Settings.automation == True and no_automation == False:
            self.full_analysis()
            
        print_instructions()

    def create_report(self):
        self.Settings.move_markdown_files()
        print("Install quarto from: https://quarto.org/docs/get-started/")
        print("You also might to install Jupyter. In general, I recommend to use VS code for this purpose.")
        print(f"This function should only be called after the initial automation, since it relies on the file structure in {self.plot_path}\nAlso do not change the names of the directories or files")
        print("Test Quarto")
        subprocess.run(["quarto", "check"])
        print("Render qmd file")
        create_quarto(self.experiment, self.plot_path, self.binding_data, self.experiments_list)
        subprocess.run(["quarto", "render", os.path.join(self.plot_path, self.experiment + ".qmd")])
        print(f"You can find the html file to the report at: {self.plot_path}/{self.experiment}.qmd")



    def change_java_heap_size(self, new_size):
        self.java_heap_size = new_size
    
    def export_sequencing_report(self):
        self.sequencing_report.to_csv(os.path.join(self.experiment_path, "sequencing_report_processed.csv"))
        print(f"Lowest peptide sequence length: {self.sequencing_report['aaSeqCDR3'].str.len().min()}")
        print(f"lowest read count: {self.sequencing_report['readCount'].min()}")
        
    

    def chat(self,):
        """
        :return: starts a conversation with your data. Be creative! You can ask any question and even create plots :)
        """
        if self.Automation == None:
            if self.binding_data is not None:
                merged_reprot = self.merge_bind_seq_report()
                print("Your data looks like: \n")
                merged_reprot.head(10)
                self.Automation = AumotativeReport(merged_reprot,
                                                   self.Report.origin_seq_report,
                                                   self.global_params)
            else:
                print("Your data looks like: \n")
                self.sequencing_report.head(10)
                self.Automation = AumotativeReport(self.sequencing_report,
                                                self.Report.origin_seq_report,
                                                self.global_params)
            
        response = self.Automation.chat_()
        return response

    def save_in_plots(self, enter_filename):
        plt.savefig(fname = os.path.join(self.plot_path, enter_filename + f".png"), dpi = 300,  format="png",  bbox_inches='tight')

    def get_best_binder(self):
        if self.binding_data is None:
            return None
        found_antigens = self.binding_data.columns.to_list()[1:]
        best_binding = {}

        for i in self.experiments_list:
            sub_report = self.sequencing_report[self.sequencing_report["Experiment"] == i]
            merged_report = sub_report.merge(self.binding_data, how = "inner")
            clone_fractions = merged_report["cloneFraction"]
            clone_fractions.fillna(0, inplace = True)
            best_binding_local = {}
            for antigen in found_antigens:
                binding_values = merged_report[antigen]
                binding_values.fillna(0, inplace = True)
                significance_binding = (binding_values * clone_fractions).mean()
                if antigen in best_binding_local.keys():
                    if best_binding_local.get(antigen) < significance_binding:
                        best_binding_local[antigen] = significance_binding
                    else:
                        pass
                else:
                    best_binding_local[antigen] = significance_binding

            sorted_keys = sorted(best_binding_local, key=best_binding_local.get, reverse=True)

            # Get the best and second best keys
            best_key = sorted_keys[0]
            if len(found_antigens) >= 3:
                second_best_key = sorted_keys[1] if len(sorted_keys) > 1 else None
                best_binding[i] = [best_key, second_best_key]
            else:
                best_binding[i] = [best_key]
        return best_binding
    
    
    def model_protein():
        print("Visit one of these pages and explore modelling proteins for free!"
              f"ESM FOLD: https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/ESMFold.ipynb#scrollTo=CcyNpAvhTX6q")

    def full_analysis(self):
        """
        :return: returns some plots about quality in /my_experiments/YOUR_EXPERIMENT_NAME/plots. Is called automatically when launching the pipeline.
        """
        self.plot_path, self.mixcr_plots_path, self.experiment_path, self.report_path = self.Settings.check_dirs_automation(self.experiment,
                                                                                                                                self.region_of_interest)
        self.plot_path = os.path.join(self.module_dir,
                                      "my_experiments",
                                      self.experiment,
                                      "plots",
                                      self.region_of_interest)
        report_seq_cluster = os.path.join(self.plot_path, "sequence_cluster", "reports")
        jaccard_report = os.path.join(self.report_path, "jaccard_identity" + ".xlsx")
        if not os.path.isdir(self.plot_path):
            os.mkdir(self.plot_path)
        if not os.path.isdir(os.path.join(self.plot_path, "length_distributions")):
            os.mkdir(os.path.join(self.plot_path, "length_distributions"))
            
        if not os.path.isdir(os.path.join(self.plot_path, "clone_fraction")):
            os.mkdir(os.path.join(self.plot_path, "clone_fraction"))
        if not os.path.isdir(os.path.join(self.plot_path, "rarefraction_curves")):
            os.mkdir(os.path.join(self.plot_path, "rarefraction_curves"))
        if not os.path.isdir(os.path.join(self.plot_path, "logo_plots")):
            os.mkdir(os.path.join(self.plot_path, "logo_plots"))
        if not os.path.isdir(os.path.join(self.plot_path, "sequence_cluster")):
            os.mkdir(os.path.join(self.plot_path, "sequence_cluster"))

        if not os.path.isdir(os.path.join(self.plot_path, "sequence_embedding")):
            os.mkdir(os.path.join(self.plot_path, "sequence_embedding"))
            
        if not os.path.isdir(os.path.join(self.plot_path, "sequence_embedding", "sgt")):
            os.mkdir(os.path.join(self.plot_path, "sequence_embedding", "sgt")) 
            
        if not os.path.isdir(os.path.join(self.plot_path, "sequence_embedding", "protbert")):
            os.mkdir(os.path.join(self.plot_path, "sequence_embedding", "protbert")) 
            
        if not os.path.isdir(os.path.join(self.plot_path, "sequence_embedding", "T5")):
            os.mkdir(os.path.join(self.plot_path, "sequence_embedding", "T5")) 
        
        if not os.path.isdir(report_seq_cluster):
            os.mkdir(report_seq_cluster)
            
        #plot 
        try:

            mh_save_path = os.path.join(self.report_path, "morosita_horn_identity" + ".xlsx")
            if not os.path.isfile(os.path.join(self.plot_path, "morosita_horn.png")):
                print("Start preparing morosita horn matrix.")
                self.morosita_horn(matrix_save_path = mh_save_path)
                self.save_in_plots("morosita_horn")
        except:
            print("Creating Morosita Horn matrix failed")
        #plot
        try:
            if not os.path.isfile(os.path.join(self.plot_path, "sorensen.png")):
                print("Start preparing sorensen matrix.")

                sorensen_save_path = os.path.join(self.report_path, "sorensen_identity" + ".xlsx")
                self.sorensen(matrix_save_path = sorensen_save_path)
                self.save_in_plots("sorensen")
        except:
            print("Creating sorense matrix failed")
        #plot
        try:
            if not os.path.isfile(os.path.join(self.plot_path, "jaccard.png")):
                print("Start preparing jaccard matrix.")

                jaccard_save_path = os.path.join(self.report_path, "jaccard_identity" + ".xlsx")
                self.jaccard(matrix_save_path = jaccard_save_path)
                self.save_in_plots("jaccard")
        except:
            print("Creating sorense matrix failed")

        # clone fract plot
        for experiment in self.experiments_list:
            try:
                self.rel_seq_abundance(experiment)
                self.save_in_plots(os.path.join("clone_fraction", experiment))
            except:
                print(f"Clone fraction for {experiment} failed")
        # diversity plot
        
        for method in ["Shannon", "InverseSimpson"]:
            try:
                self.sample_diversity(method)
                self.save_in_plots(os.path.join("diversity_" + method))
            except:
                print(f"Diversity for {method} plot failed")        
        #plot
        for experiment in self.experiments_list:
            try:
                self.lengthDistribution_single(experiment)
                self.save_in_plots(os.path.join("length_distributions", experiment))
            except:
                print(f"Length Distribution for {experiment} failed")
        #mixcr plots 
        try:
            print("Begin with creating mixcr plots. You can find a description about them here:\n https://mixcr.com/mixcr/reference/mixcr-exportPlots/")
            self.export_mixcr_quality()
            self.mixcr_explain_diversity()
            self.mixcr_vdj_usage()
          #  self.mixcr_overlap()
        except:
            print(f".clns files could not be found. Most certainly, you have uploaded the sequencing report.")
        # plot
        try:
            self.alignment_quality()
            self.save_in_plots("alignment_quality")
        except:
            print("Alignment reports could not be found")
            

        #plot
        try:
            if not os.path.isfile(os.path.join(self.plot_path, "rarefraction_all.png")):
                print("create rarefraction curves of all samples")
                self.rarefraction_curves(self.experiments_list)
                self.ControlFigure.fig.tight_layout()
                self.save_in_plots("rarefraction_all")
        except:
            print("creating rarefraction plot of all samples failed")
            
        if not os.path.isdir(os.path.join(self.plot_path, "rarefraction_curves")):
            os.mkdir(os.path.join(self.plot_path, "rarefraction_curves"))

        for i in self.experiments_list:
            if not os.path.isfile(os.path.join("rarefraction_curves", "rarefraction_curves_"+i + ".png")):
                try:
                    self.rarefraction_curves(samples = [i])
                    self.save_in_plots(os.path.join("rarefraction_curves", "rarefraction_curves_"+i))
                except:
                    pass
            
        # plot
        for experiment in self.experiments_list:
            try:
                if not os.path.isfile(os.path.join(self.plot_path, "logo_plots", experiment + ".png")):
                    self.logoPlot_single(sample = experiment)
                    self.save_in_plots(os.path.join("logo_plots", experiment))
            except:
                print(f"Logo Plot for {experiment} failed")


        print("Cluster NGS sequences in dendrograms and network plots using levenshtein distance of 2 for network plots and ls-distance of 1 for dendrogram. Batch size is set to 1000.")        
        #plot
        try:
            print(f"Create clustering for all samples for top 50% of clone fraction and for maximum of 100 reads")
            self.connect_samples()
            self.save_in_plots(os.path.join(self.plot_path, "ls_connection_all"))
        except:
            print("Creation of cluster for all samples failed.")
        for single_experiment in self.experiments_list:
            try:
                if not os.path.isfile(os.path.join(self.plot_path, "sequence_cluster", single_experiment + "ls_dendro.png")):
                    self.levenshtein_dendrogram(sample = single_experiment, max_cluster_dist = 1)
                    self.save_in_plots(os.path.join("sequence_cluster", single_experiment + "ls_dendro"))
                if not os.path.isfile(os.path.join(self.plot_path,
                                                    "sequence_cluster",
                                                    single_experiment + "ls_cluster.png")):
                    self.basic_cluster(samples = [single_experiment],
                                        max_ld = 2,
                                        batch_size=1000,
                                        label_type = None,
                                        save_report_path = os.path.join(report_seq_cluster, f"{single_experiment}" + "ls_cluster" + ".xlsx"))
                    self.save_in_plots(os.path.join("sequence_cluster", single_experiment + "ls_cluster"))
            except:
                print(f"Sequence clustering failed at: {single_experiment}")
        threshold_identity = 0.2
        
        while True:
            overlapping_samples = read_matrix(threshold=threshold_identity, path_matrix=mh_save_path)
            # if None path_matrix does not exist
            if overlapping_samples == None:
                if mh_save_path == jaccard_report:
                    overlapping_samples = None
                    break
                else:  
                    mh_save_path = jaccard_report
                    continue
            else:
                mean_length = sum(len(v) for v in overlapping_samples.values()) / len(overlapping_samples)
                if mean_length > 2.5:
                    break
                elif threshold_identity < 0.05:
                    break
                else:
                    threshold_identity -= 0.01
                    

        print("Start sequence embedding of samples with 1000 iterations for tsne and a batch size of 300")
        if overlapping_samples != None:
            for single_experiment in list(overlapping_samples.keys()):
                if not os.path.isfile(os.path.join(self.plot_path, "sequence_embedding","sgt", single_experiment + "embedding_tsne.png")):
                    try:
                        self.embedding_tsne(samples = overlapping_samples[single_experiment],model = "sgt", strands = False, batch_size = 300, iterations_tsne = 1000)
                        self.save_in_plots(os.path.join("sequence_embedding","sgt", single_experiment + "embedding_tsne"))
                    except:
                        print(f"Clustering with TSNE failed for {single_experiment}.")
        
                        
        if overlapping_samples != None:
            for single_experiment in list(overlapping_samples.keys()):
                if not os.path.isfile(os.path.join(self.plot_path, "sequence_embedding","T5", single_experiment + "embedding.png")):
                    try:
                        self.embedding_tsne(samples = overlapping_samples[single_experiment], strands = False, batch_size = 300)
                        self.save_in_plots(os.path.join("sequence_embedding","T5", single_experiment + "embedding_tsne"))
                    except:
                        print(f"Clustering with TSNE failed for {single_experiment}.")   
                
        else: 
            print("Neither Morosita Horn nor Jaccard matrix was generated, thus no sequence embedding can be created.")
        if self.region_of_interest not in self.binding_data.columns.tolist():
            print(f"There is no binding data available for {self.region_of_interest}")
            best_binder = None
        else:
            best_binder = self.get_best_binder()    
            
        if best_binder == None:
            print("No binding data was found. Thus, no binding plots can be created.")

        else:
            clustering_antigens_path = os.path.join(self.plot_path, "clustering_antigens")
            if not os.path.isdir(clustering_antigens_path):
                os.mkdir(clustering_antigens_path)
            experiment_keys = list(best_binder.keys())
            print("Starting clustering of binding data for antigens and sequences using PCA and t-SNE.\nPerplexity is set to 25.\n70 principal components are taken for t-SNE.\n1000 iterations are used for tsne. Batch size is set to 300.")
            report_tsne_cluster = os.path.join(self.plot_path, "clustering_antigens", "reports")
            if not os.path.isdir(report_tsne_cluster):
                os.mkdir(report_tsne_cluster)
                
            #plot
            for experiment in experiment_keys:
                try:
                    if not os.path.isfile(os.path.join(self.plot_path, "clustering_antigens", experiment + "tsne_cluster_AG.png")):
                        self.cluster_binding_data(samples = [experiment],
                                                  batch_size = 300,
                                                    antigens = best_binder[experiment],
                                                    show_antigen_names = False,
                                                    iterations_tsne=1000,
                                                    save_report_path = os.path.join(report_tsne_cluster, experiment + f"_best_binder{experiment}" + ".xlsx"))
                        self.save_in_plots(os.path.join( "clustering_antigens", experiment + "tsne_cluster_AG"))
                except:
                    print(f"Cluster based on sequence embedding combined with binding data failed for: {experiment}")
            dendro_path = os.path.join(self.plot_path,"clustering_antigens", "dendro_binding")
            
            if not os.path.isdir(dendro_path):
                os.mkdir(dendro_path)
            ls_cluster_path = os.path.join(self.plot_path, "clustering_antigens", "ls_binding_cluster")
            if not os.path.isdir(ls_cluster_path):
                os.mkdir(ls_cluster_path)
            print("Create dendrogram with binding data and levenshtein distance of 2 and batch size of 300 for the dendrogram and 1000 for the network plots.")
            
            report_ls_cluster = os.path.join(self.plot_path, "clustering_antigens","ls_binding_cluster", "reports")
            if not os.path.isdir(report_ls_cluster):
                os.mkdir(report_ls_cluster)
            #plot
            for experiment in experiment_keys:
                try:
                    plt.close('all')
                    fig = self.dendro_bind(sample = experiment,
                                     antigens = best_binder[experiment],
                                     batch_size=300
                                     )
                    if fig == False:
                        continue
                    else:
                        plt.close(fig)
                        self.save_in_plots(os.path.join("clustering_antigens","dendro_binding", experiment + "cluster_dendrogram"))
                        try:
                            plt.close('all')
                            self.cluster_one_AG(antigen = best_binder[experiment][0],
                                                specific_experiments = [experiment],
                                                batch_size = 1000,
                                                max_ld = 2,
                                                preferred_cmap = "Reds",
                                                label_type = None,
                                                save_report_path = os.path.join(report_ls_cluster, experiment + f"_ls_binding_cluster" + ".xlsx"))
                            self.save_in_plots(os.path.join("clustering_antigens","ls_binding_cluster", experiment + f"_ls_binding_cluster"))
                        except:
                            print("Could not create levenshtein distance cluster for best binder")
                except:
                    print(f"Dendrogram with binding data failed for {experiment}")


        print(f"The pipeline has created some plots in: {self.plot_path}")
                    
        print("You can create and automatic report of the generated analysis by installing Quarto(https://quarto.org/docs/get-started/) and type: plot.create_report()")

    def change_preferred_antigen(self, antigen=None):
        if self.binding_data is not None:
            columns_bind_data = self.binding_data.columns.to_list()
            if antigen == None:
                antigen = columns_bind_data[1]
            self.preferred_antigen = antigen
        else:
            self.preferred_antigen = None

    def change_preferred_sample(self, sample):
        assert sample in self.experiments_list, "The provided sample name is not in your sequencing report. Please check the spelling or use the print_samples function to see the names of your samples"
        self.preferred_sample = sample

    def change_filter(self, divisible_by=3, length_threshold_aa=6, min_read_count=2):
        """
        :param divisible_by: filter all sequences which are not divisible by this value
        :param length_threshold: filter out all sequences which are shorter than this value
        :min_read_count: filter out all sequences which have less clones/reads than this value
        """
        self.Report.prepare_seq_report(self.region_string, divisible_by, length_threshold_aa, min_read_count)
        self.sequencing_report = self.Report.sequencing_report

    def change_region(self, region = None):
        """
        :return: changes the region you want to analyse
        """
        intermediate = self.region_of_interest
        possible_regions = self.avail_regions 
        if region == None:
            region_string = input(
                f"Which region do you want to plot? ExpoSeq could find the following regions: {self.avail_regions}.")
        else:
            region_string = region
        if region_string in possible_regions:
            if not region_string == "targetSequences":
                region_string = region_string.replace("nSeq", "")
                self.Report.prepare_seq_report(region_string, divisible_by=3, length_threshold=9, min_read_count=0, remove_gaps = self.remove_gaps)
            else:
                region_string = "targetSequences"
                self.Report.prepare_seq_report(region_string, divisible_by=3, length_threshold=9, min_read_count=0, remove_gaps = self.remove_gaps )
            self.sequencing_report = self.Report.sequencing_report
            self.region_of_interest = "aaSeq" + region_string
        else:
            print(f"The region you want to plot is not valid. The options are: {self.avail_regions}")
            self.region_of_interest = intermediate
        plot_dir = os.path.dirname(self.plot_path)
        self.plot_path = os.path.join(plot_dir, self.region_of_interest)

    def discard_samples(self, samples_to_discard):
        """

        :param samples_to_discard: the name of the samples in a list which you want to delete from your further analysis
        :return: updates the sequencing report where the values for the given samples do not exist anymore

        """
        assert type(samples_to_discard) == list, "You have to give a list with the samples you want to discard"
        self.sequencing_report = self.sequencing_report[~self.sequencing_report['Experiment'].isin(samples_to_discard)]

    def add_binding_data(self):
        """"
        :return: adds binding data to your analysis. You can add mutliple files with the filechooser or a given path to the file. For more information about the necessary file structure, check the supplementary information.
        """
        self.binding_data = collect_binding_data(binding_data=self.binding_data)

        self.change_preferred_antigen(self.binding_data.columns.to_list()[1])

        # self.binding_data.to_csv("binding_data.csv")

    def merge_bind_seq_report(self, merge_technique = "left"):
        """
        :param: merge_technique: default "left". You can choose between "left", "right", "inner" and "outer" for the merge technique.
        :return: a python variable which contains the merged sequencing report and binding report which can be saved with report.to_csv("filename.csv")
        """
        assert self.binding_data is not None, "You have not given binding data. You can add it with the add_binding_data function"

        merged_reports = pd.merge(self.sequencing_report, self.binding_data, on="aaSeqCDR3", how=merge_technique)
        merged_reports.fillna(0, inplace=True)
        return merged_reports

    def print_antigens(self):
        """
        :return: prints the antigens (columns) of your binding data
        """
        try:
            print(self.binding_data.columns.to_list()[1:])
        except:
            print("You have not given binding data.")

    def print_samples(self):
        """
        :return: prints the names of your samples, so you can insert them in lists or similar for the analysis with some plots
        """
        print(self.experiments_list)

    def save(self, enter_filename, dpi=300, format="png"):
        """
        :return: can be used to save your plots. It will ask you automatically for the directory where you want to save it
        """
        print("Choose the directory where you want to save your plot with the filechooser")
        try:
            from tkinter import filedialog
            save_dir = filedialog.askdirectory()
        except:
            save_dir = ""

        plt.savefig(fname=os.path.join(save_dir, enter_filename + f".{format}"), dpi=dpi, format=format)

    def close(self):
        plt.close()

    def change_experiment_names(self, specific=None):
        """
        :param specific: optional parameter. You can use this function to change the names of a specific sample.
        :return:You can use this function to change the name of your samples. Thus, you can change the labels of your plots.
        """

        self.expreiments_list = self.sequencing_report["Experiment"].unique().tolist()
        experiments_dic = dict(zip(self.experiments_list, self.experiments_list))
        if specific != None:
            assert type(
                specific) == str, "You have to give a string (inside: "") as input for the specific sample you want to change"

        if specific == None:
            print(
                "You will no go through all your samples and you will be able to change their names. If you want to skip a samples, just press enter.")
            for key in self.experiments_list:
                print(f"Current value for {key}: {experiments_dic[key]}")
                new_value = input("Enter a new value or press any key to skip")
                if len(new_value) > 1:
                    experiments_dic[key] = new_value
                else:
                    pass
        else:
            new_value = input("Enter the new name for " + specific)
            experiments_dic[specific] = new_value
        self.Report.renew_exp_names_origin(experiments_dic, self.experiment_path)
        self.experiments_list = list(experiments_dic.values())
        self.sequencing_report["Experiment"] = self.sequencing_report["Experiment"].map(experiments_dic)

        
        
    def alignment_quality(self, log_transformation=False):
        """

        :param log_transformation: optional parameter. You can set it to True if you want to log transform your data
        :return: Creates a figure where you can see the Overall Reads per sample and the number of reads which could be aligned.
        """
        self.ControlFigure.check_fig()
        self.ControlFigure.plot_type = "single"
        assert log_transformation in [True, False], "You have to give True or False as input for the log transformation"
        self.ControlFigure.clear_fig()
        try:
            barplot.barplot(self.ControlFigure.ax,
                            self.alignment_report,
                            self.sequencing_report,
                            self.font_settings,
                            self.legend_settings,
                            apply_log=log_transformation)

            self.ControlFigure.update_plot()
            self.style = plot_styler.PlotStyle(self.ControlFigure.ax,
                                               self.ControlFigure.plot_type)
        except:
            print("Currenlty you do not have loaded an alignment report.")

    def aa_distribution(self, sample=None, region=[3, 7], protein=True):
        """
        :param sample: The sample you would like to analyze
        :param region: The input is a list and you can specify with the first value the beginning region and with the second value the end region within your sequences for which you want to know the amino acid composition. For instance: [3,7]
        :param protein: Default True. If you would like to analyze nucleotide sequences, set it to False
        :return: Returns a plot which shows the amino acid distribution in the given sequence range.
        """
        if sample == None:
            sample = self.preferred_sample
        self.ControlFigure.check_fig()
        self.ControlFigure.plot_type = "single"
        assert sample in self.experiments_list, "The provided sample name is not in your sequencing report. Please check the spelling or use the print_samples function to see the names of your samples"
        assert type(sample), "You have to give a string as input for the sample"
        assert type(
            region), "You have to give a list with the start and end position of the region you want to analyze. For instance: [3,7]"
        assert protein in [True, False], "You have to give True or False as input for the protein parameter"
        self.ControlFigure.clear_fig()
        stacked_aa_distribution.stacked_aa_distr(self.ControlFigure.ax,
                                                 self.sequencing_report,
                                                 sample,
                                                 region,
                                                 protein,
                                                 self.font_settings,
                                                 self.region_of_interest)

        self.ControlFigure.update_plot()
        self.style = plot_styler.PlotStyle(self.ControlFigure.ax,
                                           self.ControlFigure.plot_type)

    def rarefraction_curves(self, samples=None):
        """
        :param samples: you insert a list which contains the sample names
        :return: Shows your the rarefraction curves for the given samples.
        """
        if samples == None:
            samples = [self.experiments_list[0]]
        incorrect_samples = [x for x in samples if x not in self.experiments_list]
        assert not incorrect_samples, f"The following sample(s) are not in your sequencing report: {', '.join(incorrect_samples)}. Please check the spelling or use the print_samples function to see the names of your samples"

        self.ControlFigure.check_fig()
        self.ControlFigure.plot_type = "single"
        if type(samples) != list:
            print(
                "You have to give the sample names in the list, also if it is only one! A list is a container which is marked through:[] . Please try again.")
        else:
            self.ControlFigure.clear_fig()
            rarefraction_curves.RarefractionCurves(
                              sequencing_report=self.sequencing_report,
                              samples=samples,
                              region_of_interest=self.region_of_interest,
                              ax=self.ControlFigure.ax,
                              font_settings=self.font_settings,
                              legend_settings=self.legend_settings,
                              )

            self.ControlFigure.update_plot()
            self.style = plot_styler.PlotStyle(self.ControlFigure.ax,
                                               self.ControlFigure.plot_type)

    def logoPlot_single(self,
                        sample=None,
                        highlight_specific_pos=None,
                        method="proportion",
                        chosen_seq_length=None,
                        color_scheme = "chemistry",
                        **kwargs):
        """
        :param sample: insert the sample name
        :param highlight_specific_pos: optional. you can highlight a specific position. For instance if you want to highlight the 3rd position, you insert 3.
        :param method: You can specify whether you want to have on your y axis the frequency of the amino acids or the information content in bits. The default is proportion. If you want to have the information content, insert "bits"
        :param chosen_seq_length: Max per default. If you plot the length distribution it will take the sequences with the highest bar there. You always analyze only one sequence length! You can enter an integer to analyze other lengths.
        :param color_scheme: Default is skylign_protein. You can choose between skylign_protein, chemistry, weblogo_protein, hydrophobicity, dmslogo_charge, charge, NajafabadiEtAl2017.
        :return: A logo Plot which shows you the composition of aminoacids per position
        """
        logo_schemes = ["skylign_protein", "chemistry", "weblogo_protein", "hydrophobicity", "dmslogo_charge", "charge", "NajafabadiEtAl2017"]
        assert color_scheme in logo_schemes, f"You have to give a valid color scheme. The options are: {logo_schemes}"
        if sample == None:
            sample = self.preferred_sample
        self.ControlFigure.check_fig()
        self.ControlFigure.plot_type = "single"
        assert sample in self.experiments_list, "The provided sample name is not in your sequencing report. Please check the spelling or use the print_samples function to see the names of your samples"
        assert type(sample) == str, "You have to give a string as input for the sample"
        if highlight_specific_pos != None:
            assert type(
                highlight_specific_pos) == int, "You have to give an integer as input for the specific position you want to highlight"
        if chosen_seq_length != None:
            assert type(chosen_seq_length) == int, "You have to give an integer as input for the sequence length you want to analyze"
        self.ControlFigure.clear_fig()

        logo_plot.LogoPlot(self.ControlFigure.ax,
                           self.sequencing_report,
                            self.region_of_interest,
                           sample,
                            highlight_specific_pos,
                            self.font_settings,
                           chosen_seq_length,
                           method,
                           color_scheme,
                           **kwargs
                           )

        self.ControlFigure.update_plot()
        self.style = plot_styler.PlotStyle(self.ControlFigure.ax,
                                           self.ControlFigure.plot_type)

    def logoPlot_multi(self,
                       samples="all",
                       chosen_seq_length=1,
                       method="proportion",
                       ):
        """
        :param samples: You analyze all samples per default. If you want to analyze specific samples it has to be a list with the corresponding sample names
        :param chosen_seq_length: 16 per default. You always analyze online one sequence length! You can change it if you would like
        :return: Gives you in one figure one logoPlot per sample.
        """
        self.ControlFigure.check_fig()
        self.ControlFigure.plot_type = "multi"
        if samples != "all":
            assert type(samples) == list, "You have to give a list with the samples you want to analyze"
            incorrect_samples = [x for x in samples if x not in self.experiments_list]
            assert not incorrect_samples, f"The following sample(s) are not in your sequencing report: {', '.join(incorrect_samples)}. Please check the spelling or use the print_samples function to see the names of your samples"
        if chosen_seq_length != "max":
            assert type(chosen_seq_length) == int, "You have to give an integer as input for the sequence length you want to analyze"
        self.ControlFigure.clear_fig()

        logo_plot.plot_logo_multi(self.ControlFigure.fig,
                                  self.sequencing_report,
                                  samples,
                                  self.font_settings,
                                  self.region_of_interest,
                                  method,
                                  chosen_seq_length,
                                  )
        self.ControlFigure.update_plot()
        self.style = plot_styler.PlotStyle(self.ControlFigure.ax,
                                           self.ControlFigure.plot_type)
      #  self.ControlFigure.tighten()

    def lengthDistribution_single(self, sample=None):
        """

        :param sample: name of the sample from which you would like to see the length distribution
        :return: Shows you the length Distribution of your sequences of the given sample
        """
        if sample == None:
            sample = self.preferred_sample
        self.ControlFigure.check_fig()
        self.ControlFigure.plot_type = "single"
        assert type(sample) == str, "You have to give a string as input for the sample"
        assert sample in self.experiments_list, "The provided sample name is not in your sequencing report. Please check the spelling or use the print_samples function to see the names of your samples"
        self.ControlFigure.clear_fig()
        length_distribution.length_distribution_single(self.ControlFigure.fig,
                                                       self.ControlFigure.ax,
                                                       self.sequencing_report,
                                                       sample,
                                                       self.font_settings,
                                                       self.region_of_interest)

        self.ControlFigure.update_plot()
        self.style = plot_styler.PlotStyle(self.ControlFigure.ax,
                                           self.ControlFigure.plot_type)


    def lengthDistribution_multi(self, samples="all"):
        """
        :param samples: You analyze all samples per default. If you want to analyze specific samples it has to be a list with the corresponding sample names
        :return: Outputs one figure with one subplot per sample which shows you the distribution of the sequence length
        """
        self.ControlFigure.check_fig()
        self.ControlFigure.plot_type = "multi"
        if samples != "all":
            assert type(samples) == list, "You have to give a list with the samples you want to analyze"
            incorrect_samples = [x for x in samples if x not in self.experiments_list]
            assert not incorrect_samples, f"The following sample(s) are not in your sequencing report: {', '.join(incorrect_samples)}. Please check the spelling or use the print_samples function to see the names of your samples"
        self.ControlFigure.clear_fig()
        length_distribution.length_distribution_multi(self.ControlFigure.fig,
                                                      self.sequencing_report,
                                                      samples,
                                                      self.font_settings,
                                                      self.region_of_interest,
                                                      test_version=self.is_test)
        self.ControlFigure.update_plot()
        self.style = plot_styler.PlotStyle(self.ControlFigure.ax,
                                           self.ControlFigure.plot_type)
      #  self.ControlFigure.tighten()

    def rel_seq_abundance(self, sample=None, visualize_sequences = True, prefered_cmap = "Reds", top_clone_fraction = 0.75, seqs_viz_fraction = 0.3, alpha_val = 0.75, pad_rectangles = True, force_reducing = 500, limit_seq_filter = 10):
        """
        :param samples: For a qualitative analysis choose samples from the same panning experiment. Input is a list
        :param visualize_sequences: Default is True. If you want to visualize a fraction of the sequences in the rectangles, set it to True
        :param prefered_cmap: Default is Reds. You can choose any colormap from matplotlib
        :param top_clone_fraction: Default is 0.95. You can choose the top clone fraction you want to visualize
        :param seqs_viz_fraction: Default is 0.75. You can choose the fraction of sequences you want to visualize
        :param alpha_val: Default is 0.5. You can choose the alpha value of the rectangles
        :param pad_rectangles: Default is True. If you want to have a padding between the rectangles, set it to True
        :param force_reducing: Default is 1000. If you have more than 1000 sequences after filtering by taking the top 95 % clones, this threshold automatically limits the sequences to 1000. You can set it to None if you want to avoid this filtering.
        """
        if sample == None:
            sample = self.preferred_sample
        self.ControlFigure.check_fig()
        self.ControlFigure.plot_type = "single"
        assert type(sample) == str, "You have to give a list with the samples you want to analyze"
        assert type(visualize_sequences) == bool, "You have to give True or False as input for the visualize_sequences parameter"
        assert type(prefered_cmap) == str, "You have to give a string as input for the prefered_cmap parameter"
        assert type(top_clone_fraction) == float, "You have to give a float as input for the top_clone_fraction parameter"
        assert top_clone_fraction < 1, "The top_clone_fraction has to be smaller than 1"
        assert top_clone_fraction > 0, "The top_clone_fraction has to be larger than 0"
        assert type(seqs_viz_fraction) == float, "You have to give a float as input for the seqs_viz_fraction parameter"
        assert seqs_viz_fraction < 1, "The seqs_viz_fraction has to be smaller than 1"
        assert seqs_viz_fraction > 0, "The seqs_viz_fraction has to be larger than 0"
        assert type(alpha_val) == float, "You have to give a float as input for the alpha_val parameter"
        assert alpha_val < 1, "The alpha_val has to be smaller than 1"
        assert alpha_val > 0, "The alpha_val has to be larger than 0"
        assert type(pad_rectangles) == bool, "You have to give True or False as input for the pad_rectangles parameter"
        self.ControlFigure.clear_fig()
        clone_fraction.VisFrac(self.sequencing_report, 
                               self.ControlFigure.ax,
                               sample,
                               visualize_sequences,
                               prefered_cmap,
                               top_clone_fraction,
                               seqs_viz_fraction,
                               alpha_val,
                               pad_rectangles,
                               force_reducing,
                               limit_seq_filter,
                               self.font_settings, 
                               self.region_of_interest)

        self.ControlFigure.update_plot()
        self.style = plot_styler.PlotStyle(self.ControlFigure.ax,
                                           self.ControlFigure.plot_type)

    def save_cluster_report(self,cluster_report, path=None):
        if path == None:
            while True:
                cluster_report_save = input("Do you want to save the generated data? (Y/n)")
                if cluster_report_save in ["Y", "y", "N", "n"]:
                    break
                else:
                    print("Please enter Y or n")
            if cluster_report_save.lower() in ["Y", "y"]:
                while True:
                    filename_matrix = input(f"Enter a name for the file. The file will be saved in {self.report_path} in your IDE.")
                    if filename_matrix.endswith(".csv"):
                        pass
                    else:
                        filename_matrix = filename_matrix + ".csv"
                    if not os.path.isfile(filename_matrix):
                        path = os.path.join(self.report_path, filename_matrix)
                        cluster_report.to_excel(path)
                        break
                    else:
                        print("This file already exists. Please choose another name.")
        else:
            cluster_report.to_excel(path)


    def basic_cluster(self, samples = None, batch_size = 1000, max_ld = 1, min_ld = 0, second_figure = False, label_type = "numbers", save_report_path = None):
        """
        :param samples: A list containing the samples you would like to analyze. Analyze just one sample with: ["My_sampel_name"].
        :param max_ld: Maximum allowed levenshtein distance between sequences within one cluster. The higher the distance the larger the clusters.
        :param min_ld: Minimum allowed levenshtein distance between sequences within one cluster.
        :param second_figure: Default is False. Creates a second figure with the levenshtein distance as bars
        :param label_type: Default is numbers. This will label the nodes in the plot with the corresponding identifier in the output report. You can type sequences for labeling the nodes with the sequneces. If you do not want to have labels set it to None.
        :param save_report_path: Default is None which saves your report in my_experiments/reports_pipeline. If you want to change it somewhere else you need to insert the full path with filename.
        """
        if samples == None:
            samples = [self.preferred_sample]
        self.ControlFigure.check_fig()
        self.ControlFigure.plot_type = "single"
      #  assert sample in self.experiments_list, "The provided sample name is not in your sequencing report. Please check the spelling or use the print_samples function to see the names of your samples"
        assert type(samples) == list, "You have to give a string as input for the sample"
        assert type(max_ld) == int, "You have to give an integer as input for the maximum levenshtein distance"
        assert type(min_ld) == int, "You have to give an integer as input for the minimum levenshtein distance"
        assert type(second_figure) == bool, "You have to give True or False as input for the second figure"
        self.ControlFigure.clear_fig()
        cluster_report = levenshtein_clustering.clusterSeq(
            self.ControlFigure.ax,
            self.sequencing_report,
            samples,
            max_ld,
            min_ld,
            batch_size,
            self.font_settings,
            second_figure,
            self.region_of_interest,
            label_type)

        self.ControlFigure.update_plot()
        self.style = plot_styler.PlotStyle(self.ControlFigure.ax,
                                           self.ControlFigure.plot_type)
        self.save_cluster_report(cluster_report = cluster_report,
                                 path = save_report_path)
        if second_figure == True:
            print("close the second window before you continue")

    def show(self):
        self.ControlFigure.fig.show()

    def cluster_binding_data(self, samples=None, strands=True, model = 'Rostlab/prot_t5_xl_half_uniref50-enc', pca_components=80, perplexity=30, iterations_tsne=2500, batch_size=1000, antigens = None, show_antigen_names = False, extra_figure = False, save_report_path = None):
        """
        :param samples: type is list. The samples you would like to compare towards their sequences
        :param strands: Default is True. It means that you will plot a batch of the strands in your plot
        :param pca_components: Default is 80. Has to be applied for better accuracy of t-SNE. You can indirectly change the described variance with this.
        :param model: The model you would like to choose for the embedding. Default is 'Rostlab/prot_t5_xl_half_uniref50-enc'. You can choose all models on https://huggingface.co/Rostlab and sgt.
        :param perplexity: Default is 30. It roughly determines the number of nearest neighbors that are considered in the embedding. A higher perplexity value results in a more global structure in the low-dimensional embedding, while a lower perplexity value emphasizes local structure. The optimal perplexity value for a given dataset depends on the dataset's intrinsic dimensionality, and it is usually determined by trial and err
        :param iterations_tsne: Default is 2500. number of times that the algorithm will repeat the optimization process for reducing the cost function. The optimization process aims to minimize the difference between the high-dimensional and low-dimensional representations of the data. More iterations result in a more optimized low-dimensional representation, but also increases the computational cost.
        :param batch_size: Default is 1000. The size of the sample which is chosen. The higher it is, the more computational intense.
        :param antigens: Datatype is a list. Is the names of the antigens in your binding data you would like to analyze. If you do not give any input, it will take the preferred antigen.
        :param show_antigen_names: Default is False. If you want to show the names of the antigens in the plot, set it to True.
        :param extra_figure: Default is False. If you want to have an extra figure with the indexes of the sequences in the vector space, set it to True.
        :return: Returns a plot where the sequences of the input samples are transformed in a vector space. Dimension reduction such as PCA and following t-SNE is used to plot it on a two dimensional space. The different colors indicate the different samples.
        """
        if samples == None:
            samples = [self.preferred_sample]
        self.ControlFigure.check_fig()
        self.ControlFigure.plot_type = "multi"
        assert self.binding_data is not None, "You have not given binding data. You can add it with the add_binding_data function"
        if antigens == None:
            antigens = [self.preferred_antigen]
        model_types_adv = ["Rostlab/ProstT5_fp16", "Rostlab/prot_t5_xl_uniref50", "Rostlab/prot_t5_base_mt_uniref50", "Rostlab/prot_bert_bfd_membrane", "Rostlab/prot_t5_xxl_uniref50", "Rostlab/ProstT5", "Rostlab/prot_t5_xl_half_uniref50-enc", "Rostlab/prot_bert_bfd_ss3", "Rostlab/prot_bert_bfd_localization", "Rostlab/prot_electra_generator_bfd", "Rostlab/prot_t5_xl_bfd", "Rostlab/prot_bert", "Rostlab/prot_xlnet", "Rostlab/prot_bert_bfd", "Rostlab/prot_t5_xxl_bfd"]
        models_all = model_types_adv + ["sgt"]
        assert model in models_all, f"Please enter a valid model name which are\n{models_all}. You can find nearly all of the models at: https://huggingface.co/Rostlab"
        if samples == None:
            samples = [self.experiments_list[0]]
        self.ControlFigure.check_fig()
        self.ControlFigure.plot_type = "single"
        incorrect_samples = [x for x in samples if x not in self.experiments_list]
        assert not incorrect_samples, f"The following sample(s) are not in your sequencing report: {', '.join(incorrect_samples)}. Please check the spelling or use the print_samples function to see the names of your samples"
        assert type(samples) == list, "You have to give a list with the samples you want to analyze"
        assert type(strands) == bool, "You have to give True or False as input for the strands parameter"
        assert type(pca_components) == int, "You have to give an integer as input for the pca_components"
        assert type(perplexity) == int, "You have to give an integer as input for the perplexity"
        assert type(iterations_tsne) == int, "You have to give an integer as input for the iterations_tsne"
        assert type(batch_size) == int, "You have to give an integer as input for the batch_size"
        assert batch_size > pca_components, "The batch_size has to be larger than the pca_components"
        assert batch_size > perplexity, "The batch_size has to be larger than the perplexity"
        assert pca_components > perplexity, "The pca_components has to be larger than the perplexity"
        assert type(extra_figure) == bool, "You have to give True or False as input for the extra_figure parameter"
        assert type(antigens) == list, "You have to give a list with the antigens you want to analyze"
        assert show_antigen_names in [True, False], "You have to give True or False as input for the show_antigen_names parameter"
        assert type(save_report_path) == str or save_report_path == None, "You have to give a string as input for the save_report_path parameter"

        self.ControlFigure.clear_fig()
        EmbeddingPlot = protein_embedding.Plot_Embedding(self.ControlFigure.ax, 
                                    self.sequencing_report,
                                    model,
                                    samples,
                                    strands,
                                    add_clone_size = 300,
                                    batch_size = batch_size,
                                    pca_components = pca_components,
                                    perplexity = perplexity,
                                    iterations_tsne = iterations_tsne,
                                    font_settings=self.font_settings,
                                    legend_settings = self.legend_settings,
                                    region_of_interest= self.region_of_interest,
                                    binding_data = self.binding_data,
                                    colorbar_settings = self.colorbar_settings,
                                    antigens = antigens, 
                                    toxin_names = show_antigen_names,
                                    extra_figure = extra_figure
                                    )
        
        self.ControlFigure.update_plot()
        self.style = plot_styler.PlotStyle(self.ControlFigure.ax,
                                           self.ControlFigure.plot_type)

        self.save_cluster_report(EmbeddingPlot.tsne_results, path = save_report_path)


    def cluster_one_AG(self, antigen=None, max_ld=1, min_ld=0, batch_size=1000, specific_experiments=False, preferred_cmap = "Blues", label_type = "numbers", save_report_path = None):
        """
        :param antigen: is the name of the antigen you would like to analyze
        :param max_ld: optional Parameter where its default is 1. Is the maximum Levenshtein distance you allow per cluster
        :param min_ld: optional Parameter where its default is 0. Is the minimum Levenshtein distance between sequences you allow
        :param batch_size: optional Parameter where its default is 1000. Is the batch size you want to use for the analysis
        :param specific_experiments: optional Parameter. You can give the names of specific samples in a list if you want
        :param preferred_cmap: optional Parameter. Define the colormap you would like to use.
        :param save_report_path: Default is None which saves your report in my_experiments/reports_pipeline. If you want to change it somewhere else you need to insert the full path with filename.
        :return: Creates a figure where sequences are clustered based on Levenshtein distance. Additionally the binding data of the sequences against a specific antigen is given.
        """

        self.ControlFigure.check_fig()
        self.ControlFigure.plot_type = "multi"
        assert self.binding_data is not None, "You have not given binding data. You can add it with the add_binding_data function"
        if antigen == None:
            antigen = self.preferred_antigen
        assert type(antigen) == str, "You have to give a string as input for the antigen"
        assert type(max_ld) == int, "You have to give an integer as input for the maximum levenshtein distance"
        assert type(min_ld) == int, "You have to give an integer as input for the minimum levenshtein distance"
        if specific_experiments != False:
            assert type(specific_experiments) == list, "You have to give a list with the samples you want to analyze"
        self.ControlFigure.clear_fig()
        cluster_report = levenshtein_clustering.cluster_single_AG(self.ControlFigure.fig,
                                                 self.sequencing_report,
                                                 antigen,
                                                 self.binding_data,
                                                 max_ld,
                                                 min_ld,
                                                 batch_size,
                                                 self.region_of_interest,
                                                 preferred_cmap,
                                                 specific_experiments,
                                                 label_type
                                                 )
        self.ControlFigure.update_plot()
        self.style = plot_styler.PlotStyle(self.ControlFigure.ax,
                                           self.ControlFigure.plot_type)
        self.save_cluster_report(cluster_report = cluster_report,
                                 path = save_report_path)

    def tsne_cluster_AG(self, sample=None, antigen=None, antigen_names=True, pca_components=70, perplexity=25,
                        iterations_tsne=2500, save_report_path = None):
        """
        :param sample: the sample you would like to analyze
        :param antigen: the toxins you would like to cluster
        :param antigen_names: Default is True. Prints the name of the toxin for the corresponding embedded sequence in the plot
        :param pca_components: optional. Default is 70
        :param perplexity: optional. Default 25
        :param iterations_tsne: optional. Default is 2500
        :return: It first embeds the sequences in a vector space and then reduces the dimensions and clusters them with PCA and TSNE. The sequences with the binding data are processed with the input sequences, to enable the plotting of the binding data.
        """
        if sample == None:
            sample = self.preferred_sample
        self.ControlFigure.check_fig()
        self.ControlFigure.plot_type = "multi"
        assert self.binding_data is not None, "You have not given binding data. You can add it with the add_binding_data function"
        if antigen == None:
            antigen = [self.preferred_antigen]
        assert type(sample) == str, "You have to give a string as input for the sample"
        assert sample in self.experiments_list, "The provided sample name is not in your sequencing report. Please check the spelling or use the print_samples function to see the names of your samples"
        assert antigen_names in [True, False], "You have to give True or False as input for the antigen_names parameter"
        assert type(pca_components) == int, "You have to give an integer as input for the pca_components"
        assert type(perplexity) == int, "You have to give an integer as input for the perplexity"
        assert type(iterations_tsne) == int, "You have to give an integer as input for the iterations_tsne"
        assert type(antigen) == list, "You have to give a list with the antigens you want to analyze"

        self.ControlFigure.clear_fig()
        tsne_results = embedding_with_binding.cluster_toxins_tsne(self.ControlFigure.fig,
                                                                  self.sequencing_report,
                                                                  sample,
                                                                  antigen,
                                                                  self.binding_data,
                                                                  antigen_names,
                                                                  pca_components,
                                                                  perplexity,
                                                                  iterations_tsne,
                                                                  self.font_settings,
                                                                  self.colorbar_settings,
                                                                  False,
                                                                  self.region_of_interest)
        self.ControlFigure.update_plot()
        self.style = plot_styler.PlotStyle(self.ControlFigure.ax,
                                           self.ControlFigure.plot_type)

        self.save_cluster_report(tsne_results, path = save_report_path)
        
    def embedding_network(self, samples = None,model = 'Rostlab/prot_t5_xl_half_uniref50-enc', batch_size = 200, cmap = "Blues", nodesize = 300, threshold_distance = 5):
        """
        :param samples: type is list. The samples you would like to compare towards their sequences
        :param model: The model you would like to choose for the embedding. You can choose "Rostlab/ProstT5_fp16",\n"Rostlab/prot_t5_xl_uniref50",\n"Rostlab/prot_t5_base_mt_uniref50",\n"Rostlab/prot_bert_bfd_membrane",\n"Rostlab/prot_t5_xxl_uniref50", \n"Rostlab/ProstT5", \n"Rostlab/prot_t5_xl_half_uniref50-enc", \n"Rostlab/prot_bert_bfd_ss3", \n"Rostlab/prot_bert_bfd_localization", \n"Rostlab/prot_electra_generator_bfd", \n"Rostlab/prot_t5_xl_bfd", \n"Rostlab/prot_bert", \n"Rostlab/prot_xlnet", \n"Rostlab/prot_bert_bfd", \n"Rostlab/prot_t5_xxl_bfd"
        :param batch_size: Default is 200. The size of the sample which is chosen. The higher it is, the more computational intense.
        :param nodesize: Default is 300. This value will be multiplied with the clone fractions of the sequences. If you set it to None, all nodes will have the same size.
        :threshold_distance: Default is 5. This value is multiplied with the summed batch_size of the samples. This defines the cluster size. The higher the values, the more complex the clusters.
        """
        model_types_adv = ["Rostlab/ProstT5_fp16", "Rostlab/prot_t5_xl_uniref50", "Rostlab/prot_t5_base_mt_uniref50", "Rostlab/prot_bert_bfd_membrane", "Rostlab/prot_t5_xxl_uniref50", "Rostlab/ProstT5", "Rostlab/prot_t5_xl_half_uniref50-enc", "Rostlab/prot_bert_bfd_ss3", "Rostlab/prot_bert_bfd_localization", "Rostlab/prot_electra_generator_bfd", "Rostlab/prot_t5_xl_bfd", "Rostlab/prot_bert", "Rostlab/prot_xlnet", "Rostlab/prot_bert_bfd", "Rostlab/prot_t5_xxl_bfd"]
        assert model in model_types_adv, f"Please enter a valid model name which are\n{model_types_adv}. You can find nearly all of the models at: https://huggingface.co/Rostlab"
        if samples == None:
            samples = [self.experiments_list[0]]
        self.ControlFigure.check_fig()
        self.ControlFigure.plot_type = "single"
        incorrect_samples = [x for x in samples if x not in self.experiments_list]
        assert type(samples) == list, "You have to give a list with the samples you want to analyze"
        assert type(batch_size) == int, "You have to give an integer as input for the batch_size"
        assert type(nodesize) == int, "Nodesize has to be an integer"
        assert nodesize > 0 & nodesize < 1000, "Please enter a value between 0 and 1000"
        assert type(threshold_distance) == int, "Please enter an integer for the threashold distance"
        if model in model_types_adv:
            protein_network_embedding.Network_Embedding(self.ControlFigure.ax,
                                                     self.sequencing_report,
                                                     samples,
                                                     model,
                                                     batch_size, 
                                                     self.font_settings,
                                                     cmap,
                                                     nodesize,
                                                     threshold_distance,
                                                     self.region_of_interest)
        
        self.ControlFigure.update_plot()
        self.style = plot_styler.PlotStyle(self.ControlFigure.ax,
                                           self.ControlFigure.plot_type)
        
        

    def embedding_tsne(self,
                       samples=None,
                       strands=True,
                       model = 'Rostlab/prot_t5_xl_half_uniref50-enc',
                       pca_components=80,
                       perplexity=30,
                       iterations_tsne=2500,
                       batch_size=1000):
        """
        :param samples: type is list. The samples you would like to compare towards their sequences
        :param strands: Default is True. It means that you will plot a batch of the strands in your plot
        :param pca_components: Default is 80. Has to be applied for better accuracy of t-SNE. You can indirectly change the described variance with this.
        :param model: The model you would like to choose for the embedding. Default is 'Rostlab/prot_t5_xl_half_uniref50-enc'. You can choose all models on https://huggingface.co/Rostlab and sgt.
        :param perplexity: Default is 30. It roughly determines the number of nearest neighbors that are considered in the embedding. A higher perplexity value results in a more global structure in the low-dimensional embedding, while a lower perplexity value emphasizes local structure. The optimal perplexity value for a given dataset depends on the dataset's intrinsic dimensionality, and it is usually determined by trial and err
        :param iterations_tsne: Default is 2500. number of times that the algorithm will repeat the optimization process for reducing the cost function. The optimization process aims to minimize the difference between the high-dimensional and low-dimensional representations of the data. More iterations result in a more optimized low-dimensional representation, but also increases the computational cost.
        :param batch_size: Default is 1000. The size of the sample which is chosen. The higher it is, the more computational intense.
        :return: Returns a plot where the sequences of the input samples are transformed in a vector space. Dimension reduction such as PCA and following t-SNE is used to plot it on a two dimensional space. The different colors indicate the different samples.
        """
        model_types_adv = ["Rostlab/ProstT5_fp16", "Rostlab/prot_t5_xl_uniref50", "Rostlab/prot_t5_base_mt_uniref50", "Rostlab/prot_bert_bfd_membrane", "Rostlab/prot_t5_xxl_uniref50", "Rostlab/ProstT5", "Rostlab/prot_t5_xl_half_uniref50-enc", "Rostlab/prot_bert_bfd_ss3", "Rostlab/prot_bert_bfd_localization", "Rostlab/prot_electra_generator_bfd", "Rostlab/prot_t5_xl_bfd", "Rostlab/prot_bert", "Rostlab/prot_xlnet", "Rostlab/prot_bert_bfd", "Rostlab/prot_t5_xxl_bfd"]
        models_all = model_types_adv + ["sgt"]
        assert model in models_all, f"Please enter a valid model name which are\n{models_all}. You can find nearly all of the models at: https://huggingface.co/Rostlab"
        if samples == None:
            samples = [self.experiments_list[0]]
        self.ControlFigure.check_fig()
        self.ControlFigure.plot_type = "single"
        incorrect_samples = [x for x in samples if x not in self.experiments_list]
        assert not incorrect_samples, f"The following sample(s) are not in your sequencing report: {', '.join(incorrect_samples)}. Please check the spelling or use the print_samples function to see the names of your samples"
        assert type(samples) == list, "You have to give a list with the samples you want to analyze"
        assert type(strands) == bool, "You have to give True or False as input for the strands parameter"
        assert type(pca_components) == int, "You have to give an integer as input for the pca_components"
        assert type(perplexity) == int, "You have to give an integer as input for the perplexity"
        assert type(iterations_tsne) == int, "You have to give an integer as input for the iterations_tsne"
        assert type(batch_size) == int, "You have to give an integer as input for the batch_size"
        assert batch_size > pca_components, "The batch_size has to be larger than the pca_components"
        assert batch_size > perplexity, "The batch_size has to be larger than the perplexity"
        self.ControlFigure.clear_fig()
        if model in model_types_adv:
            protein_embedding.Plot_Embedding(self.ControlFigure.ax, 
                                              self.sequencing_report,
                                              model,
                                              samples,
                                              strands,
                                              add_clone_size = 300,
                                              batch_size = batch_size,
                                              pca_components = pca_components,
                                              perplexity = perplexity,
                                              iterations_tsne = iterations_tsne,
                                                font_settings=self.font_settings,
                                              legend_settings = self.legend_settings,
                                              region_of_interest= self.region_of_interest
                                              )
        elif model == "sgt":
            cluster_embedding.show_difference(self.sequencing_report,
                                            samples,
                                            strands,
                                            batch_size,
                                            pca_components,
                                            perplexity,
                                            iterations_tsne,
                                            self.region_of_interest,
                                            self.ControlFigure.ax,
                                            self.legend_settings,
                                            self.font_settings)
        else: 
            return
        self.ControlFigure.update_plot()
        self.style = plot_styler.PlotStyle(self.ControlFigure.ax,
                                           self.ControlFigure.plot_type)

    def morosita_horn(self,cmap = "Blues", annotate_cells=False, specific_experiments=False, matrix_save_path = None):
        """
        :param annotate_cells: Default is False. If you want to see the values of the matrix, set it to True.
        :param specific_experiments: you can give a list with specific experiments
        :param matrix_save_path: Path where you want to save the matrix
        :return: Returns a matrix of the identity between your samples based on the Morosita Horn Index
        """
        self.ControlFigure.check_fig()
        self.ControlFigure.plot_type = "single"
        if specific_experiments != False:
            assert type(specific_experiments) == list, "You have to give a list with the samples you want to analyze"
            incorrect_samples = [x for x in specific_experiments if x not in self.experiments_list]
            assert not incorrect_samples, f"The following sample(s) are not in your sequencing report: {', '.join(incorrect_samples)}. Please check the spelling or use the print_samples function to see the names of your samples"
        self.ControlFigure.clear_fig()
        matrix = plt_heatmap.plot_heatmap(self.sequencing_report,
                                          True,
                                          "morosita_horn",
                                          self.ControlFigure.ax,
                                          self.colorbar_settings,
                                          self.font_settings,
                                          annotate_cells,
                                          self.region_of_interest,
                                          cmap,
                                          specific_experiments=specific_experiments,
                                          )
        self.ControlFigure.update_plot()
        self.style = plot_styler.PlotStyle(self.ControlFigure.ax,
                                           self.ControlFigure.plot_type)
        save_matrix(matrix, matrix_save_path)

    def jaccard(self,cmap = "Blues",  annotate_cells=False, specific_experiments=False, matrix_save_path = None):
        """
        :param annotate_cells: Default is False. If you want to see the values of the matrix, set it to True.
        :param specific_experiments: give a list with specific samples you would like to analyze
        :return: Returns a matrix of the identity between your samples based on the Jaccard Index
        """
        self.ControlFigure.check_fig()
        self.ControlFigure.plot_type = "single"
        if specific_experiments != False:
            assert type(specific_experiments) == list, "You have to give a list with the samples you want to analyze"
            incorrect_samples = [x for x in specific_experiments if x not in self.experiments_list]
            assert not incorrect_samples, f"The following sample(s) are not in your sequencing report: {', '.join(incorrect_samples)}. Please check the spelling or use the print_samples function to see the names of your samples"
        self.ControlFigure.clear_fig()
        matrix = plt_heatmap.plot_heatmap(self.sequencing_report,
                                          True,
                                          "jaccard",
                                          self.ControlFigure.ax,
                                          self.colorbar_settings,
                                          self.font_settings,
                                          annotate_cells,
                                          self.region_of_interest,
                                          cmap,
                                          specific_experiments=specific_experiments,
                                          )
        self.ControlFigure.update_plot()
        self.style = plot_styler.PlotStyle(self.ControlFigure.ax,
                                           self.ControlFigure.plot_type)

        save_matrix(matrix, matrix_save_path)

    def sample_diversity(self, method = "InverseSimpson"):
        """
        :param method: Default is InverseSimpson. You can choose between InverseSimpson and Shannon which are two available indices to measure the diversity of the sequences in each of your samples based on the clone fraction.
        :return: Returns a barplot where the diversity of the sequences in each of your samples is shown.
        """
        assert method in ["InverseSimpson", "Shannon"], "Please enter a valid method. You can choose between InverseSimpson and Shannon"
        self.ControlFigure.check_fig()
        self.ControlFigure.plot_type = "single"
        self.ControlFigure.clear_fig()
        diversity_plot.DiversityPlot(self.sequencing_report,
                                     self.ControlFigure.ax,
                                     self.region_of_interest, 
                                     self.font_settings,
                                     method = method)
        self.ControlFigure.update_plot()
        self.style = plot_styler.PlotStyle(self.ControlFigure.ax,
                                           self.ControlFigure.plot_type)

    def sorensen(self,cmap = "Blues",  annotate_cells=False, specific_experiments=False, matrix_save_path = None):
        """
        :param annotate_cells: Default is False. If you want to see the values of the matrix, set it to True.
        :param specific_samples: give a list with specific samples you would like to analyze
        :return: Returns a matrix of the identity between your samples based on the Sorensen Dice Index
        """
        self.ControlFigure.check_fig()
        self.ControlFigure.plot_type = "single"
        if specific_experiments != False:
            assert type(specific_experiments) == list, "You have to give a list with the samples you want to analyze"
            incorrect_samples = [x for x in specific_experiments if x not in self.experiments_list]
            assert not incorrect_samples, f"The following sample(s) are not in your sequencing report: {', '.join(incorrect_samples)}. Please check the spelling or use the print_samples function to see the names of your samples"
        self.ControlFigure.clear_fig()
        matrix = plt_heatmap.plot_heatmap(self.sequencing_report,
                                          True,
                                          "sorensen",
                                          self.ControlFigure.ax,
                                          self.colorbar_settings,
                                          self.font_settings,
                                          annotate_cells,
                                          self.region_of_interest,
                                          cmap,
                                          specific_experiments=specific_experiments,
                                          )
        self.ControlFigure.update_plot()
        self.style = plot_styler.PlotStyle(self.ControlFigure.ax,
                                           self.ControlFigure.plot_type)
        save_matrix(matrix, matrix_save_path)

    def relative(self, cmap = "Blues",  annotate_cells=False, specific_experiments=False):
        """
        :param annotate_cells: Default is False. If you want to see the values of the matrix, set it to True.
        :param specific_samples: give a list with specific samples you would like to analyze
        :return: Returns a matrix of the identity between your samples based on the proportion of identical sequences
        """
        self.ControlFigure.check_fig()
        self.ControlFigure.plot_type = "single"
        if specific_experiments != False:
            assert type(specific_experiments) == list, "You have to give a list with the samples you want to analyze"
            incorrect_samples = [x for x in specific_experiments if x not in self.experiments_list]
            assert not incorrect_samples, f"The following sample(s) are not in your sequencing report: {', '.join(incorrect_samples)}. Please check the spelling or use the print_samples function to see the names of your samples"
        self.ControlFigure.clear_fig()
        matrix = plt_heatmap.plot_heatmap(self.sequencing_report,
                                          True,
                                          "relative",
                                          self.ControlFigure.ax,
                                          self.colorbar_settings,
                                          self.font_settings,
                                          annotate_cells,
                                          self.region_of_interest,
                                          cmap,
                                          specific_experiments=specific_experiments,
                                          )
        self.ControlFigure.update_plot()
        self.style = plot_styler.PlotStyle(self.ControlFigure.ax,
                                           self.ControlFigure.plot_type)
        save_matrix(matrix)

    def levenshtein_dendrogram(self, sample=None, max_cluster_dist=2, batch_size=1000):
        """
        :params sample: the sample you would like to analyze
        :max_cluster_dist: Default is 2. Maximum levenshtein distance between sequences within a cluster.
        :params batch_size: Default is 1000. The size of the sample which is chosen.
        :return: Returns a dendrogram of the sequences based on the Levenshtein distance
        """
        if sample == None:
            sample = self.preferred_sample
        assert type(sample) == str, "You have to give a string as input for the sample"
        assert sample in self.experiments_list, "The provided sample name is not in your sequencing report. Please check the spelling or use the print_samples function to see the names of your samples"
        assert type(
            max_cluster_dist) == int, "You have to give an integer as input for the maximum levenshtein distance"
        assert type(batch_size) == int, "You have to give an integer as input for the batch size"
        self.ControlFigure.check_fig()
        self.ControlFigure.plot_type = "single"
        self.ControlFigure.clear_fig()
        linked  = hist_lvst_dist.levenshtein_dend(self.ControlFigure.ax,
                                        self.sequencing_report,
                                        sample,
                                        batch_size,
                                        max_cluster_dist,
                                        self.font_settings,
                                        self.region_of_interest
                                        )
        self.ControlFigure.update_plot()
        self.style = plot_styler.PlotStyle(self.ControlFigure.ax,
                                           self.ControlFigure.plot_type)
        

    def dendro_bind(self, sample=None, antigens=None, max_cluster_dist=2, batch_size=1000, ascending=True):
        """
        :params sample: the sample you would like to analyze
        :params antigens: the antigens you would like to analyze. The input is a list.
        :max_cluster_dist: Default is 2. Maximum levenshtein distance between sequences within a cluster. The higher this number the bigger the dendrogram
        :params batch_size: Default is 1000. The number of sequences starting with the highest fractions which are used for the analysis from the sample
        :params ascending: Default is True. If you want to see the highest fractions first, set it to False
        :return: Creates a dendrogram which shows the realtionship between your sequences in the sample and your sequences with binding data based on levenshtein distance. Additionally a barplot with the binding data of the found sequences is given.
        """
        if sample == None:
            sample = self.preferred_sample
        assert self.binding_data is not None, "You have not given binding data. You can add it with the add_binding_data function"
        if antigens == None:
            antigens = [self.preferred_antigen]
        assert type(sample) == str, "You have to give a string as input for the sample"
        assert sample in self.experiments_list, "The provided sample name is not in your sequencing report. Please check the spelling or use the print_samples function to see the names of your samples"
        assert type(
            max_cluster_dist) == int, "You have to give an integer as input for the maximum levenshtein distance"
        assert type(batch_size) == int, "You have to give an integer as input for the batch size"
        assert type(ascending) == bool, "You have to give True or False as input for the ascending parameter"

        self.ControlFigure.check_fig()
        self.ControlFigure.plot_type = "multi"
        self.ControlFigure.clear_fig()
        fig2 = hist_lvst_dist.dendo_binding(self.ControlFigure.fig,
                                     self.sequencing_report,
                                     self.binding_data,
                                     sample,
                                     antigens,
                                     batch_size,
                                     max_cluster_dist,
                                     self.font_settings,
                                     self.region_of_interest,
                                     ascending
                                     )
        self.ControlFigure.update_plot()
        self.style = plot_styler.PlotStyle(self.ControlFigure.ax,
                                           self.ControlFigure.plot_type)
      #  self.ControlFigure.tighten()
        return fig2

    def connect_samples(self, summed_clonefraction = 0.5, max_num_reads = 100,color_lines = "peachpuff", max_weight_lines = 100 ):
        """
        :param summed_clonefraction: You will take the top 50% of clones per sample per default. You can change that to take more or less crones
        :param max_num_reads: You will limit the maximum number of clones to 100 per default. If you want to have more your can increase that
        :param max_weight_lines: If you have a levenshtein distance of 0 between two sequences this value will be the thickness of the line. The thickness of other distances is indicated through: max_weight_lines * (1 / ls_distance ** 2)
        :result: A network plot which shows the relation between the samples for the chosen clones based on levenshtein distance. The thickness of the lines indicates the inverse levenshtein distance.
        """
        self.ControlFigure.check_fig()
        self.ControlFigure.plot_type = "single"
        self.ControlFigure.clear_fig()
        sample_cluster.ClusterExperiment(self.sequencing_report,
                                         self.ControlFigure.ax,
                                         self.region_of_interest,
                                         summed_clonefraction=summed_clonefraction,
                                         max_num_reads=max_num_reads,
                                         edge_color = color_lines,
                                         max_weight_lines=max_weight_lines)
        self.ControlFigure.update_plot()
        self.style = plot_styler.PlotStyle(self.ControlFigure.ax,
                                           self.ControlFigure.plot_type)

    def validate_mixcr_path(self):
        path_to_mixcr = self.global_params["mixcr_path"]
        settings_dir = os.path.join(self.module_dir, "settings", "global_vars.txt")
        check_mixcr(path_to_mixcr, self.global_params, settings_dir)

    def create_parser(self):
        self.validate_mixcr_path()
        path_to_mixcr = self.global_params["mixcr_path"]
        commands = []
        commands.extend(["java", f"-Xms{self.java_heap_size}M", "-jar"])  # enable change of para
        commands.extend([path_to_mixcr])
        return commands
    

    def export_mixcr_quality(self,save_dir = None, specific_chain=None, metric=None, plot_type=None):
        if save_dir == None:
            save_dir = self.mixcr_plots_path
        output_file = os.path.join(save_dir, self.experiment + "_alignQc.pdf")
        path_to_tables = os.path.join(self.module_dir, "my_experiments", self.experiment, "clones_result", "*.clns")
        export_plots_commands = self.create_parser()
        
        export_plots_commands.extend(["exportQc"])
        export_plots_commands.extend(["align"])
        export_plots_commands.extend([path_to_tables])
        export_plots_commands.extend([output_file])
        export_plots_commands.extend(["--force-overwrite"])
        subprocess.run(export_plots_commands)
        print(f"You can find the result file at: {output_file}")

    #    for i in glob.glob(path_to_tables):

     #       name=os.path.basename(i).split(".")[0]
      #      output_file = os.path.join(save_dir, name + "_tags.pdf")
       #     export_plots_commands = self.create_parser()
        #    export_plots_commands.extend(["exportQc"])
         #   export_plots_commands.extend(["tags"])
         #   export_plots_commands.extend([i])
         #   export_plots_commands.extend([output_file])
         #   export_plots_commands.extend(["--force-overwrite"])
          #  subprocess.run(export_plots_commands)
            #print(f"You can find the result file at: {output_file}")
        #    name = os.path.basename(i).split(".")[0]
         #   output_file = os.path.join(save_dir, name + "chainUsage.png")
          #  export_plots_commands = self.create_parser()
          #  export_plots_commands.extend(["exportQc"])
           # export_plots_commands.extend(["chainUsage"])
           # export_plots_commands.extend(["--hide-non-functional"])
           # export_plots_commands.extend([i])
           # export_plots_commands.extend([output_file])
           # export_plots_commands.extend(["--force-overwrite"])
           # subprocess.run(export_plots_commands)
           # print(f"You can find the result file at: {output_file}")
        
        output_file = os.path.join(save_dir, self.experiment + "_coverage.png")
        export_plots_commands = self.create_parser()
        export_plots_commands.extend(["exportQc"])
        export_plots_commands.extend(["coverage"])
        export_plots_commands.extend([path_to_tables])
        export_plots_commands.extend([output_file])
        export_plots_commands.extend(["--force-overwrite"])
        subprocess.run(export_plots_commands)
        print(f"You can find the result file at: {output_file}")

        if specific_chain != None:
            export_plots_commands.extend(["--chains", specific_chain])
        if metric != None:
            export_plots_commands.extend(["--metrics", metric])
        if plot_type != None:
            export_plots_commands.extend(["--plotType", plot_type])


    def mixcr_explain_diversity(self,plot_name = "cdr3_metrics", plot_save_dir = None, metric = "cdr3metrics",chains = None, plot_type = "boxplot", output_type = "pdf"):
        """
        :param: metric: Default is cdr3metrics. You can choose between cdr3metrics and diversity. If you choose diversity you will focus on the clone fractions whereas with cdr3metrics you will focus on certain aspects of these chains such as length (https://mixcr.com/mixcr/reference/mixcr-postanalysis/#diversity-measures)
        :param: chains: Default is None. Possible inputs are : IGH, IGL, IGK, TRA, TRB, TRG, TRD, IG. You will choose a specific chain for the analysis.
        :param: plot_type: Default is violin. You can choose between: boxplot, boxplot-bindot, boxplot-jitter, violin, violin-bindot, barplot, barplot-stacked, lineplot, lineplot-jitter, lineplot-bindot, scatter
        :param: output_type: Default is pdf. Allows you to control the file type of your output. Possible values are: pdf, jpeg, png, svg, eps,
        """
        if metric in ["cdr3metrics", "diversity"]:
            pass
        else:
            print("Please choose a correct value for the metric. Pipeline continues with:\ncdr3metrics and cdr3lenAA for metric_type")
            metric = "cdr3metrics"
            metric_type = "cdr3lenAA"

        chains_values = [None, "IGH", "IGL", "IGK", "TRA", "TRB", "TRG", "TRD", "IG"]
        if chains not in chains_values:
            print("Please enter a correct shortcut for the chain preferred chain. These are: ")
            [print(i) for i in chains_values]
        plot_type_values = ["boxplot", "boxplot-bindot", "boxplot-jitter", "violin", "violin-bindot", "barplot", "barplot-stacked", "lineplot", "lineplot-jitter", "lineplot-bindot", "scatter"]
        if plot_type not in plot_type_values:
            print("Please enter a correct shortcut for the chain plot type. These are: ")
            [print(i) for i in plot_type_values]


        print(f"To be able to understand the following processing go to: https://mixcr.com/mixcr/reference/mixcr-postanalysis/#overlap-postanalysis")
        plot_name = plot_name

        if plot_save_dir == None:
            plot_save_dir = self.mixcr_plots_path
        save_dir = os.path.join(self.module_dir, "my_experiments", self.experiment)
        json_dir = os.path.join(save_dir, "pa_individual.json")
        path_to_tables = os.path.join(self.module_dir, "my_experiments", self.experiment, "clones_result", "*.clns")
        if not os.path.isfile(json_dir):
            export_plots_commands = self.create_parser()
            export_plots_commands.extend(["postanalysis"])
            export_plots_commands.extend(["individual"])
            export_plots_commands.extend(["--default-downsampling", "none"])
            export_plots_commands.extend(["--default-weight-function", "none"])
            export_plots_commands.extend([path_to_tables])
            export_plots_commands.extend([json_dir])
            subprocess.run(export_plots_commands)
        export_plots_commands = self.create_parser()
        export_plots_commands.extend(["exportPlots"])
        export_plots_commands.extend([metric])
      #  export_plots_commands.extend(["--metric",metric_type])
      #  export_plots_commands.extend(["-f"])
        if not chains == None:
            export_plots_commands.extend(["--chains", chains])
        export_plots_commands.extend(["--plot-type", plot_type])
        export_plots_commands.extend([json_dir])
        export_plots_commands.extend([os.path.join(plot_save_dir, plot_name + "." + output_type)])
        subprocess.run(export_plots_commands)



    def mixcr_vdj_usage(self,):
        save_dir = self.mixcr_plots_path

        json_dir = os.path.join(self.module_dir, "my_experiments", self.experiment, "pa_individual.json")
        path_to_tables = os.path.join(self.module_dir, "my_experiments", self.experiment, "clones_result", "*.clns")
        if not os.path.isfile(json_dir):
            export_plots_commands = self.create_parser()
            export_plots_commands.extend(["postanalysis"])
            export_plots_commands.extend(["individual"])
            export_plots_commands.extend(["--default-downsampling", "none"])
            export_plots_commands.extend(["--default-weight-function", "none"])
            export_plots_commands.extend([path_to_tables])
            export_plots_commands.extend([json_dir])
            subprocess.run(export_plots_commands)

        export_plots_commands = self.create_parser()
        export_plots_commands.extend(["exportPlots"])
        export_plots_commands.extend(["vUsage"])
      #  export_plots_commands.extend(["-f"])
        export_plots_commands.extend([json_dir])
        export_plots_commands.extend([os.path.join(save_dir, "vUsage_heatmap.png")])
        subprocess.run(export_plots_commands)
        export_plots_commands = self.create_parser()
        export_plots_commands.extend(["exportPlots"])
        export_plots_commands.extend(["jUsage"])
       # export_plots_commands.extend(["-f"])
        export_plots_commands.extend([json_dir])
        export_plots_commands.extend([os.path.join(save_dir, "jUsage_heatmap.png")])
        subprocess.run(export_plots_commands)
        export_plots_commands = self.create_parser()
        export_plots_commands.extend(["exportPlots"])
        export_plots_commands.extend(["isotypeUsage"])
       # export_plots_commands.extend(["-f"])
        export_plots_commands.extend([json_dir])
        export_plots_commands.extend([os.path.join(save_dir, "isotypeUsage_heatmap.png")])
        subprocess.run(export_plots_commands)


    def mixcr_overlap(self):
        save_dir = self.plot_path
        path_to_tables = os.path.join(self.module_dir, "my_experiments", self.experiment, "clones_result")
        json_dir = os.path.join(self.module_dir, "my_experiments", self.experiment, "pa_overlap.json")
        if not os.path.isfile(json_dir):
            export_plots_commands = self.create_parser()
            export_plots_commands.extend(["postanalysis"])
            export_plots_commands.extend(["overlap"])
            export_plots_commands.extend(["--default-downsampling", "none"])
            export_plots_commands.extend(["--default-weight-function", "none"])
            export_plots_commands.extend([path_to_tables])
            export_plots_commands.extend([json_dir])
            subprocess.run(export_plots_commands)
        possible_metrics = ["SharedClonotypes", "RelativeDiversity", "F1Index", "F2Index", "JaccardIndex", "Pearson", "PearsonAll"]

        for metric in possible_metrics:
            export_plots_commands = self.create_parser()
            export_plots_commands.extend(["exportPlots"])
            export_plots_commands.extend(["overlap"])
            export_plots_commands.extend(["--metric", metric])
            export_plots_commands.extend(["--color-key", "Patient"])
            export_plots_commands.extend([json_dir])
            export_plots_commands.extend([os.path.join(save_dir, f"{metric}.png")])
            subprocess.run(export_plots_commands)
       # export_plots_commands = self.create_parser()
       # export_plots_commands.extend(["exportPlots"])
       # export_plots_commands.extend(["overlap"])
   #     export_plots_commands.extend(["--metric", metric])
       # export_plots_commands.extend(["--color-key", "Patient"])
       # export_plots_commands.extend([json_dir])
       # export_plots_commands.extend([os.path.join(save_dir, f"overlap.png")])
       # subprocess.run(export_plots_commands)





def automation(make_report = False):
    plot = PlotManager()
    regions = plot.avail_regions
    if "targetSequences" in regions:
        regions.remove("targetSequences")
    for region in regions:
        plot.change_region(region = region)
        plot.full_analysis()
        if make_report == True:
            plot.create_report()
    
        