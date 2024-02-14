from textwrap import wrap
from ExpoSeq.plots.tidy_protbert_embedding import TransformerBased
import pandas as pd
import matplotlib.pyplot as plt
import warnings


class PrepareData:
    
    @staticmethod
    def logical_check(batch_size, perplexity, pca_components, model, iterations_tsne):
        # in the pipeline is another test, where the sample name is checked which has to be in the sequencing report 
        assert batch_size > perplexity, "The batch_size has to be larger than the perplexity"
        assert batch_size > pca_components, "The batch_size has to be larger than the pca_components"
        assert batch_size != 0, "batch_size value must not be 0" # indirectly tests that perplexity and pca components must not be 0 as well
        # more models will follow
        models_all = ["Rostlab/ProstT5_fp16", "Rostlab/prot_t5_xl_uniref50", "Rostlab/prot_t5_base_mt_uniref50", "Rostlab/prot_bert_bfd_membrane", "Rostlab/prot_t5_xxl_uniref50", "Rostlab/ProstT5", "Rostlab/prot_t5_xl_half_uniref50-enc", "Rostlab/prot_bert_bfd_ss3", "Rostlab/prot_bert_bfd_localization", "Rostlab/prot_electra_generator_bfd", "Rostlab/prot_t5_xl_bfd", "Rostlab/prot_bert", "Rostlab/prot_xlnet", "Rostlab/prot_bert_bfd", "Rostlab/prot_t5_xxl_bfd"]
        assert model in models_all, f"Please enter a valid model name which are\n{models_all}. You can find the models at: https://huggingface.co/Rostlab"
        assert iterations_tsne != 0, "the number of iterations for the tsne algorithm must not be 0."
        
    @staticmethod
    def datatype_check(samples, pca_components, perplexity, iterations_tsne, batch_size, model):
        assert type(samples) == list, "You have to give a list with the samples you want to analyze"
        assert type(pca_components) == int, "You have to give an integer as input for the pca_components"
        assert type(perplexity) == int, "You have to give an integer as input for the perplexity"
        assert type(iterations_tsne) == int, "You have to give an integer as input for the iterations_tsne"
        assert type(batch_size) == int, "You have to give an integer as input for the batch_size"
        assert type(model) == str, "The input for model must be a string"
           
    
    @staticmethod
    def check_warnings(sequences, pca_components, perplexity):
        if len(sequences) < pca_components:
            warnings.warn(f"The number of sequences you have is {len(sequences)} but you need to have more sequences than principal components which is: {pca_components}") 
            print(f"Number of principal components is set to number of sequences ({len(sequences)})")
            pca_components = len(sequences) 
        if pca_components < perplexity:
            warnings.warn("The number of reduced dimensions you have is " + str(pca_components) + "but you need to have more than perplexity which is: " + str(perplexity)) 
            print(f"Perplexity is set to the half of reduced dimensions ({pca_components//2})")
            perplexity = pca_components//2
            if perplexity < 1:
                perplexity = 1
        return pca_components, perplexity
    
    @staticmethod
    def filter_binding_data(binding_data, region_of_interest, antigens):
        if binding_data is not None:
            merged_columns = [region_of_interest] + antigens
            binding_data = binding_data[merged_columns]
    
    def tidy(self, sequencing_report, list_experiments, region_of_interest, antigens = None, batch_size = 300, pca_components = 70,
             perplexity = 25, iterations_tsne = 1000, model_choice = "Rostlab/prot_bert",binding_data = None,cf_column_name = "cloneFraction", sample_column_name = "Experiment"):
        for sample in list_experiments:
            assert sample in list(sequencing_report[sample_column_name].unique()), f"{sample} does not exist"
        self.logical_check(batch_size,
                           perplexity,
                           pca_components,
                           model_choice,
                           iterations_tsne)
        self.datatype_check(list_experiments,
                            pca_components,
                            perplexity,
                            iterations_tsne, 
                            batch_size,
                            model_choice)
        self.filter_binding_data(binding_data,
                                 region_of_interest, 
                                 antigens)
        Transformer = TransformerBased(choice = model_choice)
        sequences,sequences_filtered, selected_rows = Transformer.filter_sequences(sequencing_report,
                                                                               batch_size,
                                                                               list_experiments,
                                                                               binding_data,
                                                                               region_of_interest=region_of_interest, 
                                                                               cf_column_name=cf_column_name,
                                                                               sample_column_name=sample_column_name)
        pca_components, perplexity = self.check_warnings(sequences, pca_components, perplexity)
        X = Transformer.do_pca(sequences, batch_size, pca_components)
        peptides = selected_rows[region_of_interest].to_list()
        self.clones = selected_rows[cf_column_name]
        tsne_results = Transformer.do_tsne(X, perplexity, iterations_tsne)
        return peptides, selected_rows, tsne_results
    
    @staticmethod
    def make_csv(selected_rows, tsne_results):
        selected_rows["tsne1"] = tsne_results["tsne1"]
        selected_rows["tsne2"] = tsne_results["tsne2"]
        selected_rows.to_csv("tsne_results.csv")
        



class PlotEmbedding:
    def __init__(self,sequencing_report, model_choice, list_experiments, region_of_interest, strands,add_clone_size, batch_size, pca_components, 
                 perplexity, iterations_tsne, antigens = None, font_settings = {}, legend_settings = {}, ax = None, binding_data = None, # antigens mutual attribute 
                 colorbar_settings = None,  toxin_names = None, extra_figure = False):
        self.ax = ax
        self.binding_data = binding_data
        self.data_prep  = PrepareData()
        peptides, selected_rows, tsne_results = self.data_prep.tidy(sequencing_report=sequencing_report,
                                                            list_experiments=list_experiments,
                                                            region_of_interest=region_of_interest,
                                                            antigens = antigens,
                                                            batch_size = batch_size,
                                                            pca_components=pca_components,
                                                            perplexity=perplexity,
                                                            iterations_tsne=iterations_tsne,
                                                            model_choice=model_choice,
                                                            binding_data=binding_data)
        if self.ax != None:
            if self.binding_data is not None:
                self.create_binding_plot(selected_rows, antigens, region_of_interest, toxin_names, colorbar_settings)
                if extra_figure == True and font_settings != {}:
                    self.create_second_bind_plot(font_settings)
                    title = "\n".join(wrap(f"t-SNE embedding for {antigens}", 40))
                    self.ax.set_title(title, pad= 12, **font_settings)
            else:
                self.create_plot(tsne_results, selected_rows)
                if legend_settings != {}:
                    self.add_legend(list_experiments, legend_settings)
                title = "\n".join(wrap("t-SNE embedding for given samples", 40))
                if font_settings != {}:
                    self.ax.set_title(title, pad= 12, **font_settings)
            if font_settings != {}:
                    self.ax.set_xlabel("t-SNE1", **font_settings) # add font_settings
            self.ax.set_ylabel("t-SNE2", **font_settings)

            if strands == True:
                self.add_seq_anotation(peptides, tsne_results)
            if add_clone_size != None:
                self.add_size(add_clone_size)
            
    
        
    def create_plot(self, tsne_results, selected_rows):
        experiments_batch = selected_rows["Experiment"]
        self.tsne_plot = self.ax.scatter(tsne_results.tsne1,
                                    tsne_results.tsne2,
                                    c = pd.factorize(experiments_batch)[0],
                                    alpha = 0.5,
                                    )


            
    def add_size(self, add_clone_size):
        size_points = self.data_prep.clones * add_clone_size
        self.tsne_plot.set_sizes(size_points)
        
    def add_legend(self, list_experiments, legend_settings):
        self.ax.legend(handles=self.tsne_plot.legend_elements()[0],
            labels=list_experiments,
            **legend_settings)

        

    def add_seq_anotation(self,peptides, tsne_results):
        x = tsne_results["tsne1"].values.tolist()
        y = tsne_results["tsne2"].values.tolist()
        for i in range(0, len(x), 10):
            self.ax.annotate(peptides[i],
                        (x[i][0], y[i][0]),
                        fontsize = 5,
                        )
            

    def return_binding_results(self,  selected_rows, antigens, region_of_interest):
        kds = selected_rows[antigens].max(axis = 1) # if there are multiple values for the same sequence this will find the highest one 
        ids = selected_rows[antigens].idxmax(axis = 1)
        aminoacids = selected_rows[region_of_interest].to_list()
        experiments_batch = selected_rows["Experiment"]
        unique_experiments_num = pd.factorize(experiments_batch)[0]
        self.tsne_results["experiments_factorized"] = list(unique_experiments_num)
        self.tsne_results["experiments_string"] = list(experiments_batch)
        self.tsne_results["binding"] = list(kds)
        self.tsne_results["sequences"] = list(aminoacids)
        self.tsne_results['sequence_id'] = pd.Series(range(self.tsne_results.shape[0]))
        return kds, ids
    
    def create_second_bind_plot(self, font_settings):

        fig2 = plt.figure()
        ax2 = fig2.gca()
        ax2.scatter(self.tsne_results.tsne1, self.tsne_results.tsne2, alpha = 0.0)
        ax2.set_xlabel('t-SNE 1', **font_settings)
        ax2.set_ylabel('t-SNE 2', **font_settings)
        n = 0
        for j, row in self.tsne_results.iterrows():
            if row["binding"] > 1:
                    ax2.text(row['tsne1'], row['tsne2'], row['sequence_id'], fontsize=10, weight = "bold")

            else:
                if n == 6:
                    n = 0
                    ax2.text(row['tsne1'], row['tsne2'], row['sequence_id'], fontsize=8)
            n += 1
    
    def create_binding_plot(self, selected_rows, antigens, region_of_interest, toxin_names, colorbar_settings):
        kds, ids = self.return_binding_results(selected_rows, antigens, region_of_interest)
        self.tsne_plot = self.ax.scatter(self.tsne_results.tsne1,
                            self.tsne_results.tsne2,
                            c = self.tsne_results.binding,
                            alpha = 1,
                            cmap = "magma"
                            )
        x_cor = list(self.tsne_results.tsne1.iloc[:, 0])
        y_cor = list(self.tsne_results.tsne2.iloc[:, 0])
        if toxin_names == True:
            for i, txt in enumerate(list(ids)):
                if list(kds)[i] > 0:
                    self.ax.annotate(txt, (x_cor[i], y_cor[i]))
        else:
            pass
        plt.colorbar(self.tsne_plot, **colorbar_settings)
        


#sequencing_report_path = r"src/ExpoSeq/software_tests/test_files/test_show/sequencing_report.csv"
#sequencing_report = pd.read_csv(sequencing_report_path)
#sequencing_report["cloneFraction"] = sequencing_report["readFraction"]
#peptides, selected_rows, tsne_results = PrepareData().tidy(sequencing_report, ["GeneMind_1"], region_of_interest = "aaSeqCDR3", batch_size = 50)
#print(peptides)
#print(selected_rows)
#print(tsne_results)