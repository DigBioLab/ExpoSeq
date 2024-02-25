from textwrap import wrap
from .tidy_protbert_embedding import TransformerBased
import pandas as pd
import matplotlib.pyplot as plt
import warnings
import peptides
from .contents import hytdrophobicity_metrics, isoelectric_point_scales

class GetProteinProperty:

    def __init__(self, sequences):
        """This class allows you to return a dictionary with the sequences and the sequence property of interest

        Args:
            sequences (list): list of peptide sequences
        """
        self.sequences = sequences
        self.sequence_property_interest = {}
        self.attribute_funcs = {"aliphatic_index": self.get_aliphatic_index,
        "isoelectric_point": self.get_isoelectric_point,
        "hydrophobicity": self.get_hydrophobicity,
        "weight": self.get_weight,
        "mass_charge_ratio": self.mass_charge_ratio,
        "length": self.calc_len,
     }
        
    @staticmethod
    def calc_len(peptide_object):
        return len(peptide_object.sequence)
    
    
        
    @staticmethod
    def create_peptide_object(single_sequence):
        peptide = peptides.Peptide(single_sequence)
        return peptide
    
    @staticmethod
    def get_aliphatic_index(peptide_object):
        return peptide_object.aliphatic_index()
    
    @staticmethod
    def get_isoelectric_point(peptide_object, scale = "EMBOSS"):
        assert scale in isoelectric_point_scales.scale, f"Please choose on of the metrics from: {isoelectric_point_scales.scale}"
        return peptide_object.isoelectric_point(pKscale = scale)
    
    @staticmethod
    def get_hydrophobicity(peptide_object, scale = "KyteDoolittle"):
        assert scale in hytdrophobicity_metrics.metrics, f"Please choose on of the metrics from: {hytdrophobicity_metrics.metrics}"
        return peptide_object.hydrophobicity(scale = scale)
    
    @staticmethod
    def get_weight(peptide_object):
        return peptide_object.molecular_weight()
    
    @staticmethod
    def mass_charge_ratio(peptide_object):
        return peptide_object.mz()
    
    def calc_attribute(self, attribute, **kwargs):
        assert attribute in self.attribute_funcs.keys(), f"Please enter a value for attribute which is in {self.attribute_funcs.keys()}"
        for sequence in self.sequences:
            PepObj = self.create_peptide_object(single_sequence=sequence)
            func = self.attribute_funcs.get(attribute)
            attribute_value = func(PepObj, **kwargs)
            self.sequence_property_interest[sequence] = attribute_value
            
    def add_attributes_to_table(self, table, attribute):
        table[attribute] = list(self.sequence_property_interest.values())
      
    



class PrepareData:
    def __init__(self):
        
        self.umap_results = None
        self.clones = None
    
    @staticmethod
    def logical_check(batch_size, pca_components, model, n_neighbors, characteristic, binding_data):
        # in the pipeline is another test, where the sample name is checked which has to be in the sequencing report 
        assert batch_size > pca_components, "The batch_size has to be larger than the pca_components"
        assert batch_size != 0, "batch_size value must not be 0" # indirectly tests that perplexity and pca components must not be 0 as well
        # more models will follow
        models_all = ["Rostlab/ProstT5_fp16", "Rostlab/prot_t5_xl_uniref50", "Rostlab/prot_t5_base_mt_uniref50", "Rostlab/prot_bert_bfd_membrane", "Rostlab/prot_t5_xxl_uniref50", "Rostlab/ProstT5", "Rostlab/prot_t5_xl_half_uniref50-enc", "Rostlab/prot_bert_bfd_ss3", "Rostlab/prot_bert_bfd_localization", "Rostlab/prot_electra_generator_bfd", "Rostlab/prot_t5_xl_bfd", "Rostlab/prot_bert", "Rostlab/prot_xlnet", "Rostlab/prot_bert_bfd", "Rostlab/prot_t5_xxl_bfd"]
        assert model in models_all, f"Please enter a valid model name which are\n{models_all}. You can find the models at: https://huggingface.co/Rostlab"
        assert n_neighbors > 1
        if characteristic != None:
            assert binding_data == None, "You cannot combine a binding data analysis and a sequence attribute analysis"
        
    @staticmethod
    def datatype_check(samples, pca_components, n_neighbors, random_seed, batch_size, model, densmap, characteristic):
        assert type(samples) == list, "You have to give a list with the samples you want to analyze"
        assert type(pca_components) == int, "You have to give an integer as input for the pca_components"
        assert type(n_neighbors) == int, "You have to give an integer as input for the perplexity"
        assert type(random_seed) == int, "You have to give an integer as input for the iterations_tsne"
        assert type(batch_size) == int, "You have to give an integer as input for the batch_size"
        assert type(model) == str, "The input for model must be a string"
        assert type(densmap) == bool, "densmap must be boolean"
        if characteristic != None:
            assert characteristic in list(GetProteinProperty().attribute_funcs.keys()), f"Please enter a valid characteristic from: {list(GetProteinProperty().attribute_funcs.keys())}"
        
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
    
    
    def return_binding_results(self,  selected_rows, antigens, region_of_interest):
        if antigens is not None:
            kds = selected_rows[antigens].max(axis = 1) # if there are multiple values for the same sequence this will find the highest one 
            self.umap_results["binding"] = list(kds)
            ids = selected_rows[antigens].idxmax(axis = 1)
        else:
            kds = None
            ids = None
        aminoacids = selected_rows[region_of_interest].to_list()
        experiments_batch = selected_rows["Experiment"]
        unique_experiments_num = pd.factorize(experiments_batch)[0]
        #self.tsne_results["experiments_factorized"] = list(unique_experiments_num)
        #self.tsne_results["experiments_string"] = list(experiments_batch)
        self.umap_results["sequences"] = list(aminoacids)
        self.umap_results['sequence_id'] = pd.Series(range(self.umap_results.shape[0]))
        return kds, ids
    
    @staticmethod
    def filter_binding_data(binding_data, region_of_interest, antigens):
        if binding_data is not None:
            merged_columns = [region_of_interest] + antigens
            binding_data = binding_data[merged_columns]
        return binding_data
    
    def label_sequence_characteristic(self, characteristic, sequences, region_of_interest):
        if characteristic != None:
            Property = GetProteinProperty(sequences)
            Property.calc_attribute(attribute=characteristic)
            Property.add_attributes_to_table(self.umap_results, attribute = characteristic, region_of_interest = region_of_interest)
        else:
            pass
            

            
    
    def tidy(self, sequencing_report, list_experiments, region_of_interest, antigens = None, batch_size = 2000, pca_components = 50,
             n_neighbors = 15, min_dist = 0.2, random_seed = 42, densmap = True, characteristic = None, model_choice = "Rostlab/prot_bert",binding_data = None,cf_column_name = "cloneFraction", sample_column_name = "Experiment"):
        """creates umap_results as class object which is a table containing all necessary data for the final plot

        Args:
            sequencing_report (pd.DataFrame): report which contains all sequencing information from all samples after the processing with mixcr
            list_experiments (list): List containing the names of the samples which should be used for the final plot. 
            region_of_interest (str): A string which indicates the column name from which the amino acid sequences should be taken from.
            antigens (list, optional): list containing the names of the antigens which are the headers of the binding data. From these antigens the binding data will be taken for the plot. Defaults to None.
            batch_size (int, optional): Equals to the number of sequences which are drawn from the sequencing report for the embedding and the final plot. Defaults to 300.
            pca_components (int, optional): Equals to the number of principal components which are used as input for tsne. This helps to remove noise and complexity before using UMAP. Defaults to 70.
            n_neighbors (int, optional): Larger values will result in more global structure being preserved at the loss of detailed local structure. In general this parameter should be between 5 to 50.
            min_dist (float, optional): controls how tightly the points will be set to each other. It should be the minimum distance points are allowed to be apart from each other in the low dimensional representation
            random_ssed (int, optional): Set a certain seed for reprodubility
            densmap (bool, optional): This parameter allows you to visualize points more densily which are also more dense in all dimensions to each other. You can have an idea about this here: https://umap-learn.readthedocs.io/en/latest/densmap_demo.html
            model_choice (str, optional): Is the final model you choose to embed your sequences. Defaults to "Rostlab/prot_bert".
            binding_data (pd.DataFrame, optional): Dataframe which contains the sequences and the binding values to the antigens. Defaults to None.
            cf_column_name (str, optional): Name of the column which contains the clone fraction in the sequencing report. Defaults to "cloneFraction".
            sample_column_name (str, optional): Name of the column which contains the sample names in the sequencing report. Defaults to "Experiment".

        Returns:
            _type_: _description_
        """
        for sample in list_experiments:
            assert sample in list(sequencing_report[sample_column_name].unique()), f"{sample} does not exist"
        if binding_data is not None:
            assert region_of_interest in binding_data.columns.to_list(), f"You must have sequences for {region_of_interest} in your binding data"    
        
        
        self.logical_check(batch_size,
                           pca_components,
                           model_choice,
                           n_neighbors,
                           characteristic,
                           binding_data
                           )
        self.datatype_check(list_experiments,
                            pca_components,
                            n_neighbors,
                            random_seed, 
                            batch_size,
                            model_choice,
                            densmap,
                            characteristic)
        binding_data = self.filter_binding_data(binding_data,
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
    #    pca_components, perplexity = self.check_warnings(sequences, pca_components, perplexity)
        X = Transformer.do_pca(sequences, batch_size, pca_components)
        peptides = selected_rows[region_of_interest].to_list()
        self.clones = selected_rows[cf_column_name]
        self.umap_results = Transformer.do_umap(X, n_neighbors, min_dist, random_seed, densmap)
        kds, ids = self.return_binding_results(selected_rows, antigens, region_of_interest)
        self.umap_results["cloneFraction"] = self.clones
        self.label_sequence_characteristic(characteristic, sequences)
        return peptides, selected_rows, kds, ids
    
    def make_csv(self):
        self.umap_results.to_csv("umap_results.csv")
        



class PlotEmbedding:
    def __init__(self,sequencing_report, model_choice, list_experiments, region_of_interest, strands,add_clone_size, batch_size, pca_components, 
                 n_neighbors, min_dist, random_seed, densmap = True, antigens = None, font_settings = {}, legend_settings = {}, ax = None, binding_data = None, # antigens mutual attribute 
                 colorbar_settings = None,  toxin_names = None, extra_figure = False):
        self.ax = ax
        self.binding_data = binding_data
        self.data_prep  = PrepareData()
        peptides, selected_rows, kds, ids = self.data_prep.tidy(sequencing_report=sequencing_report,
                                                            list_experiments=list_experiments,
                                                            region_of_interest=region_of_interest,
                                                            antigens = antigens,
                                                            batch_size = batch_size,
                                                            pca_components=pca_components,
                                                            n_neighbors = n_neighbors,
                                                            min_dist = min_dist,
                                                            random_seed = random_seed,
                                                            densmap=densmap,
                                                            model_choice=model_choice,
                                                            binding_data=binding_data)
        self.umap_results = self.data_prep.umap_results
        if self.ax != None:
            if self.binding_data is not None:
                self.create_binding_plot(kds, ids, toxin_names, colorbar_settings)
                if extra_figure == True and font_settings != {}:
                    self.create_second_bind_plot(font_settings)
                    title = "\n".join(wrap(f"UMAP embedding for {antigens}", 40))
                    self.ax.set_title(title, pad= 12, **font_settings)
            else:
                self.create_plot(selected_rows)
                if legend_settings != {}:
                    self.add_legend(list_experiments, legend_settings)
                title = "\n".join(wrap("UMAP embedding for given samples", 40))
                if font_settings != {}:
                    self.ax.set_title(title, pad= 12, **font_settings)
            if font_settings != {}:
                    self.ax.set_xlabel("UMAP_1", **font_settings) # add font_settings
            self.ax.set_ylabel("UMAP_2", **font_settings)

            if strands == True:
                self.add_seq_anotation(peptides)
            if add_clone_size != None:
                self.add_size(add_clone_size)
            
    
        
    def create_plot(self, selected_rows):

        experiments_batch = selected_rows["Experiment"]
        self.umap_plot = self.ax.scatter(self.umap_results["UMAP_1"],
                                    self.umap_results["UMAP_2"],
                                    c = pd.factorize(experiments_batch)[0],
                                    alpha = 0.5,
                                    )


            
    def add_size(self, add_clone_size):
        size_points = self.data_prep.clones * add_clone_size
        self.umap_plot.set_sizes(size_points)
        
    def add_legend(self, list_experiments, legend_settings):
        self.ax.legend(handles=self.umap_plot.legend_elements()[0],
            labels=list_experiments,
            **legend_settings)

        

    def add_seq_anotation(self,peptides):
        x = self.umap_results["UMAP_1"].values.tolist()
        y = self.umap_results["UMAP_2"].values.tolist()
        for i in range(0, len(x), 10):
            self.ax.annotate(peptides[i],
                        (x[i][0], y[i][0]),
                        fontsize = 5,
                        )
            

    
    def create_second_bind_plot(self, font_settings):
        self.fig2 = plt.figure(100)
        self.ax2 = self.fig2.gca()
        self.ax2.scatter(self.umap_results["UMAP_1"], self.umap_results["UMAP_2"], alpha = 0.0)
        self.ax2.set_xlabel('UMAP_1', **font_settings)
        self.ax2.set_ylabel('UMAP_2', **font_settings)
        n = 0
        for j, row in self.umap_results.iterrows():
            if row["binding"] > 1:
                    self.ax2.text(row['UMAP_1'], row['UMAP_1'], row['sequence_id'], fontsize=10, weight = "bold")

            else:
                if n == 6:
                    n = 0
                    self.ax2.text(row['UMAP_1'], row['UMAP_2'], row['sequence_id'], fontsize=8)
            n += 1
    
    def create_binding_plot(self, kds, ids, toxin_names, colorbar_settings):
        self.umap_plot = self.ax.scatter(self.umap_results["UMAP_1"],
                            self.umap_results["UMAP_2"],
                            c = self.umap_results.binding,
                            alpha = 1,
                            cmap = "magma"
                            )
        x_cor = list(self.umap_results.UMAP_1.iloc[:, 0])
        y_cor = list(self.umap_results.UMAP_2.iloc[:, 0])
        if toxin_names == True:
            for i, txt in enumerate(list(ids)):
                if list(kds)[i] > 0:
                    self.ax.annotate(txt, (x_cor[i], y_cor[i]))
        else:
            pass
        plt.colorbar(self.umap_plot, **colorbar_settings)
        


#sequencing_report["cloneFraction"] = sequencing_report["readFraction"]
#peptides, selected_rows, tsne_results = PrepareData().tidy(sequencing_report, ["GeneMind_1"], region_of_interest = "aaSeqCDR3", batch_size = 50)
#print(peptides)
#print(selected_rows)
#print(tsne_results)