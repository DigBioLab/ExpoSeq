from textwrap import wrap
from .tidy_protbert_embedding import TransformerBased
import pandas as pd
import matplotlib.pyplot as plt
import warnings
from .contents.simple_protein_property import GetProteinProperty
import numpy as np
import os
import pickle


class PrepareData:
    def __init__(self):

        self.umap_results = pd.DataFrame([])
        self.clones = None

    @staticmethod
    def logical_check(
        batch_size, pca_components, model, n_neighbors, characteristic, binding_data
    ):
        # in the pipeline is another test, where the sample name is checked which has to be in the sequencing report
        assert (
            batch_size > pca_components
        ), "The batch_size has to be larger than the pca_components"
        assert (
            batch_size != 0
        ), "batch_size value must not be 0"  # indirectly tests that perplexity and pca components must not be 0 as well
        # more models will follow
        models_all = [
            "Rostlab/ProstT5_fp16",
            "Rostlab/prot_t5_xl_uniref50",
            "Rostlab/prot_t5_base_mt_uniref50",
            "Rostlab/prot_bert_bfd_membrane",
            "Rostlab/prot_t5_xxl_uniref50",
            "Rostlab/ProstT5",
            "Rostlab/prot_t5_xl_half_uniref50-enc",
            "Rostlab/prot_bert_bfd_ss3",
            "Rostlab/prot_bert_bfd_localization",
            "Rostlab/prot_electra_generator_bfd",
            "Rostlab/prot_t5_xl_bfd",
            "Rostlab/prot_bert",
            "Rostlab/prot_xlnet",
            "Rostlab/prot_bert_bfd",
            "Rostlab/prot_t5_xxl_bfd",
        ]
        assert (
            model in models_all
        ), f"Please enter a valid model name which are\n{models_all}. You can find the models at: https://huggingface.co/Rostlab"
        assert n_neighbors > 1, "The number of neighbors must be larger than 1"
        assert (
            n_neighbors < batch_size
        ), "The number of neighbors must be smaller than the batch size"
        if characteristic == "binding":
            assert (
                binding_data is not None
            ), f"You must provide binding data for the characteristic option binding"

        if characteristic != None and characteristic != "binding":
            assert (
                binding_data == None
            ), "You cannot combine a binding data analysis and a sequence attribute analysis"

    @staticmethod
    def datatype_check(
        samples,
        pca_components,
        n_neighbors,
        random_seed,
        batch_size,
        model,
        densmap,
        characteristic,
        add_clone_size,
        metric,
        min_dist,
    ):
        assert (
            type(samples) == list
        ), "You have to give a list with the samples you want to analyze"
        assert (
            type(pca_components) == int
        ), "You have to give an integer as input for the pca_components"
        assert (
            type(n_neighbors) == int
        ), "You have to give an integer as input for the perplexity"
        assert (
            type(random_seed) == int
        ), "You have to give an integer as input for the iterations_tsne"
        assert (
            type(batch_size) == int
        ), "You have to give an integer as input for the batch_size"
        assert type(model) == str, "The input for model must be a string"
        assert type(densmap) == bool, "densmap must be boolean"
        if characteristic != None:
            assert characteristic in list(
                GetProteinProperty([]).attribute_funcs.keys()
            ), f"Please enter a valid characteristic from: {list(GetProteinProperty([]).attribute_funcs.keys())}"
        if add_clone_size != None:
            assert (
                type(add_clone_size) == int
            ), "add_clone_size must be an integer if it is not None"
        else:
            assert add_clone_size == None
        avail_metrics = [
            "euclidean",
            "manhatten",
            "chebyshev",
            "minkowski",
            "canberra",
            "braycurtis",
            "haversine",
            "mahalanobis",
            "wminkowski",
            "seuclidean",
            "cosine",
            "correlation",
        ]
        assert (
            metric in avail_metrics
        ), f"Please choose one of the metrics from this list: {avail_metrics}"
        possible_characteristics = [
            "isoelecrtric_point",
            "aliphatic_index",
            "hydrophobicity",
            "weight",
            "mass_charge_ratio",
            "length",
            None,
            "binding",
        ]
        assert (
            characteristic in possible_characteristics
        ), f"Please choose one of the values from: {possible_characteristics}"
        assert type(min_dist) == float, "min_dist must be a float"
        assert min_dist > 0, "min_dist must be larger than 0"

    @staticmethod
    def check_warnings(sequences, pca_components, perplexity):
        if len(sequences) < pca_components:
            warnings.warn(
                f"The number of sequences you have is {len(sequences)} but you need to have more sequences than principal components which is: {pca_components}"
            )
            print(
                f"Number of principal components is set to number of sequences ({len(sequences)})"
            )
            pca_components = len(sequences)
        if pca_components < perplexity:
            warnings.warn(
                "The number of reduced dimensions you have is "
                + str(pca_components)
                + "but you need to have more than perplexity which is: "
                + str(perplexity)
            )
            print(
                f"Perplexity is set to the half of reduced dimensions ({pca_components//2})"
            )
            perplexity = pca_components // 2
            if perplexity < 1:
                perplexity = 1
        return pca_components, perplexity

    @staticmethod
    def return_binding_results(
        selected_rows, antigens, region_of_interest, add_clone_size, umap_results
    ):
        """Creates the tsne results objects based on the information in selected rows.

        Args:
            selected_rows (_type_): _description_
            antigens (_type_): _description_
            region_of_interest (_type_): _description_
            add_clone_size (_type_): _description_
        """
        if antigens is not None:
            kds = selected_rows[antigens].max(
                axis=1
            )  # if there are multiple values for the same sequence this will find the highest one

            umap_results["binding"] = list(kds)
            ids = selected_rows[antigens].idxmax(axis=1)
            umap_results["highest_binder"] = list(ids)

        else:
            kds = None
            ids = None
        aminoacids = selected_rows[region_of_interest].to_list()
        experiments_batch = selected_rows["Experiment"]
        experiments_batch = experiments_batch.replace(
            0, "non-merged binding data"
        )  # this happens because of the merge in TransformerBased and you need to change the label
        unique_experiments_num = list(pd.factorize(experiments_batch)[0])
        umap_results["experiments_string"] = experiments_batch.to_list()
        umap_results["experiments_factorized"] = unique_experiments_num
        umap_results["sequences"] = list(aminoacids)
        umap_results["sequence_id"] = list(range(umap_results.shape[0]))
        if add_clone_size != None:
            umap_results["size"] = (
                np.array(selected_rows["cloneFraction"].to_list()) * add_clone_size
            )
        else:
            umap_results["size"] = (
                len(selected_rows["cloneFraction"].to_list()) * [add_clone_size]
            )
        # self.umap_results.reset_index(inplace = True, drop = True)
        return kds, ids, umap_results

    @staticmethod
    def filter_binding_data(binding_data, region_of_interest, antigens):
        if binding_data is not None:
            merged_columns = [region_of_interest] + antigens
            binding_data = binding_data[merged_columns]
        return binding_data

    @staticmethod
    def label_sequence_characteristic(
        characteristic, sequences, selected_rows, antigens, umap_results
    ):
        if characteristic != None:
            if characteristic != "binding":
                Property = GetProteinProperty(sequences)
                Property.calc_attribute(attribute=characteristic)
                umap_results = Property.add_attributes_to_table(
                    umap_results, attribute=characteristic
                )
                property_result = list(Property.sequence_property_interest.values())
            else:
                kds = selected_rows[antigens].max(axis=1)
                property_result = kds.fillna(0)
        else:
            property_result = None
            pass
        return property_result, umap_results

    @staticmethod
    def save_current_embedding(
        X_path: str,
        sequences_list: np.array,
        batch_size: int,
        antigens,
        region_of_interest: str,
        model_choice: str,
        binding_data: pd.DataFrame,
        characteristic: str,
    ):
        """This method saves the sequences_list from which contains the array with the high dimensions as npz file.
        Further, it creates in the same directory a pkl file which contains the parameters for that specific embedding.


        Args:
            X_path (str): Path where yhe array should be saved
            sequences_list (np.array): array with the data from the embedding.
            batch_size (int): Number of sequences per sample
            antigens (list): Value for antigens which should be incorporated. either none or list of strings
            region_of_interest (str): Region of interest
            model_choice (str): Model which was chosen for the embedding
            binding_data (pd.DataFrame): pandas dataframe with the embedding
        """
        np.savez_compressed(X_path, a=sequences_list)
        current_run = {
            "batch_size": batch_size,
            "antigens": antigens,
            "region_of_interest": region_of_interest,
            "model_choice": model_choice,
            "binding_data": binding_data,
            "characteristic": characteristic,
        }
        dir_X = os.path.dirname(X_path)
        np.save(os.path.join(dir_X, "current_run.npy"), current_run)

    @staticmethod
    def load_and_compare_past_run(
        X_path,
        batch_size: int,
        antigens,
        region_of_interest: str,
        model_choice: str,
        binding_data: pd.DataFrame,
        characteristic: str,
    ):
        """Loads the np array where the embedding is saved and compares the settings of the past run to the settings of the current run and if any is different res will be False and the embedding will be repeated.

        Args:
            X_path (_type_): _description_
            batch_size (int): _description_
            antigens (_type_): _description_
            region_of_interest (str): _description_
            model_choice (str): _description_
            binding_data (pd.DataFrame): _description_

        Returns:
            _type_: _description_
        """
        saved_array = np.load(X_path)
        sequences_list = saved_array["a"]
        dir_X = os.path.dirname(X_path)
        try:
            old_settings = np.load(
                os.path.join(dir_X, "current_run.npy"), allow_pickle=True
            ).item()

        except:
            res = False
            return res, sequences_list
        new_settings = {
            "batch_size": batch_size,
            "antigens": antigens,
            "region_of_interest": region_of_interest,
            "model_choice": model_choice,
            "binding_data": binding_data,
            "characteristic": characteristic,
        }
        binding_data_old = old_settings.get("binding_data")
        binding_data_new = new_settings.get("binding_data")
        new_settings.pop("binding_data")
        old_settings.pop("binding_data")
        res = all(
            (new_settings.get(k) == v for k, v in old_settings.items())
        ) and binding_data_old.equals(binding_data_new)

        return res, sequences_list

    def tidy(
        self,
        sequencing_report,
        list_experiments,
        region_of_interest,
        antigens=None,
        batch_size=2000,
        pca_components=50,
        n_neighbors=50,
        min_dist=0.2,
        random_seed=42,
        densmap=True,
        metric="cosine",
        characteristic=None,
        eps_dbscan=0.5,
        min_pts_dbscan=2,
        add_clone_size=500,
        model_choice="Rostlab/prot_bert",
        binding_data=None,
        cf_column_name="cloneFraction",
        sample_column_name="Experiment",
        X_path="temp/current_array.npz",
        number_jobs=-1,
    ):
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
            random_seed (int, optional): Set a certain seed for reprodubility
            densmap (bool, optional): This parameter allows you to visualize points more densily which are also more dense in all dimensions to each other. You can have an idea about this here: https://umap-learn.readthedocs.io/en/latest/densmap_demo.html
            metric (str, optional): You need to insert a string as input which is the distance metric for the UMAP algorithm.
            eps (float, optional): Parameter for DBSCAN. Maximum distance between two points to still form one cluster. Defaults to 3.
            min_pts (int, optional): Parameter for DBSCAN. Fewest number of points required to form a cluster. Defaults to 4.
            model_choice (str, optional): Is the final model you choose to embed your sequences. Defaults to "Rostlab/prot_bert".
            binding_data (pd.DataFrame, optional): Dataframe which contains the sequences and the binding values to the antigens. Defaults to None.
            cf_column_name (str, optional): Name of the column which contains the clone fraction in the sequencing report. Defaults to "cloneFraction".
            sample_column_name (str, optional): Name of the column which contains the sample names in the sequencing report. Defaults to "Experiment".

        Returns:
            _type_: _description_
        """
        batch_size = batch_size * len(list_experiments)
        for sample in list_experiments:
            assert sample in list(
                sequencing_report[sample_column_name].unique()
            ), f"{sample} does not exist"
        if binding_data is not None:
            assert (
                region_of_interest in binding_data.columns.to_list()
            ), f"You must have sequences for {region_of_interest} in your binding data"

        self.logical_check(
            batch_size,
            pca_components,
            model_choice,
            n_neighbors,
            characteristic,
            binding_data,
        )
        self.datatype_check(
            list_experiments,
            pca_components,
            n_neighbors,
            random_seed,
            batch_size,
            model_choice,
            densmap,
            characteristic,
            add_clone_size,
            metric,
            min_dist,
        )
        binding_data = self.filter_binding_data(
            binding_data, region_of_interest, antigens
        )

        sequences, sequences_filtered, selected_rows = (
            TransformerBased.filter_sequences(
                sequencing_report,
                batch_size / len(list_experiments),
                list_experiments,
                binding_data,
                region_of_interest=region_of_interest,
                cf_column_name=cf_column_name,
                sample_column_name=sample_column_name,
            )
        )
        peptides = selected_rows[region_of_interest].to_list()
        self.clones = selected_rows[cf_column_name].to_list()

        if os.path.exists(X_path):
            res, sequences_list = self.load_and_compare_past_run(
                X_path,
                batch_size,
                antigens,
                region_of_interest,
                model_choice,
                binding_data,
                characteristic,
            )
        else:
            res = False

        if res == False:
            Transformer = TransformerBased(choice=model_choice)
            sequences_list = Transformer.embedding_per_seq(sequences)
            self.save_current_embedding(
                X_path,
                sequences_list,
                batch_size,
                antigens,
                region_of_interest,
                model_choice,
                binding_data,
                characteristic,
            )

        X = TransformerBased.do_pca(sequences_list, pca_components)
        property_result, umap_results = self.label_sequence_characteristic(
            characteristic, sequences, selected_rows, antigens, self.umap_results
        )

        umap_results, reduced_dim = TransformerBased.do_umap(
            X,
            n_neighbors,
            min_dist,
            random_seed,
            densmap,
            y=property_result,
            metric=metric,
            number_jobs=number_jobs,
        )
        get_clusters = TransformerBased.cluster_with_hdbscan(
            umap_results, eps=eps_dbscan, min_pts=min_pts_dbscan
        ).tolist()
        if property_result != None:
            umap_results[characteristic] = property_result
        umap_results["cluster_id"] = get_clusters
        umap_results["cloneFraction"] = self.clones
        kds, ids, umap_results = self.return_binding_results(
            selected_rows, antigens, region_of_interest, add_clone_size, umap_results
        )  # add columns to report
        for sample in list_experiments:  # tests
            assert (
                umap_results[umap_results["experiments_string"] == sample].shape[0] >= 1
            ), f"After processing your data for your parameters no sequences are left for {sample}"
        assert type(umap_results["experiments_string"].tolist()) == list
        self.umap_results = umap_results
        return peptides, selected_rows, kds, ids

    def make_csv(self):
        self.umap_results.to_csv("umap_results.csv")


class PlotEmbedding:
    def __init__(
        self,
        sequencing_report,
        list_experiments,
        region_of_interest,
        strands=True,
        add_clone_size=300,
        batch_size=500,
        pca_components=50,
        n_neighbors=45,
        min_dist=0.01,
        random_seed=42,
        densmap=True,
        metric="cosine",
        model_choice="Rostlab/prot_bert",
        eps_dbscan=0.5,
        min_pts_dbscan=2,
        characteristic=None,
        antigens=None,
        font_settings={},
        legend_settings={},
        ax=None,
        binding_data=None,  # antigens mutual attribute
        colorbar_settings=None,
        toxin_names=None,
        extra_figure=False,
        prefered_cmap="inferno",
        number_jobs=-1,
    ):
        self.ax = ax
        self.binding_data = binding_data
        self.data_prep = PrepareData()
        peptides, selected_rows, kds, ids = self.data_prep.tidy(
            sequencing_report=sequencing_report,
            list_experiments=list_experiments,
            region_of_interest=region_of_interest,
            antigens=antigens,
            batch_size=batch_size,
            pca_components=pca_components,
            n_neighbors=n_neighbors,
            metric=metric,
            min_dist=min_dist,
            random_seed=random_seed,
            densmap=densmap,
            characteristic=characteristic,
            add_clone_size=add_clone_size,
            model_choice=model_choice,
            binding_data=binding_data,
            number_jobs=number_jobs,
        )
        self.umap_results = self.data_prep.umap_results

        if self.ax != None:
            if self.binding_data is not None:
                self.create_binding_plot(
                    kds, ids, toxin_names, colorbar_settings, prefered_cmap
                )
                if extra_figure == True and font_settings != {}:
                    self.create_second_bind_plot(font_settings)
                    title = "\n".join(wrap(f"UMAP embedding for {antigens}", 40))
                    self.ax.set_title(title, pad=12, **font_settings)
            else:
                sm = self.create_plot(characteristic, prefered_cmap)
                title = "\n".join(wrap("UMAP embedding for given samples", 40))
                if font_settings != {}:
                    self.ax.set_title(title, pad=12, **font_settings)
                if colorbar_settings != {} and characteristic != None:
                    self.add_colorbar(colorbar_settings, characteristic, sm)

            if font_settings != {}:
                self.ax.set_xlabel("UMAP_1", **font_settings)  # add font_settings
                self.ax.set_ylabel("UMAP_2", **font_settings)

            if strands == True:
                self.add_seq_anotation(peptides)

            if legend_settings != {}:
                self.add_legend(legend_settings)

    def create_plot(self, characteristic, prefered_cmap):
        markers = ["o", "+", "x", "s", "p", "x", "D"]
        if characteristic == None:
            self.umap_results["color"] = self.umap_results["cluster_id"]
            sm = None
        else:
            self.umap_results["color"] = self.umap_results[characteristic]
            global_min_color = self.umap_results["color"].min()
            global_max_color = self.umap_results["color"].max()
            sm = plt.cm.ScalarMappable(
                cmap=prefered_cmap,
                norm=plt.Normalize(vmin=global_min_color, vmax=global_max_color),
            )

            norm = plt.Normalize(vmin=global_min_color, vmax=global_max_color)

        unique_experiments = self.umap_results["experiments_string"].unique()

        for index, experiment in enumerate(unique_experiments):
            local_results = self.umap_results[
                self.umap_results["experiments_string"] == experiment
            ]
            umap_1_values = local_results["UMAP_1"]
            umap_2_values = local_results["UMAP_2"]
            if characteristic is not None:
                self.ax.scatter(
                    umap_1_values,
                    umap_2_values,
                    marker=markers[index],
                    s=local_results["size"],
                    c=local_results["color"],
                    alpha=0.5,
                    norm=norm,
                    cmap=prefered_cmap,
                    label=experiment,
                )

            else:
                self.ax.scatter(
                    umap_1_values,
                    umap_2_values,
                    marker=markers[index],
                    c=local_results["color"],
                    s=local_results["size"],
                    alpha=0.5,
                    cmap="tab20c",
                    label=experiment,
                )

        return sm

    def add_colorbar(self, colorbar_settings, label, sm):
        colorbar_settings["orientation"] = "horizontal"
        #    del colorbar_settings["spacing"]
        fig = self.ax.get_figure()
        fig.colorbar(
            sm,
            ax=self.ax,
            label=label,
            **colorbar_settings,
        )

    def add_legend(self, legend_settings):
        self.ax.legend(**legend_settings)

    def add_seq_anotation(self, peptides):
        x = self.umap_results["UMAP_1"].values.tolist()
        y = self.umap_results["UMAP_2"].values.tolist()
        for i in range(0, len(x), 10):
            self.ax.annotate(
                peptides[i],
                (x[i], y[i]),
                fontsize=5,
            )

    def create_second_bind_plot(self, font_settings):
        self.fig2 = plt.figure(100)
        self.ax2 = self.fig2.gca()
        self.ax2.scatter(
            self.umap_results["UMAP_1"], self.umap_results["UMAP_2"], alpha=0.0
        )
        self.ax2.set_xlabel("UMAP_1", **font_settings)
        self.ax2.set_ylabel("UMAP_2", **font_settings)
        n = 0
        for j, row in self.umap_results.iterrows():
            if row["binding"] > 1:
                self.ax2.text(
                    row["UMAP_1"],
                    row["UMAP_1"],
                    row["sequence_id"],
                    fontsize=10,
                    weight="bold",
                )

            else:
                if n == 6:
                    n = 0
                    self.ax2.text(
                        row["UMAP_1"], row["UMAP_2"], row["sequence_id"], fontsize=8
                    )
            n += 1

    def get_sm(self, column_char, prefered_cmap):
        self.umap_results["color"] = column_char
        global_min_color = self.umap_results["color"].min()
        global_max_color = self.umap_results["color"].max()
        sm = plt.cm.ScalarMappable(
            cmap=prefered_cmap,
            norm=plt.Normalize(vmin=global_min_color, vmax=global_max_color),
        )

        norm = plt.Normalize(vmin=global_min_color, vmax=global_max_color)
        return sm, norm

    def create_binding_plot(
        self, kds, ids, toxin_names, colorbar_settings, prefered_cmap="magma"
    ):
        markers = ["o", "+", "x", "s", "p", "x", "D"]
        self.umap_results["color"] = self.umap_results["binding"]
        sm, norm = self.get_sm(self.umap_results["binding"], prefered_cmap)

        unique_experiments = self.umap_results["experiments_string"].unique()

        for index, experiment in enumerate(unique_experiments):
            local_results = self.umap_results[
                self.umap_results["experiments_string"] == experiment
            ]
            umap_1_values = local_results["UMAP_1"]
            umap_2_values = local_results["UMAP_2"]

            self.ax.scatter(
                umap_1_values,
                umap_2_values,
                c=local_results["binding"],
                marker=markers[index],
                s=local_results["size"],
                alpha=0.9,
                cmap=prefered_cmap,
                label=experiment,
            )

            x_cor = local_results["UMAP_1"].tolist()
            y_cor = local_results["UMAP_2"].tolist()
            if toxin_names == True:
                for i, txt in enumerate(list(local_results["highest_binder"])):
                    if list(local_results["kds"])[i] > 0:
                        self.ax.annotate(txt, (x_cor[i], y_cor[i]))
            else:
                pass

        self.add_colorbar(colorbar_settings, "Binding", sm)


# PrepData = PrepareData()
# peptides, selected_rows, kds, ids = PrepData.tidy(sequencing_report, list_experiments, "aaSeqCDR3", batch_size = 80, characteristic = "length")
