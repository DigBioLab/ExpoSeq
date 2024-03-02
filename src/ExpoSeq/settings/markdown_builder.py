import os
import pandas as pd
import glob
import markdown
from .layout_finder import best_layout

class QuartoBuilder:
    def __init__(self, title, figure = False):
        self.content = ""
        self.title = title
        self.create_header(figure)

    def write_quarto(self, save_dir):
        save_path = os.path.join(save_dir,f"{self.title}.qmd" )
        with open(save_path, "w") as f:
            f.write(self.content)
            

    def create_header(self, figure):
        self.content += "---\n"
        self.content += f"title: {self.title}\n"
        self.content += f"format:"
        self.content += "\n     html:"
        self.content += "\n         self-contained: true"
        if figure == False:
            self.content += "\n         toc: true"
            self.content += "\n         toc-location: left"
            self.content += "\n         toc-expand: 2"
        else:
            self.content += "\n         page-layout: full"
            

        self.content += "\n"
        self.content += "---"

    def add_subplot_figures(self, files, captures):
        rows, cols = best_layout(len(files))
        self.next_row()
        inside_string = "{layout-ncol=" + str(cols) + "}"
        self.content += f"::: {inside_string}"
        self.next_row()
        if len(captures) != len(files):
            raise SystemError
        for picture , capture in zip(files, captures) :
            self.content += f"![{capture}]"
            self.content += f"({picture})"
            self.next_row()
        self.content += ":::"
        
    
        
    def end_python_code(self):
        self.next_row()
        self.content += "```"
        
    
    
    def create_headline(self, header, section = "#"):
        self.next_row()
        self.content += section + " " + header
    
    def next_row(self):
        self.content += "\n\n"
        
    def add_page(self):
        self.next_row()
        self.content += "{{< pagebreak >}}"
        self.next_row()
        
    def add_text(self, text):
        self.next_row()
        self.content += text
        

    def figure_text(self, figure_path, text,figure_description, position_text = "left"):
        self.next_row()
        self.content += ":::"
        self.content += "{layout-ncol=" + str(2) + "}"
        self.next_row()
        if position_text == "left":
            self.content += f"![]"
            self.content += f"({figure_path})"
            self.next_row()
            self.content += text
        else:
            self.content += text
            self.next_row()
            self.content += f"![]"
            self.content += f"({figure_path})"
            self.next_row()
            self.content += figure_description
            self.next_row()
            self.content += ":::"
    
                
        
    def add_figure_horizontal(self, ncol, picture_paths: list, figure_description, figure_subtitles = None):
        self.next_row()
        self.content += ":::"
        self.content += " {layout-ncol=" + str(ncol) + "}"
        self.next_row() 
        for index_pic in range(len(picture_paths)):
            if figure_subtitles != None:
                self.content += f"![{figure_subtitles[index_pic]}]"
            else:
                self.content += f"![]"
            self.content += f"({picture_paths[index_pic]})"
           # self.content += "{" + f"width={str(figure_width)}" + "%}"
            self.next_row()
        self.content += figure_description
        self.next_row()
        self.content += ":::"
        
    def add_figure_layout(self, layout_option, picture_paths: list):
        self.next_row()
        self.content += ":::"
        self.content += ' {layout=' + f'"{layout_option}"' + "}"
        self.next_row()
        for index_pic in range(len(picture_paths)):
            self.content += f"![]"
            self.content += f"({picture_paths[index_pic]})"
           # self.content += "{" + f"width={str(figure_width)}" + "%}"
            self.next_row()
        self.content += ":::"
        
    def add_figure(self, picture_path):
        self.next_row()
        self.content += ":::"
        self.content += " {layout-ncol=" + str(1) + "}"
        self.next_row()
        self.content += f"![]"
        self.content += f"({picture_path})"
        self.next_row()
        self.content += ":::"
        
    def md_to_text(self, md_file):
        self.next_row()
        with open(md_file, "r") as f:
            markdown_text = f.read()
        text = markdown.markdown(markdown_text)
        self.content += text
        
    def add_table(self, filename, top_,cols = None, max_column = None, ):
        self.next_row()
        df = pd.read_excel(filename)

        # Get the headers
        if cols != None:
            df = df[cols]
        else:
            pass
        
        if top_ != None:
            df = df.head(top_)
        else:
            df = df.nlargest(10, max_column)
            
        heads_df = df.columns.tolist()
        self.content += f"{df.to_markdown()}"
        
        

      #  headers = df.columns.tolist()

        # Convert headers to markdown format
      #  markdown_headers = ' | '.join(headers)
      #  markdown_divider = ' | '.join(['---'] * len(headers))
      #  self.content += f"{markdown_headers}\n{markdown_divider}"
        
        
        

def check_path_multiple(files:list):
    for i in files:
        if i == None:
            return False
        if os.path.isfile(i):
            pass
        else:
            return False
    return True

def find_file_with_substring(directory, substring):
    for root, dirs, files in os.walk(directory):
        for file in files:
            if substring in file:
                return os.path.join(root, file)
    return None

def find_multiple_files_with_substring(directory, substring):
    if substring == None:
        return None
    matching_files = []
    all_files = glob.glob(directory + "/*")
    for file in all_files:
        if substring in file:
            matching_files.append(os.path.join(directory, file))
    return matching_files
    

    


def create_quarto(experiment, plot_path, binding_data, samples):
    assert os.path.isdir(plot_path), f"The directory {plot_path} does not exist"
    assert os.path.isdir(os.path.join(plot_path, "sequence_cluster")), f"The directory {os.path.join(plot_path, 'sequence_cluster')} does not exist"
    assert os.path.isdir(os.path.join(plot_path, "sequence_embedding")), f"The directory {os.path.join(plot_path, 'sequence_embedding')} does not exist"
    assert os.path.isdir(os.path.join(plot_path, "sequence_cluster", "reports")), f"The directory {os.path.join(plot_path, 'sequence_cluster', 'reports')} does not exist"
    assert os.path.isdir(os.path.join(plot_path, "length_distributions")), f"The directory {os.path.join(plot_path, 'length_distributions')} does not exist"
    assert os.path.isdir(os.path.join(plot_path, "rarefraction_curves")), f"The directory {os.path.join(plot_path, 'rarefraction_curves')} does not exist"
    assert os.path.isdir(os.path.join(plot_path, "logo_plots")), f"The directory {os.path.join(plot_path, 'logo_plots')} does not exist"
    Builder = QuartoBuilder(experiment)
    region = os.path.basename(plot_path)
    Builder.create_headline(f"Basic Overview - {region}")
    Builder.create_headline("Theory", section = "##")
    if os.path.isfile(os.path.join("settings", "rarefraction_curves_desc.md")):
        Builder.md_to_text(os.path.join("settings", "rarefraction_curves_desc.md"))
    if os.path.isfile(os.path.join("settings", "alignment_plot.md")):
        Builder.md_to_text(os.path.join("settings", "alignment_plot.md"))
    if os.path.isfile(os.path.join("settings", "morosita_horn.md")):
        Builder.md_to_text(os.path.join("settings", "morosita_horn.md"))
    if os.path.isfile(os.path.join("settings", "jaccard.md")):
        Builder.md_to_text(os.path.join("settings", "logo_plot.md"))

    #Builder.
    Builder.create_headline("Identity between samples", section = "##")
    Builder.create_headline("Identity based on Morosita Horn Index", section = "###")
    morosita_horn_path = find_file_with_substring(plot_path, "morosita_horn")
    if morosita_horn_path != None:
        Builder.add_figure(morosita_horn_path)
    else: pass
    Builder.create_headline("Identity based on Jaccard Index", section = "###")
    jaccard_path = find_file_with_substring(plot_path, "jaccard")
    if jaccard_path != None:
        Builder.add_figure(jaccard_path)
    else: pass
    Builder.add_page()
    Builder.create_headline("Identity based on Sorensen Index", section = "###")
    sorensen_path = find_file_with_substring(plot_path, "sorensen")
    Builder.add_page()
    if sorensen_path != None:
        Builder.add_figure(sorensen_path)
    else: pass
    Builder.add_page()
    
    Builder.create_headline("Sequencing Quality of samples", section = "##")
    figure_alignment = os.path.join(plot_path, "alignment_quality.png")
    figure_rarefraction = os.path.join(plot_path, "rarefraction_all.png")
    
    if check_path_multiple([figure_alignment, figure_rarefraction]) != False:
        Builder.add_figure_horizontal(2, [figure_alignment, figure_rarefraction], "")
    else:
        for i in [figure_alignment, figure_rarefraction]:
            if i != None:
                Builder.add_figure(i)
            else:
                pass
    Builder.add_page()
    # Diversity Plots
    figure_diversity_shannon = os.path.join(plot_path, "diversity_Shannon.png")
    figure_diversity_simp = os.path.join(plot_path, "diversity_InverseSimpson.png")
    if check_path_multiple([figure_diversity_shannon, figure_diversity_simp]) != False:
        Builder.add_figure_horizontal(2, [figure_diversity_shannon, figure_diversity_simp], "")
    else:
        for i in [figure_diversity_shannon, figure_diversity_simp]:
            if i != None:
                Builder.add_figure(i)
            else:
                pass
    Builder.add_page()
    
    Builder.create_headline("Sample specific quality analysis", section = "##")
    length_distribution_path = os.path.join(plot_path, "length_distributions")
    rarefraction_curve_path = os.path.join(plot_path, "rarefraction_curves")
    clone_fraction_path = os.path.join(plot_path, "clone_fraction")
    
    for sample in samples:
        Builder.create_headline("Quality analysis for " + sample, section = "###")
        length_distr_single = find_file_with_substring(length_distribution_path, sample)
        rarefraction_single = find_file_with_substring(rarefraction_curve_path, sample)
        clone_single = find_file_with_substring(clone_fraction_path, sample)
        if check_path_multiple([length_distr_single, rarefraction_single, clone_single]) != False:
            Builder.add_figure_layout([[1, 1], [1]], [length_distr_single, rarefraction_single, clone_single])
        else:
            for i in [length_distr_single, rarefraction_single, clone_single]:
                if i != None:
                    Builder.add_figure(i)
                else: 
                    pass
        Builder.add_page()
        
    Builder.create_headline("Sequence Logo Plots", section = "##")
    logo_plot_path = os.path.join(plot_path, "logo_plots")
    for sample in samples:
        Builder.create_headline("Logo Plot for sequences with length 16 for " + sample, section = "###")
        logo_plot_single = find_file_with_substring(logo_plot_path, sample)
        if logo_plot_single != None:
            Builder.add_figure(logo_plot_single)
        else:
            pass
        Builder.add_page()

    ls_clusters_dir = os.path.join(plot_path, "sequence_cluster")
    Builder.create_headline("Sequence Clusters")
    Builder.create_headline("Theory", section = "##")
    if os.path.isfile(os.path.join("settings", "levenshtein_cluster.md")):
        Builder.md_to_text(os.path.join("settings", "levenshtein_cluster.md"))
    if os.path.isfile(os.path.join("settings", "dendrogram_cluster.md")):
        Builder.md_to_text(os.path.join("settings", "dendrogram_cluster.md"))
    if os.path.isfile(os.path.join("settings", "sequence_embedding.md")):
        Builder.md_to_text(os.path.join("settings", "sequence_embedding.md"))

    report_seq_cluster = os.path.join(ls_clusters_dir, "reports")
    
    Builder.create_headline("Levenshtein distance based", section="##")
    Builder.create_headline("Relation of all samples", section = "###")
    _isfile = find_file_with_substring(plot_path, "ls_connection_all")
    if _isfile != None:
        Builder.add_figure(_isfile)
    else:
        pass
    for sample in samples:
        
        plot_files = find_multiple_files_with_substring(ls_clusters_dir, sample)
        if check_path_multiple(plot_files):
            Builder.create_headline(f"Sequence Similarity based on Levenshtein Distance of {sample}", "###")
            Builder.add_figure_horizontal(ncol = 2, picture_paths=plot_files, figure_description="")
        else:
            Builder.create_headline(f"Sequence Similarity based on Levenshtein Distance of {sample}", "###")
            for i in plot_files:
                if i != None:
                    Builder.add_figure(i)
                else:
                    pass
        try:
            report_file = find_file_with_substring(report_seq_cluster, sample)
            if report_file != None:
                Builder.add_table(report_file, top_ = 10)
        except:
            pass

        Builder.add_page()

    Builder.create_headline("Sequence Similarity between samples", section = "##")
    embedding_dir = os.path.join(plot_path, "sequence_embedding", "sgt")
    Builder.create_headline("SGT Embedder", section = "###")
    for sample in samples:
        file = find_file_with_substring(embedding_dir, sample)
        if file != None:
            Builder.create_headline(f"Global Sequence similarity for {sample}", section = "####")
            Builder.add_figure(file)
            Builder.add_page()
        else:
            pass
        
    Builder.create_headline("ProtBert Embedder", section = "###")
    embedding_dir = os.path.join(plot_path, "sequence_embedding", "protbert")
    for sample in samples:
        file = find_file_with_substring(embedding_dir, sample)
        if file != None:
            Builder.create_headline(f"Global Sequence similarity for {sample}", section = "####")
            Builder.add_figure(file)
            Builder.add_page()
        else:
            pass
    
    Builder.create_headline("T5 Embedder", section = "###")
    embedding_dir = os.path.join(plot_path, "sequence_embedding", "T5")
    for sample in samples:
        file = find_file_with_substring(embedding_dir, sample)
        if file != None:
            Builder.create_headline(f"Global Sequence similarity for {sample}", section = "####")
            Builder.add_figure(file)
            Builder.add_page()
        else:
            pass
    
    Builder.write_quarto(save_dir=plot_path)
    
    if binding_data is not None:
        if os.path.basename(plot_path) in binding_data.columns.tolist():
            assert os.path.isdir(os.path.join(plot_path, "clustering_antigens")), f"The directory {os.path.join(plot_path, 'clustering_antigens')} does not exist"
            assert os.path.isdir(os.path.join(plot_path, "clustering_antigens", "ls_binding_cluster")), f"The directory {os.path.join(plot_path, 'clustering_antigens', 'ls_binding_cluster')} does not exist"
            assert os.path.isdir(os.path.join(plot_path, "clustering_antigens", "dendro_binding")), f"The directory {os.path.join(plot_path, 'clustering_antigens', 'dendro_binding')} does not exist"
            assert os.path.isdir(os.path.join(plot_path, "clustering_antigens", "reports")), f"The directory {os.path.join(plot_path, 'clustering_antigens', 'reports')} does not exist"
            Builder.create_headline("Cluster Antigens")
            Builder.create_headline("Theory", section = "##")
            if os.path.isfile(os.path.join("settings", "reports_binding.md")):
                Builder.md_to_text(os.path.join("settings", "reports_binding.md"))

            embedding_antigen_path = os.path.join(plot_path, "clustering_antigens")
            embedding_antigen_report = os.path.join(embedding_antigen_path, "reports")
            Builder.create_headline("Cluster with embedding", section = "##")
            for sample in samples:
                file = find_file_with_substring(embedding_antigen_path, sample)
                if file != None:
                    Builder.create_headline(f"Sequence embedding of {sample}", section = "###")
                    Builder.add_figure(file)
                    try:
                        report = find_file_with_substring(embedding_antigen_report, sample)
                        cols_choice = ["tsne1", "tsne2", "binding", "sequences", "sequence_id"]
                        cur_df = pd.read_excel(report)
                        df_cols = cur_df.columns.tolist()
                        incorrect_cols = [x for x in cols_choice if x not in df_cols]
                        assert not incorrect_cols, "Either the user changed the default column header or the developer for the tsne clusters with the antigen"
                        Builder.add_table(report, top_ = None, cols = cols_choice, max_column = "binding")
                    except:
                        print("Something went wrong with creating the tsne binding plots in the report")
                else:
                    pass
            ls_cluster_binding = os.path.join(plot_path,"clustering_antigens", "ls_binding_cluster")
            
            ls_cluster_binding_report = os.path.join(ls_cluster_binding, "reports")
            dendro_cluster_binding = os.path.join(plot_path,"clustering_antigens", "dendro_binding")
            Builder.create_headline("Cluster binding data based on levenshtein distance", section = "##")
            for sample in samples:
                ls_file = find_file_with_substring(ls_cluster_binding, sample)
                dendro_file = find_file_with_substring(dendro_cluster_binding, sample)
                if ls_file != None:
                    Builder.create_headline(f"LS-Distance cluster for {sample}", section = "###")
                    if dendro_file != None:
                        Builder.add_figure(ls_file)
                        Builder.add_figure(dendro_file)
                #      Builder.add_figure_horizontal(ncol = 2, picture_paths=[ls_file, dendro_file], figure_description="Cluster visualizations based on LS-distance")
                    else:
                        Builder.add_figure(ls_file)
                    try:
                        report_file = find_file_with_substring(ls_cluster_binding_report, sample)
                        Builder.add_table(report_file, top_ = 10)
                    except:
                        print("Could not create report table for ls cluster")
                else:
                    pass
        else:
            print(f"No plots with binding data were included in the report.\n If you have plots with binding data in {plot_path} please make sure that you upload the binding data to the pipeline with plot.add_binding_data() before you run plot.create_report()")
    Builder.write_quarto(save_dir=plot_path)

