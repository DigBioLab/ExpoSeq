import os
import pandas as pd
from tabulate import tabulate
import glob
import markdown


class QuartoBuilder:
    def __init__(self, title):
        self.content = ""
        self.title = title
        self.create_header()

    def write_quarto(self, save_dir):
        save_path = os.path.join(save_dir,f"{self.title}.qmd" )
        with open(save_path, "w") as f:
            f.write(self.content)
            
    def create_header(self):
        self.content += "---\n"
        self.content += f"title: {self.title}\n"
        self.content += f"format:"
        self.content += "\n     html:"
        self.content += "\n         toc: true"
        self.content += "\n         toc-location: left"
        self.content += "\n         toc-expand: 2"
        self.content += "\n"
        self.content += "---"

    
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
        figure_width = 100/len(picture_paths) 
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

def find_multiple_files_with_substring(directory, substring):
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
    Builder.create_headline("Basic Overview")
    
    Builder.md_to_text(os.path.join("settings", "rarefraction_curves_desc.md"))
    Builder.md_to_text(os.path.join("settings", "alignment_plot.md"))
    Builder.md_to_text(os.path.join("settings", "morosita_horn.md"))
    Builder.md_to_text(os.path.join("settings", "logo_plot.md"))
    #Builder.
    Builder.create_headline("Identity between samples", section = "##")
    Builder.add_figure_horizontal(2, [r"C:\Users\nilsh\my_projects\ExpoSeq\my_experiments\test_all_samples\plots\aaSeqCDR3\jaccard.png", r"C:\Users\nilsh\my_projects\ExpoSeq\my_experiments\test_all_samples\plots\aaSeqCDR3\morosita_horn.png"], "test description")
    Builder.add_page()
    Builder.create_headline("Sequencing Quality of samples", section = "##")
    figure_alignment = os.path.join(plot_path, "alignment_quality.png")
    figure_rarefraction = os.path.join(plot_path, "rarefraction_all.png")
    text_alignment = "- *If some bars are very low, it might indicate that something went wrong during the wet lab process or the sequencing.*\n- *If there is a big proportion of the orange bar, you might have chosen the wrong method for mixcr.*"
    Builder.add_figure_horizontal(2, [figure_alignment, figure_rarefraction], "Aligment Quality")
    Builder.add_text(text_alignment)

    Builder.create_headline("Sample specific quality analysis", section = "##")
    length_distribution_path = os.path.join(plot_path, "length_distributions")
    rarefraction_curve_path = os.path.join(plot_path, "rarefraction_curves")
    logo_plot_path = os.path.join(plot_path, "logo_plots")
    for sample in samples:
        Builder.create_headline("Quality analysis for " + sample, section = "###")
        length_distr_single = find_file_with_substring(length_distribution_path, sample)
        rarefraction_single = find_file_with_substring(rarefraction_curve_path, sample)
        logo_plot_single = find_file_with_substring(logo_plot_path, sample)
        if check_path_multiple([length_distr_single, rarefraction_single, logo_plot_single]) != False:
            Builder.add_figure_layout([[1, 1], [1]], [length_distr_single, rarefraction_single, logo_plot_single])
        else:
            for i in [length_distr_single, rarefraction_single, logo_plot_single]:
                if i != None:
                    Builder.add_figure(i)
        Builder.add_page()

    ls_clusters_dir = os.path.join(plot_path, "sequence_cluster")
    Builder.create_headline("Sequence Clusters")
    Builder.md_to_text(os.path.join("settings", "levenshtein_cluster.md"))
    Builder.md_to_text(os.path.join("settings", "dendrogram_cluster.md"))
    Builder.md_to_text(os.path.join("settings", "sequence_embedding.md"))
    
    report_seq_cluster = os.path.join(ls_clusters_dir, "reports")
    Builder.create_headline("Levenshtein distance based", section="##")
    for sample in samples:
        
        plot_files = find_multiple_files_with_substring(ls_clusters_dir, sample)
        if check_path_multiple(plot_files):
            Builder.create_headline(f"Sequence Similarity based on Levenshtein Distance of {sample}", "###")
            Builder.add_figure_horizontal(ncol = 2, picture_paths=plot_files, figure_description="")
        else:
            headline_made = False
            Builder.create_headline(f"Sequence Similarity based on Levenshtein Distance of {sample}", "###")
            for i in plot_files:
                if i != None:
                    Builder.add_figure(i)
        try:
            report_file = find_file_with_substring(report_seq_cluster, sample)
            if report_file != None:
                Builder.add_table(report_file, top_ = 10)
        except:
            pass

        Builder.add_page()

    Builder.create_headline("Sequence Similarity between samples", section = "##")
    embedding_dir = os.path.join(plot_path, "sequence_embedding")
    for sample in samples:
        file = find_file_with_substring(embedding_dir, sample)
        if file != None:
            Builder.create_headline(f"Global Sequence similarity for {sample}", section = "###")
            Builder.add_figure(file)
        else:
            pass
    Builder.write_quarto(save_dir=plot_path)
    if binding_data is not None:
        assert os.path.isdir(os.path.join(plot_path, "clustering_antigens")), f"The directory {os.path.join(plot_path, 'clustering_antigens')} does not exist"
        assert os.path.isdir(os.path.join(plot_path, "clustering_antigens", "ls_binding_cluster")), f"The directory {os.path.join(plot_path, 'clustering_antigens', 'ls_binding_cluster')} does not exist"
        assert os.path.isdir(os.path.join(plot_path, "clustering_antigens", "dendro_binding")), f"The directory {os.path.join(plot_path, 'clustering_antigens', 'dendro_binding')} does not exist"
        assert os.path.isdir(os.path.join(plot_path, "clustering_antigens", "reports")), f"The directory {os.path.join(plot_path, 'clustering_antigens', 'reports')} does not exist"
        Builder.create_headline("Cluster Antigens")
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
        pass
    Builder.write_quarto(save_dir=plot_path)

