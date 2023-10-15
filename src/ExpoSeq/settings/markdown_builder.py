import os
import pandas as pd


class QuartoBuilder:
    def __init__(self, title):
        self.content = ""
        self.title = title
        self.create_header()

    def write_quarto(self):
        with open(f"{self.title}.qmd", "w") as f:
            f.write(self.content)
            
    def create_header(self):
        self.content += "---"
        self.content += f"title: {self.title}"
        self.content += f"format:"
        self.content += "\n html:"
        self.content += "\n     toc: true"
        self.content += "\n     toc-location: true:"
        self.content += "\n     toc-expand: 2"
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
        self.next_row()
        self.content += f"![]"
        self.content += f"({picture_path})"
        self.next_row()
        self.content += ":::"
        
    def add_table(self, filename, top_,):
        self.next_row()
        df = pd.read_excel(filename).head(top_)
        # Get the headers
        headers = df.columns.tolist()

        # Convert headers to markdown format
        markdown_headers = ' | '.join(headers)
        markdown_divider = ' | '.join(['---'] * len(headers))
        self.content += f"{markdown_headers}\n{markdown_divider}"
        
        
        

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
    for root, dirs, files in os.walk(directory):
        for file in files:
            if substring in file:
                matching_files.append(os.path.join(root, file))
    return matching_files
    
plot_path = r"C:\Users\nilsh\my_projects\ExpoSeq\my_experiments\test_all_samples\plots\aaSeqCDR3"
# samples have to come from origin seq report
experiment = "test"
samples = ['NG-28545_Library_1_F10_2', 'NG-28545_Library_1_F11_2', 'NG-28545_Library_1_F12_2', 'NG-28545_Library_1_F1_2', 'NG-28545_Library_1_F2_2', 'NG-28545_Library_1_F3_2', 'NG-28545_Library_1_F5_2', 'NG-28545_Library_1_F6_2', 'NG-28545_Library_1_F7_2', 'NG-28545_Library_1_F8_2', 'NG-28545_Library_1_F9_2', 'NG-28545_Library_2_F11_2', 'NG-28545_Library_2_F12_2', 'NG-28545_Library_2_F2_2', 'NG-28545_Library_2_F3_2', 'NG-28545_Library_2_F5_2', 'NG-28545_Library_2_F6_2', 'NG-28545_Library_2_F8_2', 'NG-28545_Library_2_F9_2', 'NG-28545_Library_3_F11_2', 'NG-28545_Library_3_F2_2', 'NG-28545_Library_3_F3_2', 'NG-28545_Library_3_F7_2', 'NG-28545_Library_4_F1_2', 'NG-28545_Library_4_F3_2', 'NG-28545_Library_4_F5_2', 'NG-28545_Library_4_F6_2']
Builder = QuartoBuilder(experiment)
Builder.create_headline("Basic Overview")
Builder.create_headline("Identity between samples", section = "##")
Builder.add_figure_horizontal(2, [r"C:\Users\nilsh\my_projects\ExpoSeq\my_experiments\test_all_samples\plots\aaSeqCDR3\jaccard.png", r"C:\Users\nilsh\my_projects\ExpoSeq\my_experiments\test_all_samples\plots\aaSeqCDR3\morosita_horn.png"], "test description")
Builder.add_page()
Builder.create_headline("Alignment Quality of samples", section = "##")
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
report_seq_cluster = os.path.isdir(os.path.join(plot_path, "sequence_cluster", "reports"))
for sample in samples:
    Builder.create_headline(f"Sequence Similarity based on Levenshtein Distance of {sample}", "###")
    plot_files = find_multiple_files_with_substring(ls_clusters_dir, sample)
    if check_path_multiple(plot_files):
        Builder.add_figure_horizontal(ncol = 2, picture_paths=plot_files, figure_description="")
    else:
        for i in plot_files:
            if i != None:
                Builder.add_figure(i)
    try:
        report_seq_cluster = find_file_with_substring(report_seq_cluster, sample)
        Builder.add_table(report_seq_cluster, top_ = 20)
    except:
        pass
    Builder.add_page()

Builder.create_headline("Sequence Similarity between samples")
embedding_dir = os.path.join(plot_path, "sequence_embedding")
for sample in samples:
    Builder.create_headline(f"Global Sequence similarity for {sample}")
    file = find_file_with_substring()



Builder.write_quarto()
