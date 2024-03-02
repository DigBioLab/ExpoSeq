import glob
import os
from .markdown_builder import QuartoBuilder
import subprocess


class Subplotter:
    def __init__(self,  figure_title):
        self.check_base_dir()
        self.files = []
        self.figure_title = figure_title
        self.captures = []
        
    @staticmethod
    def check_base_dir():
        if os.path.isdir("tmp_quarto"):
            pass
        else:
            os.mkdir("tmp_quarto")
            
            
    def update_files(self, figure_dir):
        files_raw = glob.glob(os.path.join(figure_dir, "*png"))
        self.files = []
        for file in files_raw:
            self.files.append(os.path.abspath(file))
            
    def update_capture(self, fig_capture):
        self.captures.append(fig_capture)
        
    def add_as_subplot(self, fig, fig_capture = "", dir = "tmp_quarto" ):
        self.update_capture(fig_capture)
        if dir == "tmp_quarto":
            figure_dir = os.path.join(dir, self.figure_title)
        else:
            figure_dir = dir
            assert os.path.isdir(figure_dir)
        if os.path.isdir(figure_dir):
            pass
        else:
            os.makedirs(figure_dir)
        no = str(len(self.files))
        fig.savefig(os.path.join(figure_dir, f"{no}_.png"),
                                                     format = "png",
                                                     dpi = 300, 
                                                     bbox_inches='tight')
        self.update_files(figure_dir)
        
    def make_figure(self, dir = "tmp_quarto", ):
        Builder = QuartoBuilder(self.figure_title, figure=True)  
        Builder.add_subplot_figures(self.files, self.captures)
        
        if dir == "tmp_quarto":
            qmd_file = os.path.join(dir, self.figure_title, f"{self.figure_title}.qmd")
            Builder.write_quarto(save_dir=qmd_file)
        else:
            qmd_file = dir
            Builder.write_quarto(save_dir=qmd_file)
        subprocess.run(["quarto", "render", os.path.join(qmd_file, Builder.title + ".qmd")])




