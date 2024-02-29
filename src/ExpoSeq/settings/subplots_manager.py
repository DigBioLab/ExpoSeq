import glob
import os
from .markdown_builder import QuartoBuilder
import subprocess


class Subplotter:
    def __init__(self,  figure_title):
        self.check_base_dir()
        self.files = []
        self.figure_title = figure_title
        
    @staticmethod
    def check_base_dir():
        if os.path.isdir("tmp_quarto"):
            pass
        else:
            os.mkdir("tmp_quarto")
            
        
    def add_as_subplot(self, fig ):
        figure_dir = os.path.join("tmp_quarto", self.figure_title)
        if os.path.isdir(figure_dir):
            pass
        else:
            os.makedirs(figure_dir)
        no = str(len(self.files))
        
        fig.savefig(os.path.join(figure_dir, f"{no}_.png"),
                                                     format = "png",
                                                     dpi = 300, bbox_inches='tight')
        files_raw = glob.glob(os.path.join(figure_dir, "*png"))
        self.files = []
        for file in files_raw:
            self.files.append(os.path.abspath(file))
        
        
    def make_figure(self):
        Builder = QuartoBuilder(self.figure_title, figure=True)  
        Builder.add_subplot_figures(self.files)
        Builder.write_quarto(save_dir=os.path.join("tmp_quarto", self.figure_title, ))
        subprocess.run(["quarto", "render", os.path.join("tmp_quarto", self.figure_title, f"{self.figure_title}.qmd")])
        for png in self.files:
            os.remove(png)
        
    



