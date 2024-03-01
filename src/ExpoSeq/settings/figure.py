import matplotlib.pyplot as plt
from . import plot_styler
import os
import matplotlib
import PyQt5

class MyFigure:
    def __init__(self, test, figure_style = "seaborn-v0_8-colorblind"):
        self.fig = plt.figure(1, figsize = (12, 10))
        self.ax = self.fig.gca()
        self.plot_type = "multi"
        self.figure_style = figure_style  
        self.test = test

    def check_fig(self, ):
        if not plt.fignum_exists(1):
            self.fig = plt.figure(1)

    def clear_fig(self):
        self.fig.clear()
        if self.plot_type == "single":
            self.ax = self.fig.gca()

    def update_plot(self):
        self.ax = self.fig.gca()
        self.use_style(self.figure_style)
        self.ax_visibility()
        self.style = plot_styler.PlotStyle(self.ax, self.plot_type)
        if matplotlib.get_backend() == "Qt5Agg" and self.test != True:
            figManager = plt.get_current_fig_manager()
            figManager.window.showMaximized()
        if self.test:
            self.tighten()
            plt.show(block = False)

        

    def tighten(self):
        plt.tight_layout()
    
    def set_backend(self, backend = "Qt5Agg"):
        if self.test != True:
            plt.switch_backend(backend)
        #return
        
    @staticmethod
    def use_style(style = "seaborn-v0_8-colorblind"):
        assert style in plt.style.available, f'The style you chose is not availale. Choose one of {plt.style.available}' 
        plt.style.use(style)
        
    def ax_visibility(self):
        self.ax.spines['right'].set_visible(False) 
        self.ax.spines['top'].set_visible(False)


def save_matrix(matrix, path = None):
    if path == None:
        while True:
            save_matrix = input("Do you want to save the generated data? (Y/n)")
            if save_matrix in ["Y", "y", "N", "n"]:
                break
            else:
                print("Please enter Y or n")
        if save_matrix.lower() in ["Y", "y"]:
            while True:
                filename_matrix = input("Enter a name for the file. The file will be saved locally in your IDE.")
                if not os.path.isfile(filename_matrix):
                    matrix.to_excel(path + ".xlsx")
                    break
                else:
                    print("This file already exists. Please choose another name.")
    else:
        matrix.to_excel(path)



class Test:
    def __init__(self) -> None:
        Fig = MyFigure()
        Fig.set_backend()
        test_ax = Fig.ax
        Fig.clear_fig()
        assert Fig.ax == test_ax
        Fig.ax.scatter(x = [1,2,3], y = [3,4,5])
        Fig.plot_type = "single"
        Fig.clear_fig()
        assert Fig.ax == Fig.fig.gca()
        
