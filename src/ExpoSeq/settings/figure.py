import matplotlib.pyplot as plt
from . import plot_styler
import os

class MyFigure:
    def __init__(self):
        self.fig = plt.figure(1, figsize = (12, 10))
        self.ax = self.fig.gca()
        self.plot_type = "multi"

    def check_fig(self, ):
        if not plt.fignum_exists(1):
            self.fig = plt.figure(1)

    def clear_fig(self):
        self.fig.clear()
        if self.plot_type == "single":
            self.ax = self.fig.gca()

    def update_plot(self):
        self.ax = self.fig.gca()
        self.style = plot_styler.PlotStyle(self.ax, self.plot_type)
        figManager = plt.get_current_fig_manager()
        figManager.window.showMaximized()
        self.tighten()

    def tighten(self):
        plt.tight_layout()


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

