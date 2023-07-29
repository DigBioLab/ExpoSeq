import matplotlib.pyplot as plt
from ast import literal_eval
import os
import pkg_resources
class PlotStyle:
    def __init__(self, ax, plot_type):
        self.pkg_path = pkg_resources.resource_filename("ExpoSeq", "")
        self.module_dir = os.path.abspath("")
        self.ax = ax
        font_settings_path = os.path.join(self.pkg_path,
                                          "settings",
                                          "font_settings.txt")
        with open(font_settings_path, "r") as f:
            font_settings = f.read()
        self.font_settings = literal_eval(font_settings)
        legend_settings_path = os.path.join(self.pkg_path,
                                            "settings",
                                            "legend_settings.txt")
        with open(legend_settings_path, "r") as f:
            legend_settings = f.read()
        self.legend_settings = literal_eval(legend_settings)
        self.plot_type = plot_type
    def title_figure(self, title):
        if self.plot_type == "single":
            self.ax.set_title(title)
        else:
            print("Unfortunately, you can't change the layout of subplots. You can try to change it using the matplotlib.pyplot package")
    def title_xaxis(self, xlabel):
        if self.plot_type == "single":
            self.ax.set_xlabel(xlabel)
        else:
            print("Unfortunately, you can't change the layout of subplots. You can try to change it using the matplotlib.pyplot package")
    def title_yaxis(self, ylabel):
        if self.plot_type == "single":
            self.ax.set_ylabel(ylabel)
        else:
            print("Unfortunately, you can't change the layout of subplots. You can try to change it using the matplotlib.pyplot package")
    def xaxis_limit(self, xmin, xmax):
        if self.plot_type == "single":
            self.ax.set_xlim(xmin, xmax)
        else:
            print("Unfortunately, you can't change the layout of subplots. You can try to change it using the matplotlib.pyplot package")
    def yaxis_limit(self, ymin, ymax):
        if self.plot_type == "single":
            self.ax.set_ylim(ymin, ymax)
        else:
            print("Unfortunately, you can't change the layout of subplots. You can try to change it using the matplotlib.pyplot package")
    def xscale(self, ticks):
        if self.plot_type == "single":
            self.ax.set_xticks(ticks)
        else:
            print("Unfortunately, you can't change the layout of subplots. You can try to change it using the matplotlib.pyplot package")
    def yscale(self, ticks):
        if self.plot_type == "single":
            self.ax.set_yticks(ticks)
        else:
            print("Unfortunately, you can't change the layout of subplots. You can try to change it using the matplotlib.pyplot package")
    def xticklabels(self, labels):
        if self.plot_type == "single":
            self.ax.set_xticklabels(labels)
        else:
            print("Unfortunately, you can't change the layout of subplots. You can try to change it using the matplotlib.pyplot package")
    def yticklabels(self, labels):
        if self.plot_type == "single":
            self.ax.set_yticklabels(labels)
        else:
            print("Unfortunately, you can't change the layout of subplots. You can try to change it using the matplotlib.pyplot package")
    def xaxis_log(self):
        if self.plot_type == "single":
            self.ax.set_xscale("log")
        else:
            print("Unfortunately, you can't change the layout of subplots. You can try to change it using the matplotlib.pyplot package")
    def yaxis_log(self):
        if self.plot_type == "single":
            self.ax.set_yscale("log")
        else:
            print("Unfortunately, you can't change the layout of subplots. You can try to change it using the matplotlib.pyplot package")
    def xaxis_linear(self):
        if self.plot_type == "single":
            self.ax.set_xscale("linear")
        else:
            print("Unfortunately, you can't change the layout of subplots. You can try to change it using the matplotlib.pyplot package")
    def yaxis_linear(self):
        if self.plot_type == "single":
            self.ax.set_xscale("linear")
        else:
            print("Unfortunately, you can't change the layout of subplots. You can try to change it using the matplotlib.pyplot package")
    def grid(self, visible=True):
        if self.plot_type == "single":
            self.ax.grid(visible)
        else:
            print("Unfortunately, you can't change the layout of subplots. You can try to change it using the matplotlib.pyplot package")
    def style(self, style):
        if self.plot_type == "single":
            plt.style.use(style)
        else:
            print("Unfortunately, you can't change the layout of subplots. You can try to change it using the matplotlib.pyplot package")
    def facecolor(self, color):
        if self.plot_type == "single":
            self.ax.set_facecolor(color)
        else:
            print("Unfortunately, you can't change the layout of subplots. You can try to change it using the matplotlib.pyplot package")

    def annotate_coordinate(self, x, y, text, **kwargs):
        if self.plot_type == "single":
            self.ax.annotate(text,
                         xy=(x, y),
                         xytext=(x, y),
                         **kwargs)
        else:
            print("Unfortunately, you can't change the layout of subplots. You can try to change it using the matplotlib.pyplot package")
    def annotate_coordinate_with_arrow(self, x, y, text, **kwargs):
        if self.plot_type == "single":
            self.ax.annotate(text,
                         xy=(x, y),
                         xytext=(x, y),
                         arrowprops=dict(facecolor='red', arrowstyle="->"),
                         **kwargs)
        else:
            print("Unfortunately, you can't change the layout of subplots. You can try to change it using the matplotlib.pyplot package")
    def font_family(self, font_name, global_change = False):
        if self.plot_type == "single":
            plt.rc('font',
                family=font_name)
            if global_change == True:
                self.font_settings["family"] = font_name
        else:
            print(
                "Unfortunately, you can't change the layout of subplots. You can try to change it using the matplotlib.pyplot package")
    def font_style(self, style, global_change = False):
        if self.plot_type == "single":
            plt.rc(weight = style)
            if global_change == True:
                self.font_settings["fontweight"] = style
        else:
            print("Unfortunately, you can't change the layout of subplots. You can try to change it using the matplotlib.pyplot package")
    def font_color(self, color):
        if self.plot_type == "single":
            for item in ([self.ax.title, self.ax.xaxis.label, self.ax.yaxis.label] +
                         self.ax.get_xticklabels() + self.ax.get_yticklabels()):
                item.set_color(color)
        else:
            print("Unfortunately, you can't change the layout of subplots. You can try to change it using the matplotlib.pyplot package")
    def titlesize_figure(self, size):
        if self.plot_type == "single":
            self.ax.title.set_size(size)
        else:
            print("Unfortunately, you can't change the layout of subplots. You can try to change it using the matplotlib.pyplot package")

    def xy_titlesize(self, size, global_change = False):
        if self.plot_type == "single":
            self.ax.xaxis.label.set_size(size)
            self.ax.yaxis.label.set_size(size)
            if global_change == True:
                self.font_settings["fontsize"] = size
        else:
            print("Unfortunately, you can't change the layout of subplots. You can try to change it using the matplotlib.pyplot package")
