import matplotlib.pyplot as plt
from ast import literal_eval

class PlotStyle:
    def __init__(self, ax, ):
        self.ax = ax
        with open('font_settings.txt', "r") as f:
            font_settings = f.read()
        self.font_settings = literal_eval(font_settings)
        with open('legend_settings.txt', "r") as f:
            legend_settings = f.read()
        self.legend_settings = literal_eval(legend_settings)

    def title_figure(self, title):
        self.ax.set_title(title)

    def title_xaxis(self, xlabel):
        self.ax.set_xlabel(xlabel)
    def title_yaxis(self, ylabel):
        self.ax.set_ylabel(ylabel)
    def xaxis_limit(self, xmin, xmax):
        self.ax.set_xlim(xmin, xmax)
    def yaxis_limit(self, ymin, ymax):
        self.ax.set_ylim(ymin, ymax)
    def xscale(self, ticks):
        self.ax.set_xticks(ticks)
    def yscale(self, ticks):
        self.ax.set_yticks(ticks)
    def xticklabels(self, labels):
        self.ax.set_xticklabels(labels)
    def yticklabels(self, labels):
        self.ax.set_yticklabels(labels)
    def xaxis_log(self):
        self.ax.set_xscale("log")
    def yaxis_log(self):
        self.ax.set_yscale("log")
    def xaxis_linear(self):
        self.ax.set_xscale("linear")
    def yaxis_linear(self):
        self.ax.set_xscale("linear")
    def grid(self, visible=True):
        self.ax.grid(visible)
    def style(self, style):
        plt.style.use(style)
    def facecolor(self, color):
        self.ax.set_facecolor(color)

    def annotate_coordinate(self, x, y, text, **kwargs):
        self.ax.annotate(text,
                         xy=(x, y),
                         xytext=(x, y),
                         **kwargs)
    def annotate_coordinate_with_arrow(self, x, y, text, **kwargs):
        self.ax.annotate(text,
                         xy=(x, y),
                         xytext=(x, y),
                         arrowprops=dict(facecolor='red', arrowstyle="->"),
                         **kwargs)
    def font_family(self, font_name, global_change = False):
        plt.rc('font',
               family=font_name)
        if global_change == True:
            self.font_settings["family"] = font_name
    def font_style(self, style, global_change = False):
        plt.rc(weight = style)
        if global_change == True:
            self.font_settings["fontweight"] = style
    def font_color(self, color):
        for item in ([self.ax.title, self.ax.xaxis.label, self.ax.yaxis.label] +
                     self.ax.get_xticklabels() + self.ax.get_yticklabels()):
            item.set_color(color)
    def titlesize_figure(self, size):
        self.ax.title.set_size(size)

    def xy_titlesize(self, size, global_change = False):
        self.ax.xaxis.label.set_size(size)
        self.ax.yaxis.label.set_size(size)
        if global_change == True:
            self.font_settings["fontsize"] = size
