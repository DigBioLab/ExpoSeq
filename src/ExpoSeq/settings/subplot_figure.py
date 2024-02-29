from .layout_finder import best_layout
import matplotlib.pyplot as plt

class Subplotter:
    def __init__(self):
        self.current_figure = 1
        self.all_figures = []
        self.all_axes = []
        setattr(self, f"fig{self.current_figure}", None)
        setattr(self, f"axes_all{self.current_figure}", None)
        
    
    def create_figure(self, figure = None):
        if figure == None:
            used_figure = self.current_figure
            
        else:
            used_figure = figure
        all_axes = getattr(self, f"axes{self.current_figure}")
        rows, cols = best_layout(len(all_axes))
        fig, axes = plt.subplots(nrows=rows, ncols=cols)
        setattr(self, f"fig{used_figure}", fig)
        setattr(self, f"axes_all{used_figure}", axes)
        self.all_figures.append(getattr(self, f"fig{used_figure}"))
        self.all_axes.append(getattr(self, f"axes_all{used_figure}"))
        
    
    def set_attribute(self, figure):
        if not hasattr(self, f"axes{figure}"):
            setattr(self, f"axes{figure}", [])
        else:
            pass
        
    def add_to_subplots(self, ax, figure = 1,):
        if self.current_figure == figure:
            pass
        else:
            self.current_figure = figure
        self.set_attribute(figure)
        axes_list = getattr(self, f"axes{figure}")
        axes_list.append(ax)
        
     
     
SubPlots = Subplotter()
     
import numpy as np

# Data for the plots
x = np.linspace(0, 10, 100)
y1 = np.sin(x)
y2 = np.cos(x)
y3 = x ** 2   

        
fig1, ax1 = plt.subplots()
ax1.plot(x, y1)
ax1.set_title('Line Plot')
ax1.set_xlabel('x')
ax1.set_ylabel('y')

SubPlots.add_to_subplots(ax1)

# Plot 2: Scatter plot
fig2, ax2 = plt.subplots()
ax2.scatter(x, y2)
ax2.set_title('Scatter Plot')
ax2.set_xlabel('x')
ax2.set_ylabel('y')

SubPlots.add_to_subplots(ax2)

# Plot 3: Bar plot
fig3, ax3 = plt.subplots()
ax3.bar(x, y3)
ax3.set_title('Bar Plot')
ax3.set_xlabel('x')
ax3.set_ylabel('y')

SubPlots.add_to_subplots(ax3)
SubPlots.create_figure()
plt.show()
