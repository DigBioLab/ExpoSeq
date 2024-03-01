from .layout_finder import best_layout
import matplotlib.pyplot as plt

class Subplotter:
    def __init__(self):
        self.current_figure = 1
        self.all_figures = []
        self.all_axes = []
        setattr(self, f"fig{self.current_figure}", None)
        setattr(self, f"axes_all{self.current_figure}", None)
        
    

        
    def update_fig(self, figure, axes_list):
        rows, cols = best_layout(len(axes_list))
        fig, axes = plt.subplots(nrows=rows, ncols=cols)
        for ax in fig.axes:
            fig.delaxes(ax)

        # Add new axes from axes_list
        for ax in axes_list:
            fig._axstack.add(fig._make_key(ax), ax)
        # modify class attributes
        if hasattr(self, f"fig{figure}"):
            given_fig = getattr(self, f"fig{figure}")
        else:
            setattr(self, f"fig{figure}", fig)
            setattr(self, f"axes_all{figure}", axes)
            given_fig = getattr(self, f"fig{figure}")
            
        if given_fig in self.all_figures:
            pass
        else:
            self.all_figures.append(getattr(self, f"fig{figure}"))
   #     self.all_axes.append(getattr(self, f"axes_all{figure}"))

        
    
    
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
        self.update_fig(figure, axes_list)
        
     
     
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
