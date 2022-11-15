import logomaker
from python_scripts.plots.plot_params.open_txtfiles import openParams


logo_plot = logomaker.Logo(aa_distribution,
                            shade_below = .5,
                            fade_below = .5,
                            font_name = 'Arial Rounded MT Bold')
plot_style = openParams('plot_style.txt')
logo_plot.ax.set_ylabel("Amino Acid Frequency", **plot_style)
logo_plot.ax.set_xlabel("Sequence Position", **plot_style)
logo_plot.ax.set_xticks(range(1, longest_sequence))
logo_plot.highlight_position(p = 5, color = 'gold', alpha = .5)
logo_plot.highlight_position_range(pmin = 3, pmax = 5, color = "lightcyan")