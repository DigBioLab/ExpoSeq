import matplotlib.pyplot as plt
from numpy import unique
import networkx
import math
import pandas as pd
from ExpoSeq.tidy_data.clustering import cleaning


def clusterSeq(ax, sequencing_report, sample,max_ld, min_ld, batch_size, font_settings, second_figure):
    sample_report = sequencing_report[sequencing_report["Experiment"] == sample] ## insert test if sample not found
    sample_report = sample_report.head(batch_size)
    G, degree_sequence = cleaning(sample_report, max_ld, min_ld)
    # Create a gridspec for adding subplots of different sizes
    clone_counts = sample_report["readCount"]
    G_deg = G.degree()
    to_remove = [n for (n, deg) in G_deg if deg == 0]
    G.remove_nodes_from(to_remove)
    nodesize = []
    for i in sample_report.index:
        for g in G:
            if i == g:
                node = math.sqrt(clone_counts.iloc[i])
                nodesize.append(node)
    for y, n in enumerate(nodesize):
        if n > 0 and n < 15:
            nodesize[y] = 15
        if n == 0:
            nodesize[y] = 50
    networkx.draw_networkx(G, arrows = True, with_labels = False,ax = ax)
    original_fontsize = font_settings["fontsize"]
    font_settings["fontsize"] = 22
    ax.set_title("Connected components of " + sample, **font_settings)
    font_settings["fontsize"] = original_fontsize
    ax.set_axis_off()
    if second_figure == True:
        fig2 = plt.figure()
        ax2 = fig2.gca()
        ax2.bar(*unique(degree_sequence, return_counts=True))
        original_fontsize = font_settings["fontsize"]
        font_settings["fontsize"] = 22
        ax2.set_title(sample + ' Histogram', **font_settings)
        font_settings["fontsize"] = original_fontsize
        ax2.set_xlabel("Degree", **font_settings)
        ax2.set_ylabel("# of Nodes", font_settings)
        return fig2
    else:
        return

def cluster_single_AG(fig, sequencing_report, antigen, binding_data, max_ld, min_ld, batch_size, specific_experiments = False, ):
    report_batch = sequencing_report.groupby("Experiment").head(batch_size)
    if specific_experiments != False:
        report_batch = report_batch[report_batch['Experiment'].isin(specific_experiments)]
    else:
        pass
    experiments = report_batch["Experiment"].unique()
    Tot = experiments.shape[0]
    Cols = int(input("How many Columns for the Plots do you want?"))
    Rows = Tot // Cols
    if Tot % Cols != 0:
        Rows += 1
    # Create a Position index
    Position = range(1, Tot + 1)
    plot_number = 0

    for experiment in experiments:
        ax = fig.add_subplot(Rows, Cols, Position[plot_number])
        batch = report_batch[report_batch["Experiment"] == experiment]
        mix = pd.concat([batch, binding_data])
        mix = mix.fillna(0)
        mix = mix.reset_index()
        G, degree_sequence = cleaning(mix, max_ld, min_ld)
        G_deg = G.degree()
        to_remove = [n for (n, deg) in G_deg if deg == 0]
        G.remove_nodes_from(to_remove)
        nodesize = []
        clone_counts = mix["readCount"]
        color_values = []
        for i in mix.index:
            for g in G:
                if i == g:
                    node = math.sqrt(clone_counts.iloc[i])
                    nodesize.append(node)
                    color_values.append(mix[antigen].iloc[i])
        for y, n in enumerate(nodesize):
            if n>0 and n < 15:
                nodesize[y] = 15
            if n == 0:
                nodesize[y] = 50
        cm = plt.get_cmap("Blues")
        vmin = min(color_values)
        vmax = max(color_values)
        network_plot = networkx.draw_networkx(G,
                                              node_size = nodesize,
                                              node_color = color_values,
                                              cmap = cm,
                                              with_labels = False,
                                              arrows = True,
                                              edgecolors = "Black",
                                              ax = ax)
     #   net = networkx.draw_networkx(G, node_size=nodesize, node_color = color_values, cmap = cm,with_labels=False, edgecolors = "Black", ax = ax)
       # network_plot = networkx.draw_networkx_nodes(G, pos=networkx.spring_layout(G), node_size=nodesize, node_color = color_values, cmap = cm,edgecolors = "Black", vmin = vmin, vmax = vmax)
        ax.set_title(experiment)
     #   sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=vmin, vmax=vmax))
      #  cbar = plt.colorbar(sm)
        plot_number = plot_number+1
    fig.colorbar(network_plot,
                 ax=fig.axes)

antigens = ['Ecarpholin_S_DDEL_5ug/mL', 'Myotoxin_II_DDEL_5ug/mL',
           'PLA2_n.naja(Uniprot:P15445)_DDEL_5ug/mL', 'Ecarpholin_S_CDEL_100nM',
           'Myotoxin_II_CDEL_100nM',]



### is not in use
def cluster_antigens(sequencing_report, sample, antigens, batch_size):
    report_batch = sequencing_report.groupby("Experiment").head(batch_size)
    Tot = len(antigens)
    Cols = int(input("How many Columns for the Plots do you want?"))
    # Compute Rows required
    Rows = Tot // Cols
    #     EDIT for correct number of rows:
    #     If one additional row is necessary -> add one:
    if Tot % Cols != 0:
        Rows += 1
    # Create a Position index
    Position = range(1, Tot + 1)
    plot_number = 0
    fig = plt.figure(1)
    max_colorbar = 0
    batch = report_batch[report_batch["Experiment"] == sample]
    mix = pd.concat([batch, delfia_samples])
    mix = mix.fillna(0)
    mix = mix.reset_index()
    mix = mix.drop("index", axis=1)
    G, degree_sequence = cleaning(mix)
    G_deg = G.degree()
    to_remove = [n for (n, deg) in G_deg if deg == 0]
    G.remove_nodes_from(to_remove)
    for antigen in antigens:
        ax = fig.add_subplot(Rows, Cols, Position[plot_number])
        nodesize = []
        clone_counts = mix.readCount
        color_values = []

        for i in range(mix.shape[0]):
            for g in G:
                if i == g:
                    node = math.sqrt(clone_counts.iloc[i])
                    nodesize.append(node)
                    color_values.append(mix[antigen].iloc[i])
        for y, n in enumerate(nodesize):
            if n>0 and n<15:
                nodesize[y] = 15
            if n == 0:
                nodesize[y] = 50

        colormap = plt.get_cmap("Blues")
        vmin = min(color_values)
        vmax = max(color_values)
        if vmax > max_colorbar:
            max_colorbar = vmax
        #   net = networkx.draw_networkx(G, node_size=nodesize, node_color = color_values, cmap = cm,with_labels=False, edgecolors = "Black", ax = ax)
        net = networkx.draw_networkx_nodes(G,
                                           pos=networkx.spring_layout(G),
                                           node_size=nodesize,
                                           node_color=color_values,
                                           cmap=colormap,
                                     edgecolors="Black")
        ax.set_title(antigen)
        #   sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=0, vmax=max_colorbar))
        #  cbar = plt.colorbar(sm)
        plot_number = plot_number + 1
    fig.suptitle("Clustering for sample " + sample)
#    sm = plt.cm.ScalarMappable(cmap=colormap,
 #                              norm=plt.Normalize(vmin = vmin, vmax=max_colorbar))

    #fig.colorbar(ax, col_ax)















