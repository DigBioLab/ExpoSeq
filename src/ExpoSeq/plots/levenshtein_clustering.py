import matplotlib.pyplot as plt
from numpy import unique
import networkx as nx
import math
import pandas as pd
from ExpoSeq.tidy_data.clustering import cleaning
import community
from ExpoSeq.plots.layout_finder import best_layout

def calculate_nodesize(sample_report, region_string, G):
    nodesize = []
    label_numbers = {}
    label_sequences = {}
    for index, g in enumerate(G):
        for i in sample_report[region_string]:
            if i == g:
                node = math.sqrt(sample_report.loc[sample_report[region_string] == i, "cloneFraction"].values[0])
                nodesize.append(node)
                if node * 500 > 1 / 100 * 500:
                    label_numbers[g] = index
                    label_sequences[g] = g
                else:
                    label_numbers[g] = ""
                    label_sequences[g] = ""
                break
    for y, n in enumerate(nodesize):
        nodesize[y] = int(n * 500)

    return nodesize, label_numbers, label_sequences


def clusterSeq(ax, sequencing_report, samples: list, max_ld, min_ld, batch_size, font_settings, second_figure,
               region_string, label_type="numbers"):
    sample_report = sequencing_report.loc[sequencing_report["Experiment"].isin(samples)]
    sample_report = sample_report.groupby("Experiment").head(batch_size)
    G = cleaning(sample_report, max_ld, min_ld, region_string)
    # remove nodes that are not connected to any other nodes
    G_deg = G.degree()
    to_remove = [n for (n, deg) in G_deg if deg == 0]
    G.remove_nodes_from(to_remove)

    # Create a gridspec for adding subplots of different sizes
    nodesize, label_numbers, label_sequences = calculate_nodesize(sample_report, region_string, G)

    ## labels
    label_sequences = {}
    node_ids = list(G.nodes())
    n = 0
    for index, g in enumerate(G):
        if n == 9:
            label_sequences[g] = g
            n = 0
        else:
            label_sequences[g] = ""
        n += 1

    partition = community.best_partition(G)
    pos = nx.spring_layout(G)
    nx.draw_networkx_nodes(G, pos, node_size=nodesize, cmap=plt.cm.RdYlBu, node_color=list(partition.values()), ax=ax)
    nx.draw_networkx_edges(G, pos, alpha=0.3, ax=ax)
    if label_type == "numbers":
        nx.draw_networkx_labels(G, pos,
                                labels=label_numbers,
                                alpha=0.5,
                                font_color="grey",
                                horizontalalignment="right",
                                verticalalignment="bottom", ax=ax)
    elif label_type == "sequence":
        nx.draw_networkx_labels(G, pos,
                                labels=label_sequences,
                                alpha=0.5,
                                font_color="grey",
                                horizontalalignment="right",
                                verticalalignment="bottom", ax=ax)

    original_fontsize = font_settings["fontsize"]
    font_settings["fontsize"] = 22
    ax.set_title(f"Connected components of {samples}", **font_settings)
    font_settings["fontsize"] = original_fontsize
    ax.set_axis_off()

    ## create report

    # Convert each RGBA to nearest named color
    # color_names = [closest_color((r, g, b)) for r, g, b, _ in rgba_colors]
    # color_names = hex_list_to_names(color_names)
    cluster_report = pd.DataFrame([partition]).T
    cluster_report['Sequences'] = cluster_report.index
    cluster_report.reset_index(drop=True, inplace=True)
    cluster_report.rename(columns={0: 'Cluster No.'}, inplace=True)

    return cluster_report


def cluster_single_AG(fig, sequencing_report, antigen, binding_data, max_ld, min_ld, batch_size, region_string,
                      preferred_cmap="Blues", specific_experiments=False, label_type="numbers"):
    assert antigen in binding_data.columns, "it seems like your antigen does not exist in the binding data. Please enter the correct value"

    if specific_experiments != False:
        report_batch = sequencing_report[sequencing_report['Experiment'].isin(specific_experiments)]
    else:
        report_batch = sequencing_report
    report_batch = report_batch.groupby("Experiment").head(batch_size)
    experiments = report_batch["Experiment"].unique()
    Tot = experiments.shape[0]
    Rows, Cols = best_layout(Tot)
    Position = range(1, Tot + 1)
    plot_number = 0
    cluster_report = pd.DataFrame([])
    for experiment in experiments:
        ax = fig.add_subplot(Rows, Cols, Position[plot_number])
        batch = report_batch[report_batch["Experiment"] == experiment]
        batch = batch.dropna(subset = ["aaSeqCDR3"])
        binding_data = binding_data[[antigen, "aaSeqCDR3"]]
        binding_data = binding_data.dropna(subset = ["aaSeqCDR3"])
        mix = pd.merge(batch, binding_data, how = "left", on = "aaSeqCDR3")
        mix = mix.fillna(0)
        mix = mix.reset_index()
        G = cleaning(mix, max_ld, min_ld, region_string)
        # remove nodes that are not connected to any other nodes
        G_deg = G.degree()
        to_remove = [n for (n, deg) in G_deg if deg == 0]
        G.remove_nodes_from(to_remove)
       #fr print(mix.)
        nodesize = []
        clone_counts = mix["readCount"]
        color_values = []
        label_numbers = {}
        label_sequences = {}
        for index, g in enumerate(G):
            for i in mix["aaSeqCDR3"]:
                if i == g:
                    node = math.sqrt(mix.loc[mix[region_string] == i, "cloneFraction"].values[0])
                    nodesize.append(node)
                    color_values.append(mix.loc[mix[region_string] == i, antigen].values[0])
                    if node * 500 > 1 / 100 * 500:
                        label_numbers[g] = index
                        label_sequences[g] = g
                    else:
                        label_numbers[g] = ""
                        label_sequences[g] = ""
                    break
        for y, n in enumerate(nodesize):
            nodesize[y] = int(n * 500)



        cm = plt.get_cmap(preferred_cmap)
        vmin = 0
        vmax = max(color_values)
        print(max(color_values))

        partition = community.best_partition(G)
        pos = nx.spring_layout(G)

        nodes = nx.draw_networkx_nodes(G,
                                       pos,
                                       node_size=nodesize,
                                       cmap=cm,
                                       node_color=color_values,
                                       ax=ax,
                                       vmin=vmin,
                                       vmax=vmax)
        nx.draw_networkx_edges(G, pos, alpha=0.3)
        if label_type == "numbers":
            nx.draw_networkx_labels(G, pos,
                                    labels=label_numbers,
                                    alpha=0.5,
                                    font_color="grey",
                                    horizontalalignment="right",
                                    verticalalignment="bottom", ax=ax)
        elif label_type == "sequence":
            nx.draw_networkx_labels(G, pos,
                                    labels=label_sequences,
                                    alpha=0.5,
                                    font_color="grey",
                                    horizontalalignment="right",
                                    verticalalignment="bottom", ax=ax)

        cluster_report_inter = pd.DataFrame([partition]).T
        cluster_report_inter[f'Sequences_{experiment}'] = cluster_report_inter.index
        cluster_report_inter.reset_index(drop=True, inplace=True)
        cluster_report_inter.rename(columns={0: f'Clusters_{experiment}_No.'}, inplace=True)
        cluster_report = pd.concat([cluster_report, cluster_report_inter])
        ax.set_title(experiment)
        sm = plt.cm.ScalarMappable(cmap=preferred_cmap,
                                   norm=plt.Normalize(vmin=vmin, vmax=vmax))
        plt.colorbar(sm, ax = ax)
        plot_number = plot_number + 1
  #  fig.colorbar(nodes)
    return cluster_report


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
    G, degree_sequence = cleaning(mix, max_ld, min_ld, region_string)
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
            if n > 0 and n < 15:
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

# fig.colorbar(ax, col_ax)















