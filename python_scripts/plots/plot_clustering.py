import matplotlib.pyplot as plt
from numpy import unique
from python_scripts.tidy_data.clustering import cleaning
from networkx import draw_networkx
import imnet

G, degree_sequence = cleaning(sample_name, sequencing_report_all)
fig = plt.figure(sample_name, figsize=(5, 5))
# Create a gridspec for adding subplots of different sizes
axgrid = fig.add_gridspec(5, 4)
ax0 = fig.add_subplot(axgrid[0:3, :])
G_deg = G.degree()
to_remove = [n for (n, deg) in G_deg if deg == 0]
G.remove_nodes_from(to_remove)
draw_networkx(G, node_size=10, with_labels=False)
ax0.set_title("Connected components of " + sample_name)
ax0.set_axis_off()
ax2 = fig.add_subplot(axgrid[3:, :])
ax2.bar(*unique(degree_sequence, return_counts=True))
ax2.set_title(sample_name + ' Histogram')
ax2.set_xlabel("Degree")
ax2.set_ylabel("# of Nodes")
plt.savefig(sample_name + '.png', dpi=300)





from python_scripts.tidy_data.clustering import cleaning
import matplotlib.pyplot as plt
from numpy import unique
from python_scripts.tidy_data.clustering import cleaning

import networkx
import pandas as pd
import math
from networkx import draw_networkx
DELFIA1 = pd.read_csv(r"C:\Users\nils\Desktop\Tropical Pharmaceutical Lab\Project\sanger\Chris_main_df.csv")
DELFIA2 = pd.read_csv(r"C:\Users\nils\Desktop\Tropical Pharmaceutical Lab\Project\sanger\Helen_project_df.csv")
DELFIA3 = pd.read_csv(r"C:\Users\nils\Desktop\Tropical Pharmaceutical Lab\Project\sanger\Line_project_df.csv")
DELFIA4 = pd.read_csv(r"C:\Users\nils\Desktop\Tropical Pharmaceutical Lab\Project\sanger\Chris_myo_project_df.csv")
delfia_samples = pd.concat([DELFIA1,DELFIA2,DELFIA3,DELFIA4])
delfia_samples = delfia_samples.rename(columns={'JUNCTION_H':'aaSeqCDR3'})

experiments = sequencing_report_all["Experiment"].unique()

#report_batch = sequencing_report_all[sequencing_report_all["Experiment"] == library]
report_batch = sequencing_report_all.groupby("Experiment").head(300)
#report_batch = report_batch[["Experiment", "cloneCount", "aaSeqCDR3"]]


def cluster_single_AG(antigens):

    Tot = experiments.shape[0]
    Cols = int(input("How many Columns for the Plots do you want?"))

    Rows = Tot // Cols

    if Tot % Cols != 0:
        Rows += 1
    # Create a Position index
    Position = range(1, Tot + 1)
    plot_number = 0
    fig = plt.figure(1, constrained_layout=True)

    for experiment in experiments:
        ax = fig.add_subplot(Rows, Cols, Position[plot_number])
        batch = report_batch[report_batch["Experiment"] == experiment]
        mix = pd.concat([batch, delfia_samples])
        mix = mix.fillna(0)
        mix = mix.reset_index()
        G, degree_sequence = cleaning(experiment, mix)
        G_deg = G.degree()
        to_remove = [n for (n, deg) in G_deg if deg == 0]

        G.remove_nodes_from(to_remove)
        nodesize = []
        clone_counts = mix.cloneCount
        color_values = []
        for i in mix.index:
            for g in G:
                if i == g:
                    node = math.sqrt(clone_counts.iloc[i])
                    nodesize.append(node)
                    color_values.append(mix["a-ctxDDEL"].iloc[i])

        for y, n in enumerate(nodesize):
            if n>0 and n<15:
                nodesize[y] = 15
            if n == 0:
                nodesize[y] = 50

        cm = plt.get_cmap("Blues")
        vmin = min(color_values)
        vmax = max(color_values)
     #   net = networkx.draw_networkx(G, node_size=nodesize, node_color = color_values, cmap = cm,with_labels=False, edgecolors = "Black", ax = ax)
        networkx.draw_networkx_nodes(G, pos=networkx.spring_layout(G), node_size=nodesize, node_color = color_values, cmap = cm,edgecolors = "Black")
        ax.set_title(experiment)
     #   sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=vmin, vmax=vmax))
      #  cbar = plt.colorbar(sm)
        plot_number = plot_number+1

    antigens = ['Ecarpholin_S_DDEL_5ug/mL', 'Myotoxin_II_DDEL_5ug/mL',
           'PLA2_n.naja(Uniprot:P15445)_DDEL_5ug/mL', 'Ecarpholin_S_CDEL_100nM',
           'Myotoxin_II_CDEL_100nM',]
def cluster_antigens(sequencing_report, sample, antigens):
    plot_number = 0
    experiment = "Library_1_F5_2"
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
    batch = sequencing_report[sequencing_report["Experiment"] == sample]
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
        clone_counts = mix.cloneCount
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
        #   sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=vmin, vmax=vmax))
        #  cbar = plt.colorbar(sm)
        plot_number = plot_number + 1
    fig.suptitle("Clustering for sample " + sample)
#    sm = plt.cm.ScalarMappable(cmap=colormap,
 #                              norm=plt.Normalize(vmin = vmin, vmax=max_colorbar))
    #col_ax = fig.add_subplot(1, 4, 6)
    #fig.colorbar(sm, col_ax)
















example = sequencing_report_all[sequencing_report_all["Experiment"] == "Library_4_F7_2"]
example = example.iloc[:200, :]
example = example.aaSeqCDR3
x = np.array(example)
y = np.array(example)
x = list(x)
y = list(y)
matrix = np.zeros((len(x), len(y)))
for i in range(len(x)):
    for j in range(len(y)):
        dist = levenshtein_distance(x[i], y[j])
        matrix[i, j] = dist


my_pcoa = pcoa(matrix, 2)

from Bio.Align.substitution_matrices import blosum62 as blosum
def score_pairwise(seq1, seq2, matrix, gap_s, gap_e, gap = True):
    for A,B in zip(seq1, seq2):
        diag = ('-'==A) or ('-'==B)
        yield (gap_e if gap else gap_s) if diag else matrix[(A,B)]
        gap = diag

matrix = np.zeros((len(x), len(y)))
for i in range(len(x)):
    for j in range(len(y)):
        dist = sum(score_pairwise(x[i], y[j], blosum, -3, -1))
        matrix[i, j] = dist

from sklearn.preprocessing import StandardScaler
x = StandardScaler().fit_transform(matrix)
from sklearn.decomposition import PCA
pca = PCA(n_components=3)
principalComponents = pca.fit_transform(x)
principalDf = pd.DataFrame(data = principalComponents
             , columns = ['principal component 1', 'principal component 2', "principal component 3"])
finalDf = pd.concat([principalDf, pd.Series(y)], axis = 1)

# Principal components correlation coefficients
loadings = pca.components_
# Number of features before PCA
n_features = pca.n_features_
# Feature names before PCA
feature_names = y
# PC names
pc_list = [f'PC{i}' for i in list(range(1, n_features + 1))]
# Match PC names to loadings
pc_loadings = dict(zip(pc_list, loadings))
# Matrix of corr coefs between feature names and PCs
loadings_df = pd.DataFrame.from_dict(pc_loadings)
loadings_df['feature_names'] = feature_names
loadings_df = loadings_df.set_index('feature_names')
loadings_df
# Get the loadings of x and y axes
xs = loadings[0]
ys = loadings[1]
# Plot the loadings on a scatterplot
for i, varnames in enumerate(feature_names):
    plt.scatter(xs[i], ys[i], s=200)
    plt.text(xs[i], ys[i], varnames)
# Define the axes
xticks = np.linspace(-0.8, 0.8, num=5)
yticks = np.linspace(-0.8, 0.8, num=5)
plt.xticks(xticks)
plt.yticks(yticks)
plt.xlabel('PC1')
plt.ylabel('PC2')
# Show plot
plt.title('2D Loading plot')
plt.show()
plt.xlim(-0.1, 0.1)
plt.ylim(-0.1, 0.1)