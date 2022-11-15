import matplotlib.pyplot as plt
from numpy import unique
from python_scripts.tidy_data.clustering import cleaning
from networkx import draw_networkx

G, degree_sequence = cleaning(sample_name = sample_name)
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