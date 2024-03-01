from src.ExpoSeq.pipeline import PlotManager
from src.ExpoSeq.plots.tidy_protbert_embedding import TransformerBased
import umap
from sklearn.naive_bayes import GaussianNB
from sklearn.model_selection import train_test_split


plot = PlotManager(experiment = "max_both", 
                   length_threshold=5,
                   min_read_count=5,
                   no_automation=True,
                   allow_binding_data=False)
plot.change_region(region = "targetSequences")
sequencing_report = plot.sequencing_report.head(100)
assert "aaSeqtargetSequences" in sequencing_report.columns.to_list()
Transform = TransformerBased()
clone_fraction = sequencing_report["cloneFraction"]
sequences = sequencing_report["aaSeqtargetSequences"].to_list()
sequences_list = Transform.embedding_per_seq(sequences)
trans = umap.UMAP(n_neighbors=5, n_components=10, random_state=42).fit(sequences_list)
print(trans)

X_train, X_test, y_train, y_test = train_test_split(trans, y, test_size=0.5, random_state=0)
gnb = GaussianNB()
y_pred = gnb.fit(X_train, y_train).predict(X_test)
print("Number of mislabeled points out of a total %d points : %d"
       % (X_test.shape[0], (y_test != y_pred).sum()))