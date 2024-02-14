import pandas as pd
from src.ExpoSeq.plots.matrices.morosita_horn_matrix import PrepareData
from src.ExpoSeq.plots.matrices.make_matrix import IdentityMatrix
import numpy as np
import matplotlib.pyplot as plt


def test_matrix():
    sequencing_report_path = r"src/ExpoSeq/software_tests/test_files/test_show/sequencing_report.csv"
    sequencing_report = pd.read_csv(sequencing_report_path)
    sequencing_report["cloneFraction"] = sequencing_report["readFraction"]
    unique_sequences, unique_experiments = PrepareData().prepare_unique_seq_table(sequencing_report, "aaSeqCDR3")
    assert type(unique_experiments) == list,"unique_experiments is not a list"
    assert isinstance(unique_sequences, pd.DataFrame), "unique sequences is not a pandas dataframe"

    matrix = PrepareData.calc_matrix_values(unique_sequences, unique_experiments)
    assert matrix.max().max() == 1, "Maximum value of Matrix is not 1"
    numpy_matrix = matrix.to_numpy()
    diagonal_elements = np.diag(numpy_matrix)
    assert np.all(diagonal_elements == 1), "Not all diagonal elements are equal to 1."
    assert matrix.shape[1] == matrix.shape[0], "matrix is not symmetrical"
    assert matrix.shape[0] == len(unique_experiments), "matrix "
    IdentityMatrix(sequencing_report, "aaSeqCDR3", "morosita_horn", {})
    fig = plt.figure(1, figsize = (12, 10))
    ax = fig.gca()
    IdentityMatrix(sequencing_report, "aaSeqCDR3","sorensen",{})
    fig = plt.figure(1, figsize = (12, 10))
    ax = fig.gca()
    font_settings = {'fontfamily': 'serif', 'fontsize': '18', 'fontstyle': 'normal', 'fontweight': 'bold'}
    colorbar_settings = {'cmap': 'inferno', 'orientation': 'vertical', 'spacing': 'proportional', 'extend': 'neither'}
    IdentityMatrix(sequencing_report, "aaSeqCDR3","jaccard", colorbar_settings=colorbar_settings)
    fig = plt.figure(1, figsize = (12, 10))
    ax = fig.gca()
    IdentityMatrix(sequencing_report, "aaSeqCDR3","relative", {})

    
    