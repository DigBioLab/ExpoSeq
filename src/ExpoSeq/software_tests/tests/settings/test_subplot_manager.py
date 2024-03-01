from src.ExpoSeq.pipeline import PlotManager
from src.ExpoSeq.settings.subplots_manager import Subplotter


def test_subplotter():
    plot = PlotManager(experiment = "test_show", module_dir = "src/ExpoSeq/software_tests/test_files", allow_binding_data=False, test_version=True,)
    plot.lengthDistribution_single()
    Sub = Subplotter( "abc")
    Sub.add_as_subplot(plot.ControlFigure.fig)
    plot.logoPlot_single()
    Sub.add_as_subplot(plot.ControlFigure.fig)
    plot.jaccard()
    Sub.add_as_subplot(plot.ControlFigure.fig)
    plot.aa_distribution()
    Sub.add_as_subplot(plot.ControlFigure.fig)

    assert len(Sub.files) == 4
    Sub.make_figure()
