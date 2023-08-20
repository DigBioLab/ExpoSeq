from .design.anno_sequence import SequenceManipulator


class DesignAB:
    def __init__(self, pipeline: PlotMananger):
        self.Pipeline = pipeline
        
    def manipulate_sequence(self, seq):
        #seq = "QMQLVQSGAEVKKPGSSVKVSCKASGGTFSSYAISWVRQAPGQGLEWMGGIIPIFGTANYAQKFQGRVTITADESTSTAYMELSSLRSEDTAVYYCARSSGYYGMDVWGQGTLVTVSS"  # Replace this with your sequence
        app = SequenceManipulator(seq)
        app.run()
        my_sequence = app.sequence
        return my_sequence

