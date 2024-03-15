from matplotlib_venn import venn2, venn3
from .matrices import morosita_horn_matrix
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

class PrepareData:
    @staticmethod
    def check_input(sequencing_report, samples):
        for sample in samples:
            assert sample in sequencing_report["Experiment"].to_list(), f"Please provide a valid sample name. The samples are: {sequencing_report['Experiment'].to_list()}"
        assert len(samples) == 2 or len(samples) == 3, "Please provide 2 or 3 samples"
        
    def tidy(self, sequencing_report:pd.DataFrame, samples:list, region_of_interest:str, ):
        self.check_input(sequencing_report, samples)
        PrepData = morosita_horn_matrix.PrepareData()
        matrix, unique_experiments = PrepData.cleaningPlot(samples, sequencing_report, region_string=region_of_interest, protein = True)
        readC_sample = sequencing_report.groupby("Experiment")["readCount"].sum() # the problem of this is that the an 80 % identity is not optimally visualized because a smaller circle (because of lower read count) can have issuesto fill 80 % of the bigger circle. 
        readC_sample = readC_sample.iloc[:len(samples)].values
        morosita_horn_base = matrix.loc[samples, samples].values
        morosita_horn_readC  = morosita_horn_base  # multiply identity with overall read count 
        morosita_horn_readC = np.around(morosita_horn_readC, 2)

        if morosita_horn_readC.shape[0] == 2:
            # Create a Venn diagram for 2 sets
            subsets=(morosita_horn_readC[0, 0], morosita_horn_readC[1, 1], morosita_horn_readC[1, 0])
        else:
            # Create a Venn diagram for 3 sets
            
            subsets=(morosita_horn_readC[0, 0], morosita_horn_readC[1, 1], morosita_horn_readC[1, 0], morosita_horn_readC[2, 2], morosita_horn_readC[2,0], morosita_horn_readC[2, 1], morosita_horn_readC[2, 1] * morosita_horn_readC[0, 1])
        
        return subsets, morosita_horn_base
    
class VennDiagram:
    def __init__(self, sequencing_report, samples, region_of_interest, font_settings = None, ax = None) :
        PrepData = PrepareData()
        subsets, morosita_horn_base = PrepData.tidy(sequencing_report, samples, region_of_interest)
        self.ax = ax
        if ax != None:
            self.create_venn(subsets, samples)

    def create_venn(self,subsets, samples):
            # Create a Venn diagram for 3 sets
        if len(samples) == 3:
            assert len(subsets) == 7, f"Yuu have {len(subsets)}"
            venn3(subsets, ax= self.ax, set_labels = samples)
        elif len(samples) == 2:
            venn2(subsets, ax= self.ax, set_labels = samples,)
            

            
    
    
morosita_horn = pd.read_csv(r"C:\Users\nilsh\my_projects\ExpoSeq\tmp_test\morosita_horn_test.csv")
sequencing_report = pd.read_csv(r"C:\Users\nilsh\my_projects\ExpoSeq\src\ExpoSeq\software_tests\test_files\test_show\sequencing_report.csv")
sequencing_report["cloneFraction"] = sequencing_report["readFraction"]
fig = plt.figure(1)
ax = fig.gca()
VennDiagram(sequencing_report, samples = ["GeneMind_1", "GeneMind_5", "GeneMind_6"], region_of_interest = "aaSeqCDR3", ax =ax)

plt.show()

