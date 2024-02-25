
class FullSequence:
    def __init__(self, avail_regions):
        pass
        self.avail_regions = avail_regions
        self.region_order = ["FR1", "CDR1", "FR2", "CDR2", "FR3", "CDR3", "FR4"]
        
    @staticmethod
    def find_connecting_seq(avail_regions):
        region_order = ["FR1", "CDR1", "FR2", "CDR2", "FR3", "CDR3", "FR4"]
        region_matrix = pd.DataFrame(0, index=region_order, columns=region_order)
        len_matrix = len(region_order)
        for region_index in range(len_matrix): # there is a more elegant way if you use the avail regions as array but this was quicker
            consecutive = 0
            for region_row in range(len_matrix):
                if region_matrix.index[region_row] in avail_regions:
                    region_matrix.iloc[region_row, region_index] = 1
                    consecutive +=1
                    max_consecutive_col = consecutive
                else:
                    region_matrix.iloc[region_row, region_index] = 0
                    consecutive = 0

