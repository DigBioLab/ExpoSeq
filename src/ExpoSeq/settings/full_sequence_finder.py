import pandas as pd




class FullSequence:
    def __init__(self, avail_regions):
        for i in avail_regions:
            assert i in self.region_order
        self.avail_regions = avail_regions
        self.region_order = ["FR1", "CDR1", "FR2", "CDR2", "FR3", "CDR3", "FR4"]
        
    
    @staticmethod
    def find_most_ones(matrix_df):
        column_sums = matrix_df.sum(axis=0)  # Sum along columns
        most_connected_column = column_sums.idxmax()  # Index of the column with the highest sum
        n = 0
        z = 0
        col_list = matrix_df[most_connected_column].to_list()
        col_series = matrix_df[most_connected_column]
        for i in range( len(col_list)):
            value = col_list[i]
            if value ==1:
                n += 1
            else:
                n = 0
            fragment = matrix_df.index[i]
            if n > z or ( fragment == "CDR3" and n == z):
                end_ind = i 
                z = n
        end_ind += 1
        start = end_ind - z 
        assert start >=0
        row_indexes = col_series.index[start:end_ind].to_list()
        return row_indexes
    
    def find_connecting_seq(self):
        region_order = ["FR1", "CDR1", "FR2", "CDR2", "FR3", "CDR3", "FR4"]
        region_matrix = pd.DataFrame(0, index=region_order, columns=region_order)
        len_matrix = len(region_order)
        for region_index in range(len_matrix): # there is a more elegant way if you use the avail regions as array but this was quicker
            consecutive = 0
            for region_row in range(len_matrix):
                index_matrix = region_matrix.index[region_row]
                if index_matrix in self.avail_regions:
                    region_matrix.iloc[region_row, region_index] = 1
                    consecutive +=1
                    max_consecutive_col = consecutive
                else:
                    region_matrix.iloc[region_row, region_index] = 0
                    consecutive = 0
        row_indexes = self.find_most_ones(region_matrix)
        return row_indexes


