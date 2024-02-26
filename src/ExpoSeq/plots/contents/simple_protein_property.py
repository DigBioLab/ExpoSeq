import peptides
from . import hytdrophobicity_metrics
from . import isoelectric_point_scales

class GetProteinProperty:

    def __init__(self, sequences):
        """This class allows you to return a dictionary with the sequences and the sequence property of interest

        Args:
            sequences (list): list of peptide sequences
        """
        self.sequences = sequences
        self.sequence_property_interest = {}
        self.attribute_funcs = {"aliphatic_index": self.get_aliphatic_index,
        "isoelectric_point": self.get_isoelectric_point,
        "hydrophobicity": self.get_hydrophobicity,
        "weight": self.get_weight,
        "mass_charge_ratio": self.mass_charge_ratio,
        "length": self.calc_len,
     }
        
    @staticmethod
    def calc_len(peptide_object):
        return len(peptide_object.sequence)
    
        
    @staticmethod
    def create_peptide_object(single_sequence):
        peptide = peptides.Peptide(single_sequence)
        return peptide
    
    @staticmethod
    def get_aliphatic_index(peptide_object):
        return peptide_object.aliphatic_index()
    
    @staticmethod
    def get_isoelectric_point(peptide_object, scale = "EMBOSS"):
        assert scale in isoelectric_point_scales.scale, f"Please choose on of the metrics from: {isoelectric_point_scales.scale}"
        return peptide_object.isoelectric_point(pKscale = scale)
    
    @staticmethod
    def get_hydrophobicity(peptide_object, scale = "KyteDoolittle"):
        assert scale in hytdrophobicity_metrics.metrics, f"Please choose on of the metrics from: {hytdrophobicity_metrics.metrics}"
        return peptide_object.hydrophobicity(scale = scale)
    
    @staticmethod
    def get_weight(peptide_object):
        return peptide_object.molecular_weight()
    
    @staticmethod
    def mass_charge_ratio(peptide_object):
        return peptide_object.mz()
    
    def calc_attribute(self, attribute, **kwargs):
        assert attribute in self.attribute_funcs.keys(), f"Please enter a value for attribute which is in {self.attribute_funcs.keys()}"
        for sequence in self.sequences:
            PepObj = self.create_peptide_object(single_sequence=sequence)
            func = self.attribute_funcs.get(attribute)
            attribute_value = func(PepObj, **kwargs)
            self.sequence_property_interest[sequence] = attribute_value
            
    def add_attributes_to_table(self, table, attribute):
        table[attribute] = list(self.sequence_property_interest.values())
        return table
      
    
