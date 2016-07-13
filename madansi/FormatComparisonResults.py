from Bio.Blast import NCBIXML

class FormatComparisonResults(object):
    """Parses the output from a BLAST comparison and extracts the information"""
    def __init__(self, comparisonfile):
        self.comparisonfile = comparisonfile
        
    def parse_records(self):
    