class Gene(object):
    """Sets up the gene object"""
    def __init__(self, orientation, start, end, node, contig, qry_start):
        self.orientation = orientation
        self.start = start
        self.end = end
        self.node = node
        self.contig = contig
        self.qry_start = qry_start
        
