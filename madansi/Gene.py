class Gene(object):
    """Sets up the gene object"""
    def __init__(self, orientation, start, end, node, contig):
        self.orientation = orientation
        self.start = start
        self.end = end
        self.node = node
        self.contig = contig
        
