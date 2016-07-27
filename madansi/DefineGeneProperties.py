from madansi.GenePresent import GenePresent

class DefineGeneProperties(self):
    
    def __init__(self, graphfile):
        self.graphfile = graphfile
        
    def open_graph_file(self):
        """Open the given graph file for searching"""
        try:
            g=nx.Graph(nx.drawing.nx_pydot.read_dot(self.graphfile))
            return g
        except IOError:
            raise IOError("Error opening this file")
    
            
    g = self.open_graph_file
    for node in g.nodes():
        def node