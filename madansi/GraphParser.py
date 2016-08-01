import networkx as nx
class GraphParser(object):
    
    def __init__(self, graphfile):
        self.graphfile = graphfile
        self.graph = self.filter_graph(self.open_graph_file())
        
    def filter_graph(self, graph):
        return graph

    def open_graph_file(self):
        try:
            return nx.Graph(nx.drawing.nx_pydot.read_dot(self.graphfile))
        except IOError:
            raise IOError("Error opening this file")