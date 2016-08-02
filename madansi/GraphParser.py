import networkx as nx
class GraphParser(object):
    
    def __init__(self, graph_file):
        self.graph_file = graph_file
        self.graph = self.filter_graph(self.open_graph_file())
        
    def filter_graph(self, graph):
        nodes_to_remove = []
        for node in graph.nodes_iter():
            if graph.degree(node) > 3:
                nodes_to_remove.append(node)
        graph.remove_nodes_from(nodes_to_remove)
        return graph

    def open_graph_file(self):
        try:
            return nx.Graph(nx.drawing.nx_pydot.read_dot(self.graph_file))
        except IOError:
            raise IOError("Error opening this file")