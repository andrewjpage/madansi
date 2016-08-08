import networkx as nx
class GraphParser(object):
    
    def __init__(self, graph_file):
        self.graph_file = graph_file
        self.graph = self.filter_graph(self.open_graph_file())
       
        
    def filter_graph(self, graph):
        edges_to_remove = []
        for node in graph.nodes_iter():
            if graph.degree(node) >= 5:
                edge_list = graph.edges(node, data=True)
                max_weight = max([d['weight'] for u,v,d in edge_list])
                for edge in edge_list:
                    if edge[2]['weight'] == max_weight:
                        edges_to_remove.append(edge)
        graph.remove_edges_from(edges_to_remove)
        return graph

    def open_graph_file(self):
        try:
            return nx.Graph(nx.drawing.nx_pydot.read_dot(self.graph_file))
        except IOError:
            raise IOError("Error opening this file")