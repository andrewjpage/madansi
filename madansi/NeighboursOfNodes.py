import networkx as nx

class NeighboursOfNodes(object): 
    """Given a list of nodes in the graph it will find the nighbours of this set"""
    def __init__(self,graph):
        self.graph = graph
        self.seen_nodes = []

    def find_neighbours(self, node_list):
    
        for node in node_list:
            for neighbour in self.graph.neighbors(node):
                if neighbour not in node_list and neighbour not in self.seen_nodes:
                    self.seen_nodes.append(neighbour)
        return self.seen_nodes