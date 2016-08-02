import networkx as nx

class ContigSpanningTree(object):
    def __init__(self,contig_graph, output_graph):
        self.contig_graph = contig_graph
        self.output_graph = output_graph
    
    def construct_spanning_tree(self):
        contig_spanning_tree = nx.minimum_spanning_tree(self.contig_graph)
        nx.drawing.nx_pydot.write_dot(contig_spanning_tree, self.output_graph)