import networkx as nx

class ContigGraph(object):
    """Given a list of tuples of pairs of contigs together with the distance between them, will generate a weighted graph showing all\
    of these connections together with the distance between the two contigs."""
    
    def __init__(self, neighbouring_contigs):
        self.neighbouring_contigs = neighbouring_contigs
        self.contig_graph = nx.Graph()
    
    def create_contig_subgraph(self):
        for entry in self.neighbouring_contigs:
            self.contig_graph.add_nodes_from([entry[0][0], entry[0][1]])
            self.contig_graph.add_edge(entry[0][0],entry[0][1], weight = entry[1])
        return self.contig_graph
    
        