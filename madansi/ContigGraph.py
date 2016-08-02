import networkx as nx

class ContigGraph(object):
    """Given a list of tuples of pairs of contigs together with the distance between them, will generate a weighted graph showing all\
    of these connections together with the distance between the two contigs."""
    
    def __init__(self, neighbouring_contigs, output_contig_graph):
        self.neighbouring_contigs = neighbouring_contigs
        self.contig_graph = nx.Graph()
        self.output_contig_graph = output_contig_graph
    
    def create_contig_subgraph(self):
        for entry in self.neighbouring_contigs:
            self.contig_graph.add_nodes_from([entry[0], entry[1]])
            self.contig_graph.add_edge(entry[0],entry[1], weight = entry[2])
        nx.drawing.nx_pydot.write_dot(self.contig_graph, self.output_contig_graph) #Should put this into a second method
        return self.contig_graph
    
    def output_contig_graph(self): 
        nx.drawing.nx_pydot.write_dot(self.contig_graph, self.output_contig_graph)
        