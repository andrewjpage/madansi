import networkx as nx

class CircularSubgraph(object):
    def __init__(self,contig_graph, output_graph, unused_contig_file):
        self.contig_graph = contig_graph
        self.output_graph = output_graph
        self.contig_subgraph = self.circularising_graph()   
        self.unused_contig_file = unused_contig_file  
    
    def circularising_graph(self):
        directed_contig_graph = self.contig_graph.to_directed()
        
        longest_cycle = []
        for cycle in nx.simple_cycles(directed_contig_graph):
            if len(cycle) > len(longest_cycle):
                longest_cycle = cycle
        
        contig_subgraph = nx.Graph()
        for i in range(len(longest_cycle) - 1):
            edge_weight = self.contig_graph.edge[longest_cycle[i]][longest_cycle[i+1]]['weight']
            contig_subgraph.add_edge(longest_cycle[i], longest_cycle[i+1], weight= edge_weight)
        
        edge_weight = self.contig_graph.edge[longest_cycle[len(longest_cycle) - 1]][longest_cycle[0]]['weight']
        contig_subgraph.add_edge(longest_cycle[len(longest_cycle) - 1], longest_cycle[0], weight = edge_weight)
        
        return contig_subgraph
    
    def output_circular_subgraph(self):
        nx.drawing.nx_pydot.write_dot(self.circularising_graph(), self.output_graph)
    
    def updating_unused_contig_file(self):
        with f = open(self.unused_contig_file, 'a'):
            for contig in [node for node in self.contig_subgraph.nodes() and not in self.contig_graph.nodes()]:
                f.write(contig)
        
        f.close()
        