import networkx as nx


class SimplifyingGraph(object):
    
    def __init__(self, contig_graph, output_graph):
        self.contig_graph = contig_graph
        self.output_graph = output_graph
        
    def simplifying_graph(self):    
        changes_made = True
        while changes_made == True:
            for contig in self.contig_graph.nodes():
                changes_made = False
                contig_edge_dict = {}
                if self.contig_graph.degree(contig)>2:
                    for edge in self.contig_graph.edges(contig):
                        contig_edge_dict[edge] = int(self.contig_graph.edge[edge[0]][edge[1]]['weight'])
                    max_weight = max([value for key, value in contig_edge_dict.items()])
                    for key in contig_edge_dict:
                        if contig_edge_dict[key] == max_weight:
                            if self.contig_graph.degree(key[0]) > 1 and self.contig_graph.degree(key[1]) > 1:
                                self.contig_graph.remove_edge(key[0], key[1])
                                changes_made = True
                                break
        return self.contig_graph
    
    def output_simplified_graph(self):
        self.simplifyinging_graph()
        nx.drawing.nx_pydot.write_dot(self.contig_graph, self.output_graph)
            
                         