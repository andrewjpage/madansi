import networkx as nx
from madansi.ConnectedComponents import ConnectedComponents

class JoiningContigComponents(object):
    
    def __init__(self, refined_graph, unrefined_graph):
        self.unrefined_graph = unrefined_graph
        self.refined_graph = refined_graph
        self.components_dictionary = ConnectedComponents(self.refined_graph, self.unrefined_graph).create_component_dictionary()
        self.contigs_to_components = ConnectedComponents(self.refined_graph, self.unrefined_graph).create_contigs_to_components_dictionary()
    
    def list_ends_of_components(self):
        ends_of_components = []
        for component in self.components_dictionary:    
            for end in self.components_dictionary[component].ends:
                ends_of_components.append(end)
        return ends_of_components
        
    def add_edges(self):
        """Looks through the ends of each component of contigs and if the degree of this contig in the unrefined graph is at most two, \
        will add connections to adjacent contigs if the other contig is also at the end of a different component or is not marked as being \
        on a component"""
        ordered_contig_graph = nx.Graph()
        ordered_contig_graph.add_edges_from(self.refined_graph.edges())
        ends_of_components = self.list_ends_of_components()
        for end in ends_of_components:
            if self.unrefined_graph.degree(end)<=2:
                for edge in self.unrefined_graph.edges(end):
                    if edge[0] and edge[1] in ends_of_components:
                        ordered_contig_graph.add_edge(edge[0], edge[1])
                    elif self.contigs_to_components[edge[0]] == None or self.contigs_to_components[edge[1]] == None:
                        ordered_contig_graph.add_edge(edge[0], edge[1])
                        
        return ordered_contig_graph
        
    