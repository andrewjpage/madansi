import networkx as nx
from madansi.Component import Component

class ConnectedComponents(object):
    
    def __init__(self, refined_graph, unrefined_graph):
        self.refined_graph = refined_graph
        self.unrefined_graph = unrefined_graph
        self.component_list = self.connected_component_list()
    
    def connected_component_list(self):
        component_list = sorted(nx.connected_components(self.refined_graph), key = len, reverse = True)
        return component_list
    
    def create_component_dictionary(self): #At some point maybe think about checking that the components do not have cycles
        """Creates a dictionary including an index for each connected component in the refined graph and the component object"""
        components_dictionary = {}
        for component in self.component_list:
            edge_list = sorted([tuple(sorted(edge)) for edge in self.refined_graph.edges(component)])
            if len(component)==1:
                ends = component
            else:
                ends = [contig for contig in component if self.refined_graph.degree(contig)==1]
            components_dictionary[self.component_list.index(component)] = Component(self.component_list.index(component), component,\
                                                                                 ends, edge_list)
        return components_dictionary
    
    def create_contigs_to_components_dictionary(self):
        """Creates a dictionary that relates the contig to the component that it is on if it lies on the refined graph and gives None if the contig is not 
        present in the refined graph"""
        contigs_to_components = {}
        for component in self.component_list:
            component_index = self.component_list.index(component)
            for contig in component:
                contigs_to_components[contig] = component_index
        for contig in self.unrefined_graph.nodes():
            if contig not in contigs_to_components:
                contigs_to_components[contig] = None        
        return contigs_to_components
                