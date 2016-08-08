import networkx as nx
from madansi.Component import Component

class ConnectedComponents(object):
    
    def __init__(self, refined_graph):
        self.refined_graph = refined_graph
        self.components_dictionary = {}
    
    def connected_components(self): #At some point maybe think about checking that the components do not have cycles
        component_list = sorted(nx.connected_components(self.refined_graph), key = len, reverse = True)
        for component in component_list:
            edge_list = sorted([tuple(sorted(edge)) for edge in self.refined_graph.edges(component)])
            if len(component)==1:
                ends = component
            else:
                ends = [contig for contig in component if self.refined_graph.degree(contig)==1]
            self.components_dictionary[component_list.index(component)] = Component(component_list.index(component), component,\
                                                                                 ends, edge_list)
        return self.components_dictionary