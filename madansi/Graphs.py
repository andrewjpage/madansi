import networkx as nx
import os
from madansi.GenePresent import GenePresent

class Error(Exception): pass

class Graphs(object):
    """Reads in a graph in dot format and will compare the genes present to those given in the dictionary of genes present"""
    
    def __init__(self,graphfile,filteredfile):
        self.graphfile= graphfile
        self.filteredfile=filteredfile
        
    def open_graph(self):
        g=nx.Graph(nx.drawing.nx_pydot.read_dot(self.graphfile))
        
    def comparison_graph_lookup_table(self):
        gp=GenePresent(self.filteredfile)
        gp_dict=gp.construct_dictionary()
        G=Graphs.open_graph(self)
        X= [v for k,v in gp_dict if v]
        print(X)
        if list(G.nodes())==[v for k,v in gp_dict if v]:
            return True
        else:
            return False