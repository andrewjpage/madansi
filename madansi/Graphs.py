import networkx as nx
import os

class Error(Exception): pass

class Graphs(object):
    """Reads in a graph in dot format and will compare the genes present to those given in the dictionary of genes present"""
    
    def __init__(self,graphfile):
        self.graphfile= graphfile
        
    def open_graph(self):
        g=nx.Graph(nx.drawing.nx_pydot.read_dot(self.graphfile))
        
    def comparison_graph_lookuptable(self):