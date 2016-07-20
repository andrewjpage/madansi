import networkx as nx
from random import choice

class Error(Exception): pass
class WalkPaths(object):
"""Given a node on the given graph, will check that the gene that it represents is present from the lookup table and considers its neighbours to see if they are also present."""    
    def __init__(self,graphfile,filteredfile):
        self.graphfile = graphfile
    
    def open_graphfile(self):
        """Open the given graph file for walking"""
        try:
            g=nx.Graph(nx.drawing.nx_pydot.read_dot(self.graphfile))
            return g
        except IOError:
            raise Error("Error opening this file")
    
    def construct_dictionary(self):
        gene_present_dict = {}
        
        try:
            f = open(self.filteredfile)
        except IOError:
            raise Error("Error opening this file")
        
        for line in f:
            bh = BlastHit(line)
            if bh.bit_score >= 150:
                gene_present_dict[bh.qry_name] = True
            else:
                gene_present_dict[bh.qry_name] = False
        
        f.close() 
        return(gene_present_dict)
        
    def add_node_attribute(self):
        """Adds node attribute to the graph based on whether the gene is given as present in the lookup table"""    
        g = self.open_graphfile()
        gene_dict = self.construct_dictionary()
        for gene in nx.nodes_iter():
            if gene_dict[gene]:
                g.node[gene]['present']=True
            else:
                g.node[gene]['present']=False
    
    def choose_node(self):
        g = self.open_graphfile()
        chosen_node = choice(g.nodes())
        while g.node[chosen_node]:
            
        
    def find_neighbours(self):
        g= self.open_graphfile()
        gene_present_dict= self.construct_dictionary()
        if not gene_present_dict[first_node]:
            g.node[gene]['mark']=True
            
    def remove_added_edge_node_attributes(self):
        g = self.open_graphfile()
        for gene in nx.nodes_iter():
            del g.node[gene]['present']
            del g.node[gene]['mark']
        
            
            