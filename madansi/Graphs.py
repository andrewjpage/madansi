import networkx as nx
import os
from madansi.BlastHit import BlastHit

class Error(Exception): pass

class Graphs(object):
    """Reads in a graph in dot format and will compare the genes present to those given in the dictionary of genes present"""
    
    def __init__(self,graphfile,filteredfile):
        self.graphfile= graphfile
        self.filteredfile=filteredfile
    
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
        
    def open_graph(self):
        """Opens the graph dot file for reading"""
        g=nx.Graph(nx.drawing.nx_pydot.read_dot(self.graphfile))
        return g
        
    def extra_genes_in_graph(self):
        """Finds the additional genes noted in the graph not marked as present in the lookup table"""
        g=self.open_graph()
        gene_dict = self.construct_dictionary()
        extra_genes=[]
        for gene in nx.nodes_iter(g):
            if not gene_dict[gene]:
                extra_genes.append(gene)
        return extra_genes

    def genes_not_in_graph(self):
        """Finds the genes listed in the lookup table that are not found as nodes in the graph"""
        g=self.open_graph()
        gene_dict = self.construct_dictionary()
        gene_not_in_graph=[]
        for gene in gene_dict:
            if gene_dict[gene]:
                if gene not in nx.nodes_iter(g):
                    gene_not_in_graph.append(gene)
        return gene_not_in_graph