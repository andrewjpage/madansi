import networkx as nx
from madansi.GenesToContig import GenesToContig
from madansi.NeighboursOfNodes import NeighboursOfNodes
import sys

class RefineContigNeighbours(object):
    
    def __init__(self, neighbouring_contigs, filtered_graph, filtered_blast_file):
        self.neighbouring_contigs = neighbouring_contigs
        self.filtered_graph = filtered_graph
        self.filtered_blast_file = filtered_blast_file
        self.genes = GenesToContig(self.filtered_blast_file).genes_to_contig()
        self.refined_neighbouring_contigs = []
    
    def refine_contig_neighbours(self):
        for neighbours in self.neighbouring_contigs:
            seen_nodes = []
            for neighbour in neighbours[2]:
                seen_nodes.append(neighbour)
            for i in range(neighbours[1] + 2):
                neighbouring_nodes = NeighboursOfNodes(self.filtered_graph).find_neighbours(seen_nodes)
                for neighbour_node in neighbouring_nodes:
                    seen_nodes.append(neighbour_node)
            contig_appearances = self.count_contig_appearance(seen_nodes)
            if contig_appearances[neighbours[0][0]] > 1 and contig_appearances[neighbours[0][1]] > 1 and len(contig_appearances)==2:
                    self.refined_neighbouring_contigs.append(neighbours)          
        return self.refined_neighbouring_contigs
    
    def count_contig_appearance(self, gene_list):
        contig_appearances = {}
        for gene in gene_list:
            if gene in self.genes:
                if self.genes[gene] not in contig_appearances:
                    contig_appearances[self.genes[gene]] = 1                
                elif self.genes[gene] in contig_appearances:
                    contig_appearances[self.genes[gene]] += 1
        return contig_appearances