import networkx as nx
import pprint
from madansi.GenesToContig import GenesToContig
from madansi.NeighboursOfNodes import NeighboursOfNodes
import sys

class RefineContigNeighbours(object):
    
    def __init__(self, neighbouring_contigs, filtered_graph, filtered_blast_file, gene_detector):
        self.neighbouring_contigs   = neighbouring_contigs
        self.filtered_graph         = filtered_graph
        self.filtered_blast_file    = filtered_blast_file
        self.genes                  = GenesToContig(self.filtered_blast_file).genes_to_contig()
        self.refined_neighbouring_contigs = []
        self.gene_detector          = gene_detector
        self.contigs                = gene_detector.contigs_to_genes()
        self.components_of_contigs  = {}
        self.contig_ends            = {}
        self.contig_orientation     = {}
    
    def add_to_contig_appearance(self, gene, contig_appearances, iteration):
        if gene in self.genes:
            if self.genes[gene] not in contig_appearances:
                contig_appearances[self.genes[gene]] = [1, {iteration:[gene]}]
            else:
                contig_appearances[self.genes[gene]][0] += 1
                if iteration in contig_appearances[self.genes[gene]][1]:
                    contig_appearances[self.genes[gene]][1][iteration].append(gene)
                    contig_appearances[self.genes[gene]][1][iteration].sort()
                else:
                    contig_appearances[self.genes[gene]][1][iteration] = [gene]
        return contig_appearances
    
    #def count_contig_appearances(self, gene_list, contig_appearances):
    #    for gene in gene_list:
    #        contig_appearances = self.add_to_contig_appearance(gene, contig_appearances)
    #    return contig_appearances
    
    def find_contig_appearances(self, neighbours):
        seen_nodes = []
        contig_appearances = {neighbours[0][0] : [0, {}], neighbours[0][1] : [0,{}]}
        
        for intersection in neighbours[2]:
            seen_nodes.append(intersection)
            self.add_to_contig_appearance(intersection, contig_appearances,0)
            
        for i in range(neighbours[1] + 2):
            neighbouring_nodes = NeighboursOfNodes(self.filtered_graph).find_neighbours(seen_nodes)
            
            for neighbour_node in neighbouring_nodes:
                seen_nodes.append(neighbour_node)
                self.add_to_contig_appearance(neighbour_node, contig_appearances, i + 1)
        
        return contig_appearances
    
    
    def setup_contig_max_iteration(self):
        """Initialises the maximum contig weight and sets each contig as being on a separate component"""
        contig_count = 0
        contig_max_iteration = {}
        for contig in self.contigs:
            contig_max_iteration[contig] = [contig_count, 0]
            self.components_of_contigs[contig_count] = [contig]
            contig_count += 1
        return contig_max_iteration
    
    def check_if_connection_is_valid(self, neighbours):
        """Only keeps the cases where there is not a gene by itself or too many genes, coming from an intersection towards the middle of the contig"""
        contig_appearances = self.find_contig_appearances(neighbours)
        if  contig_appearances[neighbours[0][0]][0] in range(min(2, len(self.gene_detector.contigs_to_genes()[neighbours[0][0]].gene_objects)), neighbours[1] + 4) and\
                         contig_appearances[neighbours[0][1]][0] in range(min(2, len(self.gene_detector.contigs_to_genes()[neighbours[0][1]].gene_objects)), neighbours[1] + 4):
            if len(contig_appearances) == 2:
                return neighbours

            elif len(contig_appearances) == 3:
                if len(self.gene_detector.contigs_to_genes()[neighbours[0][0]].gene_objects) in range(1, neighbours[1] + 3):
                    if contig_appearances[neighbours[0][0]][0] == len(self.gene_detector.contigs_to_genes()[neighbours[0][0]].gene_objects):
                        return neighbours
                elif len(self.gene_detector.contigs_to_genes()[neighbours[0][1]].gene_objects) in range(1, neighbours[1] + 3):
                    if contig_appearances[neighbours[0][1]][0] == len(self.gene_detector.contigs_to_genes()[neighbours[0][1]].gene_objects):
                        return neighbours
            
    def refine_contig_neighbours(self):
        """Removes the case when we have a short contig between two others."""
        contig_max_iteration = self.setup_contig_max_iteration()
        for neighbours in sorted(self.neighbouring_contigs, key= lambda x: x[1]):
            if self.check_if_connection_is_valid(neighbours) != None:
                neighbours = self.check_if_connection_is_valid(neighbours)
                min_component_number = min(contig_max_iteration[neighbours[0][0]][0],contig_max_iteration[neighbours[0][1]][0])
                max_component_number = max(contig_max_iteration[neighbours[0][0]][0],contig_max_iteration[neighbours[0][1]][0])
                if  min_component_number == max_component_number:
                    if neighbours[1] <= contig_max_iteration[neighbours[0][0]][1] and neighbours[1] <= contig_max_iteration[neighbours[0][1]][1]:
                        self.refined_neighbouring_contigs.append(neighbours)
                        contigs_to_add = [contig for contig in [neighbours[0][0], neighbours[0][1]] if contig not in self.components_of_contigs[min_component_number]]
                        for contig in contigs_to_add:
                            self.components_of_contigs[min_component_number].append(contig)
                else:
                    self.refined_neighbouring_contigs.append(neighbours)
                    if contig_max_iteration[neighbours[0][0]][0] == min_component_number:
                        contig_max_iteration[neighbours[0][0]] = [min_component_number, neighbours[1]]
                        contig_max_iteration[neighbours[0][1]] = [min_component_number, neighbours[1]]
                        contigs_to_add = [contig for contig in self.components_of_contigs[max_component_number]]
                        for contig in contigs_to_add:
                            self.components_of_contigs[min_component_number].append(contig)
                            contig_max_iteration[contig][0] = min_component_number
                        self.components_of_contigs.pop(max_component_number)
                    else:
                        contig_max_iteration[neighbours[0][1]] = [min_component_number, neighbours[1]]
                        contig_max_iteration[neighbours[0][0]] = [min_component_number, neighbours[1]]
                        contigs_to_add = [contig for contig in self.components_of_contigs[max_component_number]]
                        for contig in contigs_to_add:
                            self.components_of_contigs[min_component_number].append(contig)
                            contig_max_iteration[contig][0] = min_component_number
                        self.components_of_contigs.pop(max_component_number)
        return self.refined_neighbouring_contigs
    
    
    def ends_of_contigs(self):
        for neighbours in self.neighbouring_contigs:
            contig_appearances = self.find_contig_appearances(neighbours)
            for contig in [neighbours[0][0], neighbours[0][1]]:
                if contig not in self.contig_ends:
                    self.contig_ends[contig] = {}
            
            if self.check_for_one_gene(neighbours, contig_appearances, 0) != False:
                self.contig_ends[neighbours[0][0]][neighbours[0][1]] = self.check_for_one_gene(neighbours, contig_appearances, 0)
            if self.check_for_one_gene(neighbours, contig_appearances, 1) != False:
                self.contig_ends[neighbours[0][1]][neighbours[0][0]] = self.check_for_one_gene(neighbours, contig_appearances, 1)
            
    
        return self.contig_ends
    
    def check_for_one_gene(self, neighbours, contig_appearances, int):
        count = 0
        qry_starts = []
        for i in sorted(contig_appearances[neighbours[0][int]][1].keys()):
            if len(contig_appearances[neighbours[0][int]][1][i]) == 1:
                qry_starts.append(self.gene_detector.contigs_to_genes()[neighbours[0][int]].gene_objects[contig_appearances[neighbours[0][int]][1][i][0]].qry_start)
                count = count + 1   
            else:
                genes_degree_at_least_two = self.remove_gene_degree_one(neighbours, contig_appearances, int, i)
                if len(genes_degree_at_least_two) == 1:
                    qry_starts.append(self.gene_detector.contigs_to_genes()[neighbours[0][int]].gene_objects[genes_degree_at_least_two[0]].qry_start)
                    count = count + 1 
            if count == 2:
                return qry_starts
                
        return False
    
    def remove_gene_degree_one(self, neighbours, contig_appearances, int, i):
        genes_degree_at_least_2 = []
        for gene in contig_appearances[neighbours[0][int]][1][i]:
            if self.filtered_graph.degree(gene) != 1:
                genes_degree_at_least_2.append(gene)
        return genes_degree_at_least_2
                        
    
        
        
        
        