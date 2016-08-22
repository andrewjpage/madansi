import networkx as nx
import pprint
from madansi.GenesToContig import GenesToContig
from madansi.NeighboursOfNodes import NeighboursOfNodes
import sys

class RefineContigNeighbours(object):
    
    def __init__(self, neighbouring_contigs, filtered_graph, filtered_blast_file, gene_detector, sequences):
        self.neighbouring_contigs           = neighbouring_contigs
        self.filtered_graph                 = filtered_graph
        self.filtered_blast_file            = filtered_blast_file
        self.genes                          = GenesToContig(self.filtered_blast_file).genes_to_contig()
        self.refined_neighbouring_contigs   = []
        self.gene_detector                  = gene_detector
        self.contigs                        = gene_detector.contigs_to_genes()
        self.components_of_contigs          = {}
        self.contig_ends                    = {}
        self.contig_orientation             = {}
        self.sequences                      = sequences
    
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
        
    def most_occurent_contig(self, contig_appearances, neighbours):
        """Checks to see if within a """
        max_appearances = 0
        most_occurent_contig = None
        for k,v in contig_appearances.items():
            if k in [neighbours[0][0], neighbours[0][1]] and v[0] >= max_appearances or k not in [neighbours[0][0], neighbours[0][1]] and v[0] > max_appearances:
                max_appearances = v[0]
                most_occurent_contig = k
        return most_occurent_contig
        
    def check_for_small_contig(self, int, neighbours, contig_appearances):
        if len(self.gene_detector.contigs_to_genes()[neighbours[0][int]].gene_objects) in range(1, neighbours[1] + 3):
            if contig_appearances[neighbours[0][int]][0] == len(self.gene_detector.contigs_to_genes()[neighbours[0][int]].gene_objects):
                return True    
        
    def check_if_connection_is_valid(self, neighbours):
        """Only keeps the cases where there is not a gene by itself or too many genes, coming from an intersection towards the middle of the contig"""
        contig_appearances = self.find_contig_appearances(neighbours)
        most_occurent_contig = self.most_occurent_contig(contig_appearances, neighbours)
        
        if  contig_appearances[neighbours[0][0]][0] in range(min(2, len(self.gene_detector.contigs_to_genes()[neighbours[0][0]].gene_objects)), neighbours[1] + 4) and\
                         contig_appearances[neighbours[0][1]][0] in range(min(2, len(self.gene_detector.contigs_to_genes()[neighbours[0][1]].gene_objects)), neighbours[1] + 4) and\
                         most_occurent_contig in [neighbours[0][0], neighbours[0][1]]:
            if len(contig_appearances) == 2:
                return neighbours
            else:
                if self.check_for_small_contig(0,neighbours,contig_appearances) or self.check_for_small_contig(1,neighbours,contig_appearances):
                    return neighbours
            
    def refine_contig_neighbours(self):
        """Removes the case when we have a short contig between two others."""
        pprint.pprint(self.neighbouring_contigs)
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
        pprint.pprint(self.refined_neighbouring_contigs)
        return self.refined_neighbouring_contigs
    
    
    def ends_of_contigs(self):
        for neighbours in self.neighbouring_contigs:
            contig_appearances = self.find_contig_appearances(neighbours)
            for contig in [neighbours[0][0], neighbours[0][1]]:
                if contig not in self.contig_ends:
                    self.contig_ends[contig] = {}
                    
            self.contig_ends[neighbours[0][0]][neighbours[0][1]] = self.getting_query_start_tuples(neighbours, contig_appearances, 0)
            self.contig_ends[neighbours[0][1]][neighbours[0][0]] = self.getting_query_start_tuples(neighbours, contig_appearances, 1)

        return self.contig_ends
    
    
    def remove_gene_degree_one(self, neighbours, contig_appearances, int, i):
        genes_degree_at_least_2 = []
        #print(contig_appearances[neighbours[0][int]][1])
        for gene in contig_appearances[neighbours[0][int]][1][i]:
            if self.filtered_graph.degree(gene) != 1:
                genes_degree_at_least_2.append(gene)
        return genes_degree_at_least_2
                        
    
    def getting_query_start_tuples(self, neighbours, contig_appearances, int):
        qry_starts = []
        iteration_gene_dict = contig_appearances[neighbours[0][int]][1]
        iterations = sorted(iteration_gene_dict.keys())
        gene_objects = self.gene_detector.contigs_to_genes()[neighbours[0][int]].gene_objects
        #print(neighbours[0][int])
        #print(contig_appearances)
        genes_degree_at_least_2_first_iteration = self.remove_gene_degree_one(neighbours, contig_appearances, int, iterations[0])
        sequence_length = self.sequences[neighbours[0][int]][2]
        
        if len(iterations) >= 2:
            
            if len(genes_degree_at_least_2_first_iteration) == 1:
                genes_degree_at_least_2_second_iteration = self.remove_gene_degree_one(neighbours, contig_appearances, int, iterations[1])
                if len(genes_degree_at_least_2_second_iteration) == 1:
                    qry_starts.append(gene_objects[genes_degree_at_least_2_first_iteration[0]].qry_start)    
                    qry_starts.append(gene_objects[genes_degree_at_least_2_second_iteration[0]].qry_start)
                    return qry_starts
                else:
                    if  all(gene_objects[iteration_gene_dict[iterations[1]][i]].qry_start < \
                            gene_objects[iteration_gene_dict[iterations[0]][0]].qry_start \
                            for i in range(len(iteration_gene_dict[iterations[1]]))) or \
                        all(gene_objects[iteration_gene_dict[iterations[1]][i]].qry_start >= \
                            gene_objects[iteration_gene_dict[iterations[0]][0]].qry_start \
                            for i in range(len(iteration_gene_dict[iterations[1]]))):
                                qry_starts.append(gene_objects[iteration_gene_dict[iterations[0]][0]].qry_start)    
                                qry_starts.append(gene_objects[iteration_gene_dict[iterations[1]][0]].qry_start)
                                return qry_starts
                    else:
                        if gene_objects[iteration_gene_dict[iterations[0]][0]].qry_start < sequence_length/2:
                            return [gene_objects[iteration_gene_dict[iterations[0]][0]].qry_start, sequence_length]
                        else:
                            return [gene_objects[contig_appearances[neighbours[0][int]][1][iterations[0]][0]].qry_start, 1]
            else:
                
                if all(gene_objects[iteration_gene_dict[iterations[0]][i]].qry_start < \
                        sequence_length/2 for i in range(len(iteration_gene_dict[iterations[0]]))):
                        return [gene_objects[iteration_gene_dict[iterations[0]][0]].qry_start, sequence_length]
                elif all(gene_objects[iteration_gene_dict[iterations[1]][i]].qry_start >= \
                        sequence_length/2 for i in range(len(iteration_gene_dict[iterations[1]]))):
                        return [gene_objects[iteration_gene_dict[iterations[0]][0]].qry_start, 1]
                else:
                    if len(iterations) > 1:
                        if all(gene_objects[iteration_gene_dict[iterations[1]][i]].qry_start < \
                            sequence_length/2 for i in range(len(iteration_gene_dict[iterations[1]]))):
                            return [gene_objects[iteration_gene_dict[iterations[1]][0]].qry_start, sequence_length]
                        elif all(gene_objects[iteration_gene_dict[iterations[1]][i]].qry_start >= \
                            sequence_length/2 for i in range(len(iteration_gene_dict[iterations[1]]))):
                            return [gene_objects[iteration_gene_dict[iterations[1]][0]].qry_start, 1]
                        else:
                            return [None, None]
        elif len(iterations) == 1:
            if gene_objects[iteration_gene_dict[iterations[0]][0]].qry_start < sequence_length/2:
                return [gene_objects[iteration_gene_dict[iterations[0]][0]].qry_start, sequence_length]
            else:
                return [gene_objects[iteration_gene_dict[iterations[0]][0]].qry_start, 1]            
        else:
            return [None, None]
            
    
    
    
    
    
        
        
        
        