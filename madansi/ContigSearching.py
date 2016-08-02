from madansi.NeighboursOfNodes import NeighboursOfNodes

class ContigSearching(object):
    def __init__(self, gene_detector, filtered_graph):
        self.gene_detector = gene_detector
        self.filtered_graph = filtered_graph
        self.genes_in_contig_radius = {}
        self.neighbouring_contigs = []    
        self.finished_contigs = set()
        self.found_contigs = set()
    
    def contig_expansion(self):
        for sequence_name, genes in self.gene_detector.contigs.items():
            gene_names = [gene for gene in genes.gene_objects] 
            self.set_expansion(gene_names, sequence_name)
        
    def neighbourhood_expansion(self):
        for sequence_name, genes in self.genes_in_contig_radius.items():
            self.set_expansion(genes, sequence_name)
            
    def set_expansion(self, gene_names, sequence_name):
        neighbouring_nodes = NeighboursOfNodes(self.filtered_graph).find_neighbours(gene_names)
        self.genes_in_contig_radius[sequence_name] = gene_names + neighbouring_nodes
        if neighbouring_nodes == []:
            if sequence_name not in self.finished_contigs:
                self.finished_contigs.add(sequence_name)
        return self
        
    def expand_all_contigs(self):
        iteration_count = 1
        self.contig_expansion()
        self.check_intersections(iteration_count) 
        while set([sequence_name for sequence_name in self.gene_detector.contigs]) - self.finished_contigs != set() and \
        set([sequence_name for sequence_name in self.gene_detector.contigs]) - self.found_contigs != set():
            print(iteration_count)
            iteration_count += 1
            self.neighbourhood_expansion()
            self.check_intersections(iteration_count)    
        return self.neighbouring_contigs  
    
    def check_intersections(self, iteration_count):
        for sequence_name_1, gene_names_1 in self.genes_in_contig_radius.items():
            for sequence_name_2, gene_names_2 in self.genes_in_contig_radius.items():
                if gene_names_1 == [] or gene_names_2 == []:
                    continue
                if sequence_name_1 == sequence_name_2:
                    continue
                if set(gene_names_1) & set(gene_names_2) != set():
                    if not any( (sequence_name_1, sequence_name_2) == (entry[0], entry[1]) or (sequence_name_1, sequence_name_2) == (entry[1], entry[0]) for entry in self.neighbouring_contigs):
                        self.neighbouring_contigs.append((sequence_name_1, sequence_name_2,iteration_count))
                        if sequence_name_1 not in self.found_contigs:
                            self.found_contigs.add(sequence_name_1)
                        if sequence_name_2 not in self.found_contigs:
                            self.found_contigs.add(sequence_name_2)
        return self
    
                