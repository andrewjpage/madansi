import networkx as nx
from madansi.Assembly import Assembly
import pprint

class GraphToFasta(object):
    
    def __init__(self,input_fasta_file, graph, output_fasta_fname, contig_ends):
        self.input_fasta_file = input_fasta_file
        self.graph = graph
        self.output_fasta_fname = output_fasta_fname
        self.contig_ends = contig_ends
        self.contig_orientation = {}

    def find_sequences(self):
        assembly = Assembly(self.input_fasta_file)
        assembly.sequence_names()
        return assembly.sequences
        
    def find_contigs_degree_one(self, component):
        contigs_degree_one = []
        for contig in component:
            if self.graph.degree(contig) == 1:
                contigs_degree_one.append(contig)
        return contigs_degree_one
    
    def determine_orientation(self, visited, contig):
        
        if len(self.contig_ends[contig]) == 2:
            for neighbour in self.graph.neighbors(contig):
                if neighbour in visited:
                    visited_contig = neighbour
                else:
                    unvisited_contig = neighbour
            if self.contig_ends[contig][visited_contig] <= self.contig_ends[contig][unvisited_contig]:
                self.contig_orientation[contig] = 1
            else:
                self.contig_orientation[contig] = -1
        else:
            
            
            self.contig_orientation[contig] = 1
        return self.contig_orientation[contig]

    def walk_contig_graph(self, component):
        
        contigs_degree_one = self.find_contigs_degree_one(component)
        contig = sorted(contigs_degree_one)[0]
        self.contig_orientation[contig] = 1
        visited = [contig]
        ends_seen = [contig]
        
        while self.graph.degree(contig)<= 2 and sorted(visited) != sorted(component):
            for neighbour in self.graph.neighbors(contig):
                if neighbour not in visited:
                    self.determine_orientation(visited, neighbour)
                    visited.append(neighbour)
                    contig = neighbour

        return visited
    
    def create_fasta_file(self):
        
        sequences = self.find_sequences()
        f = open(self.output_fasta_fname, 'w')
        for component in sorted(nx.connected_components(self.graph), key = len, reverse=True):
            visited = self.walk_contig_graph(list(component))
            pprint.pprint(self.contig_orientation)
            for contig in visited:
                f.write('>' + contig + '\n')
                if self.contig_orientation[contig] == 1:
                    f.write(str(sequences[contig][0]) + '\n')
                else:
                    f.write(str(sequences[contig][1]) + '\n')
        f.close()
        
    def combine_contigs(self, split_components):
        sequences = self.find_sequences()
        for split_component in split_components:
            pass
            
            

            
    

        
