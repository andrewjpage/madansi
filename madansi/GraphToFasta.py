import networkx as nx
from madansi.Assembly import Assembly
from madansi.DetermineOrientation import DetermineOrientation
import pprint
import sys

class GraphToFasta(object):
    
    def __init__(self,input_fasta_file, graph, output_fasta_fname, contig_ends):
        self.input_fasta_file = input_fasta_file
        self.graph = graph
        self.output_fasta_fname = output_fasta_fname
        self.contig_ends = contig_ends
        pprint.pprint(contig_ends)
        self.contig_orientation = {}
        assembly = Assembly(self.input_fasta_file)
        assembly.sequence_names()
        self.sequences = assembly.sequences
        self.combined_contigs_dict = {}
        
    def contigs_degree_one(self, component):
        contigs_degree_one = []
        for contig in component:
            if self.graph.degree(contig) == 1:
                contigs_degree_one.append(contig)
        return contigs_degree_one 

    def find_start_contig(self, component):
        contigs_degree_one = self.contigs_degree_one(component)
        try:
            start_contig = sorted(contigs_degree_one)[0] 
            self.contig_orientation[start_contig] = DetermineOrientation().determine_orientation_start_contig(start_contig,self.contig_ends, [], self.contig_orientation, self.graph, False, self.sequences)
        except IndexError:
            start_contig = sorted(list(component))[0]
            self.contig_orientation[start_contig] = DetermineOrientation().determine_orientation_start_contig(start_contig,self.contig_ends, [], self.contig_orientation, self.graph, True, self.sequences)
        return start_contig
    
    def walk_contig_graph(self, component):
        contig = self.find_start_contig(component)
        visited = [contig]
        ends_seen = [contig]
        while self.graph.degree(contig)<= 2 and sorted(visited) != sorted(component):
            for neighbour in self.graph.neighbors(contig):
                if neighbour not in visited:
                    DetermineOrientation().determine_orientation(neighbour, self.contig_ends, visited, self.contig_orientation, self.graph, False, self.sequences)
                    visited.append(neighbour)
                    contig = neighbour
        return visited
        
    def combine_contigs(self, component, contig_count):
        combined_contig = ''
        visited = self.walk_contig_graph(component)
        count = 1
        for contig in visited:
            if self.contig_orientation[contig] == 1:
                combined_contig += str(self.sequences[contig][0])
            else:
                combined_contig += str(self.sequences[contig][1])
            if count != len(visited):
                for i in range(200):
                    combined_contig += 'N'
            count += 1
        self.combined_contigs_dict[contig_count] = visited
        pprint.pprint(self.combined_contigs_dict)
        return combined_contig
        
    def create_fasta_file_combined_contigs(self):
        f = open(self.output_fasta_fname, 'w')
        contig_count = 1
        for component in sorted(nx.connected_components(self.graph), key = len, reverse=True):
            combined_contig = self.combine_contigs(component, contig_count)
            f.write('>Contig'+ str(contig_count) + '\n')
            f.write(combined_contig + '\n')
            contig_count += 1
        f.close()
            