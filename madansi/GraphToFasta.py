import networkx as nx
from madansi.DetermineOrientation import DetermineOrientation

class GraphToFasta(object):
    
    def __init__(self, sequences, graph, output_fasta_fname, contig_ends):
        self.graph = graph
        self.output_fasta_fname = output_fasta_fname
        self.contig_ends = contig_ends
        self.contig_orientation = {}
        self.sequences = sequences
        self.combined_contigs_dict = {}
        
    def contigs_degree_one(self, component):
        contigs_degree_one = []
        for contig in component:
            if self.graph.degree(contig) == 1:
                contigs_degree_one.append(contig)
        return contigs_degree_one 
    
    def is_cycle(self, component):
        return len(self.contigs_degree_one(component)) == 0

    def find_start_contig(self, component):
        contigs_degree_one = self.contigs_degree_one(component)
        cycle = self.is_cycle(component)
        if not cycle:
            start_contig = sorted(contigs_degree_one)[0] 
            self.contig_orientation[start_contig] = DetermineOrientation().determine_orientation_start_contig(start_contig,self.contig_ends, [], self.contig_orientation, self.graph, cycle, self.sequences)
        else:
            start_contig = sorted(list(component))[0]                      
            self.contig_orientation[start_contig] = 1
        return start_contig
    
    def walk_contig_graph(self, component):
        contig = self.find_start_contig(component)
        visited = [contig]
        orientation_assigned = [contig]
        while sorted(orientation_assigned) != sorted(component):
            if len(visited) == 1 and self.is_cycle(component):
                neighbour = DetermineOrientation().determine_initial_direction_cycle(contig, self.contig_ends, visited, self.contig_orientation, self.graph, True, self.sequences)
                DetermineOrientation().determine_orientation(neighbour, self.contig_ends, visited, self.contig_orientation, self.graph, False, self.sequences)
            else:
                for neighbour in self.graph.neighbors(contig):
                    if neighbour not in visited:
                        DetermineOrientation().determine_orientation(neighbour, self.contig_ends, visited, self.contig_orientation, self.graph, False, self.sequences)
                        break
            orientation_assigned.append(neighbour)
            visited.append(neighbour)
            contig = neighbour
        return visited
        
    def combine_contigs(self, component, contig_count):
        combined_contig = ''
        visited = self.walk_contig_graph(component)
        count = 1
        
        for contig in visited:
            if self.contig_orientation[contig] == 1:
                combined_contig = combined_contig + str(self.sequences[contig][0])
            else:
                combined_contig = combined_contig + str(self.sequences[contig][1])
            if count != len(visited):
                for i in range(200):
                    combined_contig = combined_contig + 'N'
            count += 1
        self.combined_contigs_dict[contig_count] = visited
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
            
