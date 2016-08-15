import networkx as nx
from madansi.Assembly import Assembly

class GraphToFasta(object):
    
    def __init__(self,input_fasta_file, graph, output_fasta_fname):
        self.input_fasta_file = input_fasta_file
        self.graph = graph
        self.output_fasta_fname = output_fasta_fname

    def find_sequences(self):
        assembly = Assembly(self.input_fasta_file)
        assembly.sequence_names()
        return assembly.sequences
        
    def find_start_contig(self, component):

        for contig in component:
            if self.graph.degree(contig) == 1:
                return contig
    
    def walk_contig_graph(self, component):
        
        contig = self.find_start_contig(component)
        visited = [contig]
        
        while self.graph.degree(contig)<= 2 and sorted(visited) != sorted(component):
            for neighbour in self.graph.neighbors(contig):
                if neighbour not in visited:
                    visited.append(neighbour)
                    contig = neighbour
        return visited
    
    def create_fasta_file(self):
        sequences = self.find_sequences()
        f = open(self.output_fasta_fname, 'w')
        for component in sorted(nx.connected_components(self.graph), key = len, reverse=True):
            visited = self.walk_contig_graph(list(component))
            for contig in visited:
                f.write('> ' + contig + '\n')
                f.write(str(sequences[contig]) + '\n')
                print('added ' + contig)
                
        f.close()
        

            
    

        
