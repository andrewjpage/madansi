from madansi.Assembly import Assembly

class UnusedContigs(object):
    def __init__(self,gene_detector, output_fasta_file, input_fasta_file): 
        self.gene_detector = gene_detector
        self.unused_contigs_list = []
        self.output_fasta_file = output_fasta_file
        self.input_fasta_file = input_fasta_file
        assembly = Assembly(self.input_fasta_file)
        self.sequence_list = assembly.sequence_names()
        self.sequences = assembly.sequences
        
        
    def contigs_not_in_filtered_file(self):
        for contig in self.sequence_list:
            if self.gene_detector.contigs[contig].gene_objects == {}:
                if contig not in self.unused_contigs_list:
                    self.unused_contigs_list.append(contig)       
        return self.unused_contigs_list
    
        
    def contigs_not_in_filtered_graph(self,filtered_graph):
        contigs_present_in_filtered_graph = filtered_graph.nodes()
        for contig in self.sequence_list:
            if  contig  not in contigs_present_in_filtered_graph and \
                contig  not in self.unused_contigs_list:
                self.unused_contigs_list.append(contig)
        return self.unused_contigs_list
    
    def add_unused_contigs_to_end(self):
        f = open(self.output_fasta_file, 'a')
        for unused_contig in self.unused_contigs_list:
            f.write('>' + unused_contig + '\n')
            f.write(str(self.sequences[unused_contig][0]) + '\n')
        f.close()
    
    
