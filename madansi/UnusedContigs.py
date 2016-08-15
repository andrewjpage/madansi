import networkx as nx

class UnusedContigs(object):
    def __init__(self,gene_detector, sequence_list, unused_contigs_file): 
        self.gene_detector = gene_detector
        self.sequence_list = sequence_list
        self.unused_contigs_list = []
        self.unused_contigs_file = unused_contigs_file
        
    def contigs_not_in_filtered_file(self):
        for contig in self.sequence_list:
            if self.gene_detector.contigs[contig].gene_objects == {}:
                if contig not in self.unused_contigs_list:
                    self.unused_contigs_list.append(contig)       
        return self.unused_contigs_list
    
    def output_unused_contigs(self):
        f = open(self.unused_contigs_file, 'w')
        for unused_contig in self.unused_contigs_list:
            f.write(unused_contig + '\n')
        f.close()
        
    def contigs_not_in_refined_graph(self,refined_graph):
        contigs_present_in_refined_graph = refined_graph.nodes()
        for contig in self.sequence_list:
            if  contig  not in contigs_present_in_refined_graph and \
                contig  not in self.unused_contigs_list:
                self.unused_contigs_list.append(contig)
        return self.unused_contigs_list
        
    