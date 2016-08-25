class RemoveContigs(object):
    def __init__(self, graph):
        self.graph = graph
        
    
    def remove_contigs_high_degree(self):
        for contig in self.graph.nodes():
            if self.graph.degree(contig) > 2:
                self.graph.remove_node(contig)
        return self.graph
    
    def remove_extra_contigs(self):
        self.remove_contigs_high_degree()
        for contig in self.graph.nodes():
            if self.graph.degree(contig) == 0:
                self.graph.remove_node(contig)
                
        return self.graph
            
