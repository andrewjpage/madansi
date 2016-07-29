from madansi.BlastHit import BlastHit
from madansi.Gene import Gene

class Contig(object):
    
    def __init__(self, sequence_name):
        self.sequence_name = sequence_name
        self.gene_objects = {}
    
    def genes(self):
        return self.gene_objects
    
    def add_blast_hit(self, blast_hit):
        gene_name = blast_hit.ref_name
        if gene_name in self.gene_objects:
            return False
        my_gene = Gene(blast_hit.orientation(), blast_hit.ref_start, blast_hit.ref_end, None, self.sequence_name)
        self.gene_objects[gene_name] = my_gene
        return True
    