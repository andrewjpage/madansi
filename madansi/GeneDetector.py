from madansi.Contig import Contig
from madansi.Assembly import Assembly
import os
from madansi.BlastHit import file_reader, BlastHit

class GeneDetector(object):
    def __init__(self, input_assembly_file, blast_hits_file):
        self.input_assembly_file = input_assembly_file
        self.blast_hits_file = blast_hits_file
        self.assembly = Assembly(self.input_assembly_file)
        
    def parse_blast_hits(self):
        hits = []
        for hit in file_reader(self.blast_hits_file):
            hits.append(hit)
        return hits
            
    def contigs_to_genes(self):
        contigs = {}
        for sequence_name in self.assembly.sequence_names():
            contigs[sequence_name]= Contig(sequence_name)
        for hit in self.parse_blast_hits():
            sequence_name = hit.qry_name
            contigs[sequence_name].add_blast_hit(hit)
        return contigs
            
            

    
    