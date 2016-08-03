from madansi.Contig import Contig
from madansi.Assembly import Assembly
import os
from madansi.BlastHit import file_reader, BlastHit

class GeneDetector(object):
    """An object that has all the contigs present as well as the genes within each contig and their information"""
    def __init__(self, input_assembly_file, filtered_blast_hits_file):
        self.input_assembly_file = input_assembly_file
        self.filtered_blast_hits_file = filtered_blast_hits_file
        self.assembly = Assembly(self.input_assembly_file)
        self.contigs = self.contigs_to_genes()
        
    def parse_blast_hits(self):
        hits = []
        for hit in file_reader(self.filtered_blast_hits_file):
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
            
            

    
    