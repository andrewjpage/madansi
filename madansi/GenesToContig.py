from madansi.BlastHit import file_reader, BlastHit
class GenesToContig(object):
    
    def __init__(self, filtered_blast_file):
        self.filtered_blast_file = filtered_blast_file
        self.genes = {}
    
    def parse_blast_hits(self):
        hits = []
        for hit in file_reader(self.filtered_blast_file):
            hits.append(hit)
        return hits
    
    def genes_to_contig(self):
        genes = {}
        for hit in self.parse_blast_hits():
            gene_name = hit.ref_name
            genes[gene_name] = hit.qry_name
        return genes
        