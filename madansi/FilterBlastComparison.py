from madansi.BlastHit import BlastHit, file_reader, Error

class FilterBlastComparison(object):
    """Filter output of BLAST comparison"""
    def __init__(self, input_blast_file, filtered_file, percent_identity=0.0, alignment_length=0, mismatches=10000, gap_openings=10, evalue=10, bit_score=0):
        self.input_blast_file = input_blast_file
        self.filtered_file = filtered_file
        try:
            self.percent_identity = float(percent_identity)
        except percentidentity > 100 or percent_identity < 0:
            raise ValueError("Percent identity should be a float or int between 0 and 100")
        self.percent_identity = float(percent_identity)
        self.alignment_length = int(alignment_length)
        self.mismatches = int(mismatches)
        self.gap_openings = int(gap_openings)
        self.evalue = float(evalue)
        self.bit_score = float(bit_score)
          
    
    def find_gene_duplicates(self):
        """Filters the output from a BLAST comparison by neglecting any genes that appear on three or more contigs"""
        try:
            f = open(self.input_blast_file)
        except IOError:
            raise Error('Error opening this file')
        
        gene_to_contigs = {}
        gene_duplicates = []
        
        for line in f:
            bh = BlastHit(line)
            if bh.ref_name not in gene_to_contigs:
                gene_to_contigs[bh.ref_name] = [bh.qry_name]
            elif bh.qry_name not in gene_to_contigs[bh.ref_name]:
                gene_to_contigs[bh.ref_name].append(bh.qry_name)
        
        for gene in gene_to_contigs:
            if len(gene_to_contigs[gene]) > 1:
                gene_duplicates.append(gene)
        
        return gene_duplicates
        
    
    def filter(self):
        try:
            f = open(self.input_blast_file)
        except IOError:
            raise Error('Error opening this file')
        
        filtered_output = open(self.filtered_file, 'a')
        gene_duplicates = self.find_gene_duplicates()
        
        for line in f:
            bh = BlastHit(line)
            if bh.alignment_length >= self.alignment_length and bh.mismatches <= self.mismatches and bh.gap_openings <= self.gap_openings \
            and bh.e_value <= self.evalue and bh.percent_identity >= self.percent_identity and bh.bit_score >= self.bit_score and \
            bh.ref_name not in gene_duplicates:
                filtered_output.write(line)
        
        f.close()
        filtered_output.close()
        
        
        