from madansi.BlastHit import BlastHit
class Error (Exception): pass

class GenePresent(object):
    
    def __init__(self,filteredfile):
        self.filteredfile = filteredfile
    
    def construct_dictionary(self):
        """Constructs a dictionary based upon the genes that are thought to be present by filtering using bitscore"""
        gene_present_dict = {}
        
        try:
            f = open(self.filteredfile)
        except IOError:
            raise Error("Error opening this file")
        
        for line in f:
            bh = BlastHit(line)
            if bh.bit_score >= 150:
                gene_present_dict[bh.qry_name] = True
            else:
                gene_present_dict[bh.qry_name] = False
        
        f.close() 
        return(gene_present_dict)    
    
    def index_filtered_file(self):
        """Constructs a dictionary linking the name of the gene with the line of the filtered file that it is on"""
        index_lines = {}
        count = 0
        
        try:
            f = open(self.filteredfile)
        except IOError:
            raise Error("Error opening this file")
            
        for line in f:
            l = line.split("\t")
            index_lines[l[0]] = count
            count += 1
        
        f.close()
        return index_lines
        
def list_genes_present(dictionary):
    """Extract all genes that are thought to be present"""
    list_of_genes_present = [k for k,v in dictionary.items() if v]
    return(list_of_genes_present)