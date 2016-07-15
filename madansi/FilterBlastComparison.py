from madansi.BlastHit import BlastHit, file_reader, Error

class FilterBlastComparison(object):
    """Filter output of BLAST comparison"""
    def __init__(self, comparisonfile, filteredfile, percentidentity=0.0, alignmentlength=0, mismatches=10, gapopenings=10, evalue=10, bitscore=0):
        self.comparisonfile = comparisonfile
        self.filteredfile = filteredfile
        self.percentidentity = float(percentidentity)
        self.alignmentlength = int(alignmentlength)
        self.mismatches = int(mismatches)
        self.gapopenings = int(gapopenings)
        self.evalue = float(evalue)
        self.bitscore = float(bitscore)
          
    def filter(self):
        """Filters the output from a BLAST comparison by including all lines with alignment length at least the given value"""
        try:
            f = open(self.comparisonfile)
        except IOError:
            raise Error('Error opening this file')
            
        filteredoutput = open(self.filteredfile, 'a')
        
        for line in f:
            bh = BlastHit(line)
            if bh.alignment_length >= self.alignmentlength and bh.mismatches <= self.mismatches and bh.gap_openings <= self.gapopenings \
            and bh.e_value <= self.evalue and bh.percent_identity >= self.percentidentity:
                filteredoutput.write(line)   
                    
        f.close()
        filteredoutput.close()