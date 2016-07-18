from madansi.BlastHit import BlastHit, file_reader, Error

class FilterBlastComparison2(object):
    """Filter output of BLAST comparison"""
    def __init__(self, comparisonfile, filteredfile, percentidentity=0.0, alignmentlength=0, mismatches=0, gapopenings=0, evalue=10, bitscore=0):
        self.comparisonfile = comparisonfile
        self.filteredfile = filteredfile
        self.percentidentity = float(percentidentity)
        self.alignmentlength = alignmentlength
        self.mismatches = mismatches
        self.gapopenings = gapopenings
        self.evalue = float(evalue)
        self.bitscore = float(bitscore)
    
    def filter_percent_identity(self):
        """Filters the output from a BLAST comparison by including all lines with percentage identity at least the given value"""
        try:
            f = open(self.comparisonfile)
        except IOError:
            raise Error('Error opening this file')
            
        filteredoutput = open(self.filteredfile, 'a')
        
        for line in f:
            bh = BlastHit(line)
            if bh.percent_identity >= self.percentidentity:
                filteredoutput.write(line)
        
        f.close()
        filteredoutput.close()
              
    def filter_alignment_length(self):
        """Filters the output from a BLAST comparison by including all lines with alignment length at least the given value"""
        try:
            f = open(self.comparisonfile)
        except IOError:
            raise Error('Error opening this file')
            
        filteredoutput = open(self.filteredfile, 'a')
        
        for line in f:
            bh = BlastHit(line)
            if bh.alignment_length >= self.alignmentlength:
                filteredoutput.write(line)
        
        f.close()
        filteredoutput.close()
    
    def filter_mismatches(self):
        """Filters the output from a BLAST comparison by including all lines with mismatches at most the given value"""
        try:
            f = open(self.comparisonfile)
        except IOError:
            raise Error('Error opening this file')
            
        filteredoutput = open(self.filteredfile, 'a')
        
        for line in f:
            bh = BlastHit(line)
            if bh.mismatches <= self.mismatches:
                filteredoutput.write(line)
        
        f.close()
        filteredoutput.close()
    
    def filter_gapopenings(self):
        """Filters the output from a BLAST comparison by including all lines with gap openings at least the given value"""
        try:
            f = open(self.comparisonfile)
        except IOError:
            raise Error('Error opening this file')
            
        filteredoutput = open(self.filteredfile, 'a')
        
        for line in f:
            bh = BlastHit(line)
            if bh.gap_openings >= self.gapopenings:
                filteredoutput.write(line)
        
        f.close()
        filteredoutput.close()       
         
    def filter_evalue(self):
        """Filters the output from a BLAST comparison by including all lines with expect value at most the given value"""
        try:
            f = open(self.comparisonfile)
        except IOError:
            raise Error('Error opening this file')
            
        filteredoutput = open(self.filteredfile, 'a')
        
        for line in f:
            bh = BlastHit(line)
            if bh.e_value <= self.evalue:
                filteredoutput.write(line)
        
        f.close()
        filteredoutput.close()        
        
    def filter_bit_score(self):
        """Filters the output from a BLAST comparison by including all lines with bit score at least the given value"""
        try:
            f = open(self.comparisonfile)
        except IOError:
            raise Error('Error opening this file')
            
        filteredoutput = open(self.filteredfile, 'a')
        
        for line in f:
            bh = BlastHit(line)
            if bh.bit_score >= self.bitscore:
                filteredoutput.write(line)
        
        f.close()
        filteredoutput.close()