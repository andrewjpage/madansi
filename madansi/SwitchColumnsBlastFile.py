from madansi.BlastHit import Error

class SwitchColumnsBlastFile(object):
    def __init__(self,input_file, output_file):
        self.input_file = input_file
        self.output_file = output_file
    
    def switch_columns_blast_file(self):
        try:
            f = open(self.input_file)
        except IOError:
            raise Error('Error opening this file')
       
        with open(self.output_file, 'a') as handle:
            for line in f:
                handle.write(self.rejoin_string(line) + '\n')
        
        handle.close()
        f.close()
        
    
    def rejoin_string(self, line):
        
        try:
            l = line.rstrip().split('\t')
            qry_name = l[0]
            ref_name = l[1]
            percent_identity = float(l[2])
            alignment_length = int(l[3])
            mismatches = int(l[4])
            gap_openings = int(l[5])
            qry_start = int(l[6])
            qry_end = int(l[7])
            ref_start = int(l[8])
            ref_end = int(l[9])
            e_value = float(l[10])
            bit_score = float(l[11])
        except:
            raise Error('Error reading this blast line:\n' + line)
        
        s = '\t'.join(str(x) for x in
            [ref_name,
             qry_name,
             '{:.2f}'.format(percent_identity),
             alignment_length,
             mismatches,
             gap_openings,
             qry_start,
             qry_end,
             ref_start,
             ref_end,
             e_value,
             '{:g}'.format(bit_score)])

        return s
