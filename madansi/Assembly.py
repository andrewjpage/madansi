from Bio import SeqIO


class Assembly(object):
    def __init__(self, inputfile):
        self.inputfile = inputfile
    
    def sequence_names(self):
        sequence_list = []
        with open(self.inputfile, "r") as handle:
            for seq_record in SeqIO.parse( handle, "fasta"):
                sequence_list.append(seq_record.id)
            handle.close()
        return sequence_list
       
		
    
      
                

        