from Bio import SeqIO

class Assembly(object):
    """Gives a list of all the contigs present in the sample"""
    def __init__(self, input_file):
        self.input_file = input_file
    
    def sequence_names(self):
        sequence_list = []
        with open(self.input_file, "r") as handle:
            for seq_record in SeqIO.parse( handle, "fasta"):
                sequence_list.append(seq_record.id)
            handle.close()
        return sequence_list
       

    
      
                

        