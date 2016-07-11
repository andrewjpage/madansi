from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

class SwitchColumns(object):
	
    def __init__(self,inputfile,outputfile):
        self.inputfile = inputfile
        self.outputfile = outputfile
		
    def run(self):
        my_records = []
        with open(self.inputfile, "r") as handle:
            for seq_record in SeqIO.parse( handle, "fasta"):
                desc = str(seq_record.description)
                n = desc.index(" ")
                Name = desc[:n]
                group = desc[n+1:]
                rec = SeqRecord(seq_record.seq, id = group, description = group + " " + Name )
                my_records.append(rec)
            handle.close()
        SeqIO.write(my_records, self.outputfile, "fasta")
	
    
	