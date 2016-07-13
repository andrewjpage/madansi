from madansi.SwitchColumns import SwitchColumns
from Bio.Blast.Applications import NcbiblastnCommandline
import os

class RunBLAST(object):
    def __init__(self,inputreference,outputreference, outputdatabase, queryfile, finaloutput, evalue = 1e-20):
        self.inputreference = inputreference
        self.outputreference = outputreference
        self.outputdatabase = outputdatabase
        self.queryfile = queryfile
        self.finaloutput = finaloutput
        self.evalue = evalue
        
    def run_switch_columns(self):
        """Switches the columns in the reference fasta file"""
        sw = SwitchColumns(self.inputreference, self.outputreference)
        sw.run()
        
    def make_reference_database(self):
        """Makes the database from the modified reference fasta file """
        makedbcmd = 'makeblastdb -dbtype nucl -in ' +  self.outputreference + ' -out ' + self.outputdatabase
        os.system(makedbcmd)
        
    def run_BLAST(self):
        """Run BLAST with the query file and the modified reference fasta file"""
        blastn_cline = NcbiblastnCommandline(query = self.queryfile, db= self.outputdatabase, outfmt=6, out=self.finaloutput, task ='blastn', evalue= self.evalue)
        print(" ".join([self.queryfile, self.outputdatabase, self.finaloutput]))
        stdout, stderr = blastn_cline()
        print(stdout)
        print(stderr)