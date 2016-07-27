from madansi.SwitchColumns import SwitchColumns
from Bio.Blast.Applications import NcbiblastnCommandline
import os

class RunBLAST(object):
    def __init__(self,query,inputreference, outputreference, outputdatabase, blast_output, evalue = 0.01):
        self.query = query
        self.inputreference = inputreference
        self.outputreference = outputreference
        self.outputdatabase = outputdatabase
        self.blast_output = blast_output
        self.evalue = evalue
        
    def run_switch_columns_database(self):
        """Switches the columns in the reference fasta file"""
        sw = SwitchColumns(self.inputreference, self.outputreference)
        sw.run()
        
    def make_reference_database(self):
        """Makes the database from the modified reference fasta file """
        makedbcmd = 'makeblastdb -dbtype nucl -in ' +  self.outputreference + ' -out ' + self.outputdatabase
        os.system(makedbcmd)
        
    def run_BLAST(self):
        """Run BLAST with the query file and the modified reference fasta file"""
        blastn_cline = NcbiblastnCommandline(query = self.query, db= self.outputdatabase, outfmt=6, out=self.blast_output, task ='blastn', evalue= self.evalue)
        stdout, stderr = blastn_cline()
