from madansi.SwitchColumns import SwitchColumns
from Bio.Blast.Applications import NcbiblastnCommandline
import os

class RunBLAST(object):
    def __init__(self,query,input_reference, output_reference, output_database, blast_output, evalue = 0.01):
        self.query = query
        self.input_reference = input_reference
        self.output_reference = output_reference
        self.output_database = output_database
        self.blast_output = blast_output
        self.evalue = evalue
        
    def run_switch_columns_database(self):
        """Switches the columns in the reference fasta file"""
        sw = SwitchColumns(self.input_reference, self.output_reference)
        sw.run()
        
    def make_reference_database(self):
        """Makes the database from the modified reference fasta file """
        makedbcmd = 'makeblastdb -dbtype nucl -in ' +  self.output_reference + ' -out ' + self.output_database
        os.system(makedbcmd)
        
    def run_BLAST(self):
        """Run BLAST with the query file and the modified reference fasta file"""
        blastn_cline = NcbiblastnCommandline(query = self.query, db= self.output_database, outfmt=6, out=self.blast_output, task ='blastn', evalue= self.evalue)
        stdout, stderr = blastn_cline()
