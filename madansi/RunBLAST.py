from madansi.SwitchColumns import SwitchColumns
from Bio.Blast.Applications import NcbiblastnCommandline
import tempfile
import os

class RunBLAST(object):
    def __init__(self,query,input_reference, temp_dir ,evalue = 0.01):
        self.query = query 
        self.input_reference = input_reference
        self.output_reference = tempfile.NamedTemporaryFile(delete = False, dir= temp_dir)
        self.output_database = tempfile.NamedTemporaryFile(delete = False, dir= temp_dir)
        self.blast_output = tempfile.NamedTemporaryFile(delete = False, dir= temp_dir)
        self.evalue = evalue
        self.temp_dir = temp_dir
        
    def run_switch_columns_database(self):
        """Switches the columns in the reference fasta file"""
        sw = SwitchColumns(self.input_reference, self.output_reference.name)
        sw.run()
        
    def make_reference_database(self):
        """Makes the database from the modified reference fasta file """
        makedbcmd = 'makeblastdb -dbtype nucl -in ' +  self.output_reference.name + ' -out ' + self.output_database.name
        os.system(makedbcmd)
        
    def run_BLAST(self):
        """Run BLAST with the query file and the modified reference fasta file"""
        blastn_cline = NcbiblastnCommandline(query = self.query, db= self.output_database.name, outfmt=6, out=self.blast_output.name, task ='blastn', evalue= self.evalue)
        print(tempfile.gettempprefix())
        stdout, stderr = blastn_cline()
