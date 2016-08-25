import os
import shutil

class KeepTempFiles(object):
    def __init__(self, output_directory, output_reference, output_database, blast_output, filtered_blast_output,  output_fasta_file ):
        self.output_directory = output_directory
        self.output_reference = output_reference
        self.output_database = output_database
        self.blast_output = blast_output
        self.filtered_blast_output = filtered_blast_output
        self.output_fasta_file = output_fasta_file
        
    def create_new_directory(self):
        if self.output_directory != None:
            os.mkdir(self.output_directory)
            self.move_and_rename_reference()
            self.move_and_rename_database()
            self.move_and_rename_blast_output()
            self.move_and_rename_filtered_blast_output()
            self.move_and_rename_output_fasta_file()

    def move_and_rename_reference(self):
        base_name = os.path.basename(self.output_reference)
        os.rename(self.output_reference, self.output_directory+'/output_reference.fa')

    def move_and_rename_database(self):
        base_name = os.path.basename(self.output_database)
        os.rename(self.output_database+'.nhr', self.output_directory+'/output_database.nhr')
        os.rename(self.output_database+'.nsq', self.output_directory+'/output_database.nsq')
        os.rename(self.output_database+'.nin', self.output_directory+'/output_database.nin')
        os.rename(self.output_database, self.output_directory+'/output_database')
    
    def move_and_rename_blast_output(self):
        base_name = os.path.basename(self.blast_output)
        os.rename(self.blast_output, self.output_directory+'/blast_output')
        
    def move_and_rename_filtered_blast_output(self):
        base_name = os.path.basename(self.filtered_blast_output)
        os.rename(self.filtered_blast_output, self.output_directory+'/filtered_blast_output')    
    
    def move_and_rename_output_fasta_file(self):
        base_name = os.path.basename(self.output_fasta_file)
        os.rename(self.output_fasta_file, self.output_directory+'/'+ base_name)
        