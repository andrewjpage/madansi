import unittest
import os
import shutil
from madansi.KeepTempFiles import KeepTempFiles

class TestTempFile(unittest.TestCase):
    
    def test_move_and_rename_file(self):
        my_file = open('madansi/tests/data/assembly_4_sequences_ref.fa', 'w')
        my_file.close()
        ktf = KeepTempFiles('output_directory','madansi/tests/data/assembly_4_sequences_ref.fa', 'output_database', 'blast_output', 'filtered_blast_output',  'output_fasta_file ')
        ktf.create_new_directory()
        ktf.move_and_rename_reference()
        self.assertTrue(os.path.exists('output_directory/output_reference.fa'))
        shutil.rmtree('output_directory')