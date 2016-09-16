import unittest
from madansi.RunBLAST import RunBLAST
import filecmp
import os
import shutil
import tempfile

class TestRunBLAST(unittest.TestCase):
    
    def test_initialise_object_no_evalue(self):
        """Tests the initialisation if evalue not specified"""
        temp_dir = tempfile.mkdtemp()
        rb = RunBLAST('ab', 'cd', temp_dir)
        self.assertTrue(rb.query)
        self.assertTrue(rb.input_reference)
        self.assertTrue(rb.output_reference)
        self.assertTrue(rb.output_database)
        self.assertTrue(rb.blast_output)
        self.assertEqual(rb.evalue, 0.01)
        shutil.rmtree(temp_dir)
    
    def test_initialise_object_evalue(self):
        """Tests if correctly initialised with evalue specified"""
        temp_dir = tempfile.mkdtemp()
        rb = RunBLAST('ab','cd', temp_dir, evalue=0.0001 )
        self.assertTrue(rb.query)
        self.assertTrue(rb.input_reference)
        self.assertTrue(rb.output_reference)
        self.assertTrue(rb.output_database)
        self.assertTrue(rb.blast_output)
        self.assertEqual(rb.evalue, 0.0001)
        shutil.rmtree(temp_dir)
    
    def test_switch_columns(self):
        """Tests if the columns are correctly swapped in the first two columns"""
        temp_dir = tempfile.mkdtemp()
        rb = RunBLAST('madansi/tests/data/input_query.fa','madansi/tests/data/input_reference.fa', temp_dir)
        rb.run_switch_columns_database()
        self.assertTrue(filecmp.cmp(rb.output_reference.name,'madansi/tests/data/expected_output_reference.fa', shallow=False))
        shutil.rmtree(temp_dir)
      
    def test_reference_database(self):
        """Tests that the expected database is created"""
        temp_dir = tempfile.mkdtemp()
        rb = RunBLAST('madansi/tests/data/input_query.fa','madansi/tests/data/input_reference.fa', temp_dir)
        rb.run_switch_columns_database()
        rb.make_reference_database()
        self.assertTrue(os.path.isfile(rb.output_database.name+'.nin'))
        self.assertTrue(os.path.isfile(rb.output_database.name+'.nhr'))
        self.assertTrue(os.path.isfile(rb.output_database.name+'.nsq'))
        self.assertTrue(filecmp.cmp(rb.output_database.name+'.nsq', 'madansi/tests/data/expected_output_database.fa.nsq', shallow=False)) 
        self.assertTrue(filecmp.cmp(rb.output_database.name+'.nhr', 'madansi/tests/data/expected_output_database.fa.nhr', shallow=False))
        shutil.rmtree(temp_dir)
    
    def test_run_BLAST(self):
        """Tests that the output from running BLAST is as expected"""
        temp_dir = tempfile.mkdtemp()
        rb = RunBLAST('madansi/tests/data/input_query.fa','madansi/tests/data/input_reference.fa',temp_dir)
        rb.run_switch_columns_database()
        rb.make_reference_database()
        rb.run_BLAST()
        self.assertTrue(filecmp.cmp(rb.blast_output.name, 'madansi/tests/data/expected_comparison', shallow=False))
        shutil.rmtree(temp_dir)
    
        
        