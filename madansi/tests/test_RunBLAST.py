import unittest
from madansi.RunBLAST import RunBLAST
import filecmp
import os

class TestRunBLAST(unittest.TestCase):
    
    def test_initialise_object_no_evalue(self):
        """Tests the initialisation if evalue not specified"""
        rb = RunBLAST('ab', 'cd', 'ef', 'gh', 'ij')
        self.assertTrue(rb.query)
        self.assertTrue(rb.input_reference)
        self.assertTrue(rb.output_reference)
        self.assertTrue(rb.output_database)
        self.assertTrue(rb.blast_output)
        self.assertEqual(rb.evalue, 0.01)
    
    def test_initialise_object_evalue(self):
        """Tests if correctly initialised with evalue specified"""
        rb = RunBLAST('ab','cd', 'ef', 'gh', 'ij', evalue=0.0001 )
        self.assertTrue(rb.query)
        self.assertTrue(rb.input_reference)
        self.assertTrue(rb.output_reference)
        self.assertTrue(rb.output_database)
        self.assertTrue(rb.blast_output)
        self.assertEqual(rb.evalue, 0.0001)
   
    def test_switch_columns(self):
         """Tests if the columns are correctly swapped in the first two columns"""
         rb = RunBLAST('madansi/tests/data/input_query.fa','madansi/tests/data/input_reference.fa', 'output_reference' , 'output_database', 'output_file')
         rb.run_switch_columns_database()
         self.assertTrue(filecmp.cmp('output_reference','madansi/tests/data/expected_output_reference.fa', shallow=False))
         os.unlink('output_reference')
      
    def test_reference_database(self):
        """Tests that the expected database is created"""
        rb = RunBLAST('madansi/tests/data/input_query.fa','madansi/tests/data/input_reference.fa', 'madansi/tests/data/expected_output_reference.fa' , 'output_database', 'output_file')
        rb.make_reference_database()
        self.assertTrue(os.path.isfile('output_database.nin'))
        self.assertTrue(os.path.isfile('output_database.nhr'))
        self.assertTrue(os.path.isfile('output_database.nsq'))
        self.assertTrue(filecmp.cmp('output_database.nsq', 'madansi/tests/data/expected_output_database.fa.nsq', shallow=False)) 
        self.assertTrue(filecmp.cmp('output_database.nhr', 'madansi/tests/data/expected_output_database.fa.nhr', shallow=False))
        self.clean_up()
      
    def test_run_BLAST(self):
        """Tests that the output from running BLAST is as expected"""
        rb = RunBLAST('madansi/tests/data/input_query.fa','madansi/tests/data/input_reference.fa', 'madansi/tests/data/expected_output_reference.fa', 'madansi/tests/data/expected_output_database.fa',  'output_comparison')
        rb.run_BLAST()
        self.assertTrue(filecmp.cmp('output_comparison', 'madansi/tests/data/expected_comparison', shallow=False))
        self.clean_up()
        
    def clean_up(self):
        for filename in ['output_reference', 'output_database.nin', 'output_database.nhr', 'output_database.nsq', 'output_comparison']:
            if os.path.exists(filename):
                os.remove(filename)
        