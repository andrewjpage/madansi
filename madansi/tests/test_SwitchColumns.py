import unittest
from madansi.SwitchColumns import SwitchColumns
import filecmp

class TestSwitchColumns(unittest.TestCase):
    
    def test_initialise_object(self):
        sw = SwitchColumns('abc','efg')
        self.assertTrue(sw.inputfile)
        self.assertTrue(sw.outputfile)
    
    def test_switching_around_columns(self):
        sw = SwitchColumns('madansi/tests/data/input_reference.fa','output_file')
        sw.run()
        self.assertTrue(filecmp.cmp('output_file','madansi/tests/data/expected_output_reference.fa'))