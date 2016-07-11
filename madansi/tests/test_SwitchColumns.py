import unittest
from madansi.SwitchColumns import SwitchColumns
import filecmp

class TestSwitchColumns(unittest.TestCase):
    
    def test_initialise_object(self):
        sw = SwitchColumns('abc','efg')
        self.assertTrue(sw.inputfile)
        self.assertTrue(sw.outputfile)
    
    def test_switching_around_columns(self):
        sw = SwitchColumns('madansi/tests/data/inputfile.fa','outputfile')
        sw.run()
        self.assertTrue(filecmp.cmp('outputfile','madansi/tests/data/expected_output_file.fa'))