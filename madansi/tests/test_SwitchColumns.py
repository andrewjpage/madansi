import unittest
from madansi.SwitchColumns import SwitchColumns
import filecmp
import os

class TestSwitchColumns(unittest.TestCase):
    
    def test_initialise_object(self):
        sw = SwitchColumns('abc','efg')
        self.assertTrue(sw.input_file)
        self.assertTrue(sw.output_file)
    
    def test_switching_around_columns(self):
        sw = SwitchColumns('madansi/tests/data/input_reference.fa','output_file')
        sw.run()
        self.assertTrue(filecmp.cmp('output_file','madansi/tests/data/expected_output_reference.fa', shallow=False))
        self.clean_up()
    
    def clean_up(self):
        if os.path.exists('output_file'):
                os.remove('output_file')