import unittest
from madansi.SwitchColumns import SwitchColumns
import filecmp
import tempfile
import os

class TestSwitchColumns(unittest.TestCase):
    
    def test_initialise_object(self):
        sw = SwitchColumns('abc','efg')
        self.assertTrue(sw.input_file)
        self.assertTrue(sw.output_file)
    
    def test_switching_around_columns(self):
        output_file = tempfile.NamedTemporaryFile(delete = False)
        sw = SwitchColumns('madansi/tests/data/input_reference.fa',output_file.name)
        sw.run()
        self.assertTrue(filecmp.cmp(output_file.name,'madansi/tests/data/expected_output_reference.fa', shallow=False))
        os.unlink(output_file.name)