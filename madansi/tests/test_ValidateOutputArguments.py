import unittest
from madansi.ValidateOutputArguments import ValidateOutputArguments

class TestValidateOutputArguments(unittest.TestCase):
    
    def test_validate_output_file(self):
        voa = ValidateOutputArguments('madansi/tests/data/expected_output_reference.fa')
        self.assertRaises(ValueError, voa.run())
        voa = ValidateOutputArguments('madansi/tests/data/output_reference.fa')
        self.assertRaises(ValueError, voa.run())
                            
