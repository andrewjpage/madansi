import unittest
from madansi.BlastHit import BlastHit, Error, file_reader
import filecmp
import os

class TestParseBlast(unittest.TestCase):
    
    def test_init_BlastHit(self):
        '''Check that init followed by str gives the same as the original input'''

        testlines = ['Contig1\t1.1.B1.36.cap3_contig\t98.04\t51\t1\t0\t359\t409\t325\t375\t1e-22\t93.7',
                     'Contig1\t1.1.B1.36.cap3_contig\t98.04\t51\t1\t0\t359\t409\t325\t375\t1e-22\t93',
                     'Contig1\t1.1.B1.36.cap3_contig\t98.04\t51\t1\t0\t359\t409\t325\t375\t0.0\t93.7',
                     'Contig1\t1.1.B1.36.cap3_contig\t98.04\t51\t1\t0\t359\t409\t325\t375\t0.01\t93.7',
                     'Contig1\t1.1.B1.36.cap3_contig\t98.04\t51\t1\t0\t359\t409\t325\t375\t0.0\t93.7\t100\t200']

        for t in testlines:
            self.assertEqual(t, str(BlastHit(t)))
            print(BlastHit(t).qry_name)
            
        bad_testlines = ['Contig1\t1.1.B1.36.cap3_contig\tX\t51\t1\t0\t359\t409\t325\t375\t1e-22\t93.7',
                         'Contig1\t1.1.B1.36.cap3_contig\t99\t51.1\t1\t0\t359\t409\t325\t375\t1e-22\t93.7',
                         'Contig1\t1.1.B1.36.cap3_contig\t99\t51\tX\t0\t359\t409\t325\t375\t1e-22\t93.7',
                         'Contig1\t1.1.B1.36.cap3_contig\t99\t51\t1\tA\t359\t409\t325\t375\t1e-22\t93.7',
                         'Contig1\t1.1.B1.36.cap3_contig\t98.04\t51\t1\t0\t359\t409\t325\t375\t0.0\t93.7\t100']

        for l in bad_testlines:
            with self.assertRaises(Error):
                self.assertEqual(l, str(BlastHit(l)))
        
class TestFileReader(unittest.TestCase):
    def test_file_reader(self):
        '''file_reader should iterate through a blast file correctly'''
        tmp_out = 'blast_unittest.m8.tmp'

        for f in ['madansi/tests/data/blast_unittest.m8', 'madansi/tests/data/blast_unittest.m8.with_lengths']:
            blast_reader = file_reader('madansi/tests/data/blast_unittest.m8')
            fout = open(tmp_out, 'w')
            for hit in blast_reader:
                print(hit, file=fout)
            fout.close()
            self.assertTrue(filecmp.cmp('madansi/tests/data/blast_unittest.m8', tmp_out))
            os.unlink(tmp_out)
            
   # def test_blast_array(self):
   #     '''file_reader should create the correct array'''
   #     tmp_out = 'blast_unittest.m8.tmp'
   #
   #     for f in ['madansi/tests/data/blast_unittest.m8']: #, 'madansi/tests/data/blast_unittest.m8.with_lengths'
   #         blast_reader = file_reader('madansi/tests/data/blast_unittest.m8')
   #         fout = open(tmp_out, 'w')
   #         for hit in blast_reader:
   #             print(hit, file=fout)
   #         fout.close()
   #         expected_blast_array_m8 = [['Contig2','7.23.B265.14.cap3_contig','99.25','402','3','0','3','404','402','1','0.0','731'], \
   #         ['Contig2','7.23.B265.9.cap3_contig	99.50','202','1','0','3','204','402','201','4e-100','351'], \
   #         ['Contig3','7.23.B265.9.cap3_contig','99.57','231','1','0','1','231','1','231','2e-117','408'], \
   #         ['7.23.B265.14.cap3_contig','Contig2','99.25','402','3','0','1','402','404','3','0.0','731'], \
   #         ['7.23.B265.14.cap3_contig','7.23.B265.9.cap3_contig','98.51','202','3','0','201','402','201','402','3e-95','335'], \
   #         ['7.23.B265.9.cap3_contig','Contig3','99.57','231','1','0','1','231','1','231','2e-117','408'], \
   #         ['7.23.B265.9.cap3_contig','Contig2','99.50','202','1','0','201','402','204','3','4e-100','351'],\
   #         ['7.23.B265.9.cap3_contig','7.23.B265.14.cap3_contig','98.51','202','3','0','201','402','201','402','3e-95','335']]
   #         self.assertEqual(expected_blast_array_m8, blast_array)
   #         os.unlink(tmp_out)
        