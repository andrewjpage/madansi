import unittest
from madansi.Contig import Contig
from madansi.BlastHit import BlastHit

class TestContig(unittest.TestCase):
    def test_initiailisation(self):
        contig = Contig('Sequence_name')
        self.assertEqual(contig.genes(), {})
    
    def test_add_one_blast_result(self):
        blast_result = 'Contig1\tgene1\t98.04\t51\t1\t0\t359\t409\t325\t375\t1e-22\t93.7'
        blast_hit_object = BlastHit(blast_result)
        contig = Contig('Contig1')
        self.assertTrue(contig.add_blast_hit(blast_hit_object))
        self.assertCountEqual(list(contig.genes().keys()), ['gene1'])
    
    def test_add_three_blast_results(self):
        blast_results = ['Contig1\tgene1\t98.04\t51\t1\t0\t359\t409\t325\t375\t1e-22\t93.7',
                     'Contig1\tgene2\t98.04\t51\t1\t0\t359\t409\t325\t375\t1e-22\t93',
                     'Contig1\tgene3\t98.04\t51\t1\t0\t359\t409\t325\t375\t0.0\t93.7']
        contig = Contig('Contig1')
        for blast_result in blast_results:
            blast_hit_object = BlastHit(blast_result)
            self.assertTrue(contig.add_blast_hit(blast_hit_object))
        self.assertCountEqual(list(contig.genes().keys()), ['gene1', 'gene2', 'gene3'])
    
    def test_two_genes_same_name_on_same_contig(self):
        blast_results = ['Contig1\tgene1\t98.04\t51\t1\t0\t359\t409\t325\t375\t1e-22\t93.7',
                     'Contig1\tgene1\t98.04\t51\t1\t0\t359\t409\t325\t375\t1e-22\t93']
        contig = Contig('Contig1')
        blast_hit_object = BlastHit(blast_results[0])
        self.assertTrue(contig.add_blast_hit(blast_hit_object))
        
        blast_hit_object = BlastHit(blast_results[1])
        self.assertFalse(contig.add_blast_hit(blast_hit_object))
        
        self.assertCountEqual(list(contig.genes().keys()), ['gene1'])
        