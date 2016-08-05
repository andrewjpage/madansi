import unittest
from madansi.GenesToContig import GenesToContig

class TestGenesToContig(unittest.TestCase):
    
    def test_initialisation(self):
        gene = GenesToContig('')
        self.assertEqual(gene.genes, {})
    
    def test_empty_file(self):
        gene = GenesToContig('madansi/tests/data/empty_file')
        self.assertEqual(gene.genes_to_contig(), {})
    
    def test_one_blast_hit(self):
        gene = GenesToContig('madansi/tests/data/one_blast_hit')
        self.assertEqual(gene.genes_to_contig(), {'gene1':'Contig1'})
    
    def test_four_blast_hits(self):
        gene = GenesToContig('madansi/tests/data/four_blast_hits')
        self.assertEqual(gene.genes_to_contig(), {'gene1':'Contig1', 'gene2':'Contig1', 'gene3':'Contig1', 'geneA':'Contig2'})