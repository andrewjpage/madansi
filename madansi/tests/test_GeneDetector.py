import unittest
from madansi.GeneDetector import GeneDetector
from madansi.Contig import Contig
from madansi.BlastHit import BlastHit
from madansi.Gene import Gene

class TestGeneDetector(unittest.TestCase):
       
    def test_parse_blast_fits_empty_file(self):
        """Parses empty file correctly"""
        gene_detector= GeneDetector('madansi/tests/data/empty_file.fa', 'madansi/tests/data/empty_file')
        self.assertEqual(gene_detector.parse_blast_hits(), [])
    
    def test_parse_blast_hits_one(self):
        """Tests that the correct type of object is obtained from parsing one blast hit"""
        gene_detector = GeneDetector('madansi/tests/data/assembly.fa', 'madansi/tests/data/one_blast_hit')
        self.assertTrue(isinstance(gene_detector.parse_blast_hits()[0], BlastHit)) 
       
    def test_contigs_to_genes_no_hits(self):
        """Tests that the correct keys are given in the contigs object and that there is no value for each contig in this object"""
        gene_detector = GeneDetector('madansi/tests/data/assembly.fa', 'madansi/tests/data/empty_file')
        self.assertCountEqual(list(gene_detector.contigs_to_genes().keys()), ['Contig1', 'Contig2', 'Contig3'])
        self.assertEqual(gene_detector.contigs_to_genes()['Contig1'].gene_objects, {})        
       
    def test_contigs_to_genes_one_hit(self):
        """Tests that given one blast hit, will give the correct keys and values"""
        gene_detector = GeneDetector('madansi/tests/data/assembly.fa', 'madansi/tests/data/one_blast_hit')
        self.assertCountEqual(list(gene_detector.contigs_to_genes().keys()), ['Contig1', 'Contig2', 'Contig3'])
        
        self.assertEqual(gene_detector.contigs_to_genes()['Contig2'].gene_objects, {})
        self.assertEqual(gene_detector.contigs_to_genes()['Contig3'].gene_objects, {})
        self.assertCountEqual(list(gene_detector.contigs_to_genes()['Contig1'].gene_objects.keys()), ['gene1'])
    
        my_gene = Gene(-1,402,1,None, 'Contig1')
        self.assertEqual(gene_detector.contigs_to_genes()['Contig1'].gene_objects['gene1'].orientation, my_gene.orientation)
        self.assertEqual(gene_detector.contigs_to_genes()['Contig1'].gene_objects['gene1'].start, my_gene.start)
        self.assertEqual(gene_detector.contigs_to_genes()['Contig1'].gene_objects['gene1'].end, my_gene.end)
        self.assertEqual(gene_detector.contigs_to_genes()['Contig1'].gene_objects['gene1'].node, my_gene.node)
        self.assertEqual(gene_detector.contigs_to_genes()['Contig1'].gene_objects['gene1'].contig, my_gene.contig)
        
    def test_contigs_to_genes_four_hits(self):
        """Tests output for multiple blast hits"""
        gene_detector = GeneDetector('madansi/tests/data/assembly.fa', 'madansi/tests/data/four_blast_hits' )
        self.assertCountEqual(list(gene_detector.contigs_to_genes().keys()), ['Contig1', 'Contig2', 'Contig3'])
        self.assertCountEqual(list(gene_detector.contigs_to_genes()['Contig1'].gene_objects.keys()), ['gene1', 'gene2', 'gene3'])
        self.assertCountEqual(list(gene_detector.contigs_to_genes()['Contig2'].gene_objects.keys()), ['geneA'])
        self.assertCountEqual(list(gene_detector.contigs_to_genes()['Contig3'].gene_objects.keys()), [])
    

        
        
    