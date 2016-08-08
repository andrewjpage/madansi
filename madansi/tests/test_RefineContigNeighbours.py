import unittest
from madansi.RefineContigNeighbours import RefineContigNeighbours

class TestRefineContigNeighbours(unittest.TestCase):
    def test_count_contig_appearances_empty_list_and_dictionary(self):
        """Tests when both the gene list and the gene to contig dictionary is empty"""
        refine_contig_neighbours = RefineContigNeighbours('','','madansi/tests/data/empty_file')
        refine_contig_neighbours.genes = {}
        self.assertEqual(refine_contig_neighbours.count_contig_appearance([]), {})
    
    def test_count_contig_appearance_dictionary_empty(self):
        """Tests when there is a gene in the gene list not given in the dictionary"""
        refine_contig_neighbours = RefineContigNeighbours('','','madansi/tests/data/empty_file')
        refine_contig_neighbours.genes = {}
        self.assertEqual(refine_contig_neighbours.count_contig_appearance(['gene1']), {})
    
    def test_count_contig_appearance_gene_list_empty(self):
        """Tests when the dictionary is not empty but the gene list is"""
        refine_contig_neighbours = RefineContigNeighbours('','','madansi/tests/data/empty_file')
        refine_contig_neighbours.genes = {'gene1':'Contig1'}
        self.assertEqual(refine_contig_neighbours.count_contig_appearance([]), {})
        
    def test_count_contig_appearance_one_entry(self):
        """Tests when the gene list has one entry present in the dictionary"""
        refine_contig_neighbours = RefineContigNeighbours('','','madansi/tests/data/empty_file')
        refine_contig_neighbours.genes = {'gene1':'Contig1', 'gene2':'Contig1', 'gene3':'Contig2'}
        self.assertEqual(refine_contig_neighbours.count_contig_appearance(['gene1']), {'Contig1':1})
    
    def test_count_contig_appearance_three_entries(self):
        """Tests when there are multiple genes in the same contig present in the list"""
        refine_contig_neighbours = RefineContigNeighbours('','','madansi/tests/data/empty_file')
        refine_contig_neighbours.genes = {'gene1':'Contig1', 'gene2':'Contig1', 'gene3':'Contig2'}
        self.assertEqual(refine_contig_neighbours.count_contig_appearance(['gene1', 'gene2', 'gene3']), {'Contig1':2, 'Contig2':1})
        
        #Need tests for refine neighbouring contigs!
        