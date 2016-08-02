#!/usr/bin/env python3

import argparse
from madansi.ContigSearching import ContigSearching
from madansi.GeneDetector import GeneDetector
from madansi.GraphParser import GraphParser

parser = argparse.ArgumentParser(description = '')
parser.add_argument('input_assembly_file')
parser.add_argument('blast_hits_file')
parser.add_argument('graph_file')
args = parser.parse_args()


graph_parser = GraphParser(args.graph_file)
filtered_graph = graph_parser.graph
print('filtered graph')

gene_detector = GeneDetector(args.input_assembly_file, args.blast_hits_file)
gene_detector.contigs_to_genes()
print('Constructed contigs dictionary')

contig_searching = ContigSearching(gene_detector, filtered_graph)
contig_searching.expand_all_contigs()
print('contig_searching.neighbouring_contigs')