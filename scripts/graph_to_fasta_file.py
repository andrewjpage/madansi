#!/usr/bin/env python3

from madansi.GraphToFasta import GraphToFasta
import networkx as nx
import argparse


parser = argparse.ArgumentParser()
parser.add_argument('input_fasta_file', type = str)
parser.add_argument('ordered_contig_graph', type = str)
parser.add_argument('output_fasta_file', type = str)
args = parser.parse_args()

graph = nx.drawing.nx_pydot.read_dot(args.ordered_contig_graph)
graph_to_fasta = GraphToFasta(args.input_fasta_file, graph, args.output_fasta_file)
graph_to_fasta.create_fasta_file()