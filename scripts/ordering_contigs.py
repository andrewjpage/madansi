#!/usr/bin/env python3

import argparse
import os
from madansi.TemporaryDirectory             import TemporaryDirectory
from madansi.RunBLAST                       import RunBLAST
from madansi.ContigSearching                import ContigSearching
from madansi.GeneDetector                   import GeneDetector
from madansi.GraphParser                    import GraphParser
from madansi.FilterBlastComparison          import FilterBlastComparison
from madansi.UnusedContigs                  import UnusedContigs
from madansi.ProduceOrderedContigGraph      import ProduceOrderedContigGraph
from madansi.IterateJoiningContigComponents import IterateJoiningContigComponents
from madansi.GraphToFasta                   import GraphToFasta
from madansi.RemoveContigs                  import RemoveContigs  
from madansi.KeepTempFiles                  import KeepTempFiles 

parser = argparse.ArgumentParser(description =  'Script to take in the input assembly file, blast hits file and graph file and output a list of how each of the \
                                                contigs are connected to each other and the number of iterations it takes to get between the two.')
parser.add_argument('input_reference',              help ='path to the input reference fasta file',                                     type = str)
parser.add_argument('input_fasta_file',          help = 'Path to the input assembly file.',                                          type = str)
parser.add_argument('input_graph_file',             help = 'Path to the input graph file.',                                             type = str)

parser.add_argument('--output_fasta_file', '-f', help='Output filename', default = 'reordered_contigs.fa')
parser.add_argument('--output_directory', '-o', help='Output directory for all files produced')
args = parser.parse_args()

print(args.output_directory)

temp_directory = TemporaryDirectory()
temp_dir = temp_directory.create_temporary_directory()

rb = RunBLAST(args.input_fasta_file, args.input_reference, temp_dir)
rb.run_switch_columns_database()
rb.make_reference_database()
rb.run_BLAST()

fbc = FilterBlastComparison(rb.blast_output.name, temp_dir, bit_score=200)
fbc.filter()

graph_parser    = GraphParser(args.input_graph_file)
filtered_graph  = graph_parser.graph

gene_detector   = GeneDetector(args.input_fasta_file, fbc.filtered_blast_output.name)
gene_detector.contigs_to_genes()
gene_detector.assembly.sequence_names()
sequences = gene_detector.assembly.sequences

unused_contigs  = UnusedContigs(gene_detector, args.output_fasta_file, args.input_fasta_file )
unused_contigs.contigs_not_in_filtered_file()

produced_ordered_contig_graph   = ProduceOrderedContigGraph(gene_detector, filtered_graph, fbc.filtered_blast_output.name, sequences)
ordered_contig_graph            = produced_ordered_contig_graph.produce_ordered_contig_graph()
contig_ends                     = produced_ordered_contig_graph.contig_ends

remove_contigs = RemoveContigs(ordered_contig_graph)
ordered_contig_graph_filtered = remove_contigs.remove_extra_contigs()
output_filtered_graph = IterateJoiningContigComponents(ordered_contig_graph_filtered)

graph_to_fasta = GraphToFasta(sequences, ordered_contig_graph_filtered, args.output_fasta_file, contig_ends)
graph_to_fasta.create_fasta_file_combined_contigs()

unused_contigs.contigs_not_in_filtered_graph(ordered_contig_graph_filtered)
unused_contigs.add_unused_contigs_to_end()

KeepTempFiles(args.output_directory, rb.output_reference.name, rb.output_database.name, rb.blast_output.name, fbc.filtered_blast_output.name, args.output_fasta_file.name).create_new_directory()

temp_directory.remove_temporary_directory(temp_dir)




