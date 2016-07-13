#!/usr/bin/env python3

import argparse
from madansi.RunBLAST import RunBLAST

parser = argparse.ArgumentParser(description = 'Runs BLAST on a database which has had switch_columns applied to it with a query fasta file')
parser.add_argument('reference', help = 'path to the reference fasta file', type=str)
parser.add_argument('outputreference', help = 'path to the output reference fasta file', type= str)
parser.add_argument('outputdatabase', help= 'path to the output database', type = str)
parser.add_argument('query', help = 'path to the query file', type =str)
parser.add_argument('output', help='path to the output file', type=str)
args = parser.parse_args()


via = ValidateInputArguments(args.reference, args.query)
via.run()

voa = ValidateOutputArguments(args.outputreference, args.outputdatabase, args.output)
voa.run

rb = RunBLAST(args.reference, args.outputreference, args.outputdatabase, args.query, args.output)
rb.run_switch_columns()
rb.make_reference_database()
rb.run_BLAST()
