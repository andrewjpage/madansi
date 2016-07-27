#!/usr/bin/env python3

import argparse
from madansi.RunBLAST import RunBLAST
from madansi.ValidateInputArguments import ValidateInputArguments
from madansi.ValidateOutputArguments import ValidateOutputArguments

parser = argparse.ArgumentParser(description = 'Runs BLAST on a database which has had switch_columns applied to it with a query fasta file')
parser.add_argument('query', help = 'path to the query file', type =str)
parser.add_argument('inputreference', help='path to the input reference fasta file', type=str)
parser.add_argument('outputreference', help= 'path to the output reference fasta file with switched columns', type=str)
parser.add_argument('outputdatabase', help= 'path to the output database', type = str)
parser.add_argument('blastoutput', help= 'path to where the output file will be outputted', type = str)
args = parser.parse_args()


#via = ValidateInputArguments(args.query)
#via.run()
#
#voa = ValidateOutputArguments(args.outputdatabase, args.outputblast, args.finaloutput)
#voa.run()

rb = RunBLAST(args.query, args.inputreference, args.outputreference, args.outputdatabase, args.blastoutput)
rb.run_switch_columns_database()
rb.make_reference_database()
rb.run_BLAST()
