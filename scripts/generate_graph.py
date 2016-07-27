#!/usr/bin/env python3

import argparse
from madansi.WalkGraphs import WalkGraphs
from madansi.ValidateInputArguments import ValidateInputArguments
from madansi.ValidateOutputArguments import ValidateOutputArguments
from madansi.SwitchColumnsBlastFile import SwitchColumnsBlastFile

parser = argparse.ArgumentParser(description = 'Given a graph file and a file with data produces a linear subgraph')
parser.add_argument('inputgraph',     help ='path to the input graph file',    type=str)
parser.add_argument('inputdata', help='path to the input data file', type=str)
parser.add_argument('outputdata', help = 'path to the blast file with columns switched', type = str)
parser.add_argument('output',    help ='path to the output file',   type=str)
args = parser.parse_args()

via = ValidateInputArguments(args.inputgraph, args.inputdata)
via.run()

voa = ValidateOutputArguments(args.output)
voa.run()

sw = SwitchColumnsBlastFile(args.inputdata, args.outputdata)
sw.switch_columns_blast_file()

wg = WalkGraphs(args.inputgraph, args.outputdata, args.output)
wg.create_linear_subgraph()
