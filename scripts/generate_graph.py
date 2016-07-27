#!/usr/bin/env python3

import argparse
from madansi.GenerateGraph import GenerateGraph
from madansi.ValidateInputArguments import ValidateInputArguments
from madansi.ValidateOutputArguments import ValidateOutputArguments
from madansi.SwitchColumnsBlastFile import SwitchColumnsBlastFile

parser = argparse.ArgumentParser(description = 'Given a graph file and a file with data produces a linear subgraph')
parser.add_argument('inputgraph',     help ='path to the input graph file',    type=str)
parser.add_argument('inputdata', help='path to the input data file', type=str)
parser.add_argument('outputdata', help = 'path to the blast file with columns switched', type = str)
parser.add_argument('outputgraph',    help ='path to the output graph file',   type=str)
parser.add_argument('outputsequences', help = 'path to the output file containing unused sequences', type =str)
args = parser.parse_args()
#
#via = ValidateInputArguments(args.inputgraph, args.inputdata)
#via.run()
#
#voa = ValidateOutputArguments(args.output)
#voa.run()
#
#sw = SwitchColumnsBlastFile(args.inputdata, args.outputdata)
#sw.switch_columns_blast_file()

gg = GenerateGraph(args.inputgraph, args.outputdata, args.outputgraph, args.outputsequences)
gg.generate_graph()
