#!/usr/bin/env python3

import argparse
from madansi.SwitchColumns import SwitchColumns
from madansi.ValidateInputArguments import ValidateInputArguments
from madansi.ValidateOutputArguments import ValidateOutputArguments

parser = argparse.ArgumentParser(description = 'Swaps the first two entries in a FASTA file')
parser.add_argument('input',     help ='path to the input file',    type=str)
parser.add_argument('output',    help ='path to the output file',   type=str)
args = parser.parse_args()

via = ValidateInputArguments(args.input)
via.run()

voa = ValidateOutputArguments(args.output)
voa.run()

sw = SwitchColumns(args.input, args.output)
sw.run()


