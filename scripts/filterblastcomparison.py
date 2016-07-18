#!/usr/bin/env python3

import argparse
from madansi.FilterBlastComparison import FilterBlastComparison

parser = argparse.ArgumentParser(description = 'Filters the results from a BLAST comparison')
parser.add_argument('comparisonfile', help = 'Path to the input file that will be filtered', type=str)
parser.add_argument('filteredfile', help = 'Path to the filtered output file', type= str)
parser.add_argument('percentidentity', help= 'Parameter used to filter based on percent identity', type = float, default=0.0)
parser.add_argument('alignmentlength', help='Parameter used to filter based on the alignment length', type=int, default=0)
parser.add_argument('mismatches', help = 'Parameter used to filter based on the number of mismatches', type =int, default=10)
parser.add_argument('gapopenings', help='Parameter used to filter based on the number of gap openings', type=int, default=10)
parser.add_argument('evalue', help='Parameter used to filter based on the expect value', type=float, default=0.0)
parser.add_argument('bitscore', help='Parameter used to filter based on the bit score', type=float, default=0.0)
args = parser.parse_args()

via = ValidateInputArguments(args.comparisonfile)
via.run()

voa = ValidateOutputArguments(args.filteredfile)
voa.run

fbc = FilteredBlastComparison(args.comparisonfile, args.filteredfile, args.percentidentity, args.alignmentlength, args.mismatches, args.gapopenings, args.evalue, args.bitscore)
fbc.filter()