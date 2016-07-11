#!/usr/bin/env python3

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import argparse
import os.path
import sys

parser = argparse.ArgumentParser(description = 'Swaps the first two entries in a FASTA file')
parser.add_argument('input',     help ='path to the input file',    type=str)
parser.add_argument('output',    help ='path to the output file',   type=str)
args = parser.parse_args()

my_records = []
with open(args.input, "r") as handle:
    for seq_record in SeqIO.parse( handle, "fasta"):
        desc = str(seq_record.description)
        n = desc.index(" ")
        Name = desc[:n]
        group = desc[n+1:]
        rec = SeqRecord(seq_record.seq, id = group, description = group + " " + Name )
        my_records.append(rec)
    handle.close()
            
SeqIO.write(my_records, args.output, "fasta")
             