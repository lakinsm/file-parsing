#!/usr/bin/env python3

## Imports
import argparse
import random
import re

parser = argparse.ArgumentParser('python subsetFasta.py')
parser.add_argument('fasta_file', type=str, help="Input FASTA File")
parser.add_argument('output', type=str, help="Output FASTA File")
parser.add_argument('-s', '--subset_file', type=str, help="Optional file of headers for subset")
parser.add_argument('-d', '--disjunction', action="store_true", default=None, help="Take the opposite of the subset")
parser.add_argument('-a', '--ambiguous', type=str, default='', help="Replace ambiguous nucleotides with random choice or N [Y/N]")
parser.add_argument('-n', '--screen', action="store_true", default=None, help="Screen the file for protein sequences and skip if true")
parser.add_argument('-u', '--dedup', action="store_true", default=None, help="Remove duplicated entries")

ambig = {'R':('A', 'G'),
         'Y':('C', 'T'),
         'S':('G', 'C'),
         'W':('A', 'T'),
         'K':('G', 'T'),
         'M':('A', 'C'),
         'B':('C', 'G', 'T'),
         'D':('A', 'G', 'T'),
         'H':('A', 'C', 'T'),
         'V':('A', 'C', 'G'),
         'N':('A', 'C', 'G', 'T')}
aa = re.compile(r'[EFILPQ]')
args = parser.parse_args()

def fastaParse(infile):
    with open(infile, 'r') as fastaFile:
        # Skip whitespace
        while True:
            line = fastaFile.readline()
            if line is "":
                return  # Empty file or premature end of file?
            if line[0] is ">":
                break
        while True:
            if line[0] is not ">":
                raise ValueError("Records in FASTA should begin with '>'")
            header = line[1:].rstrip()
            allLines = []
            line = fastaFile.readline()
            while True:
                if not line:
                    break
                if line[0] is ">":
                    break
                allLines.append(line.rstrip())
                line = fastaFile.readline()
            yield header, "".join(allLines).replace(" ", "").replace("\r", "")
            if not line:
                return  # Stop Iteration
        assert False, "Should not reach this line"

data = [x for x in fastaParse(args.fasta_file)]
if args.dedup:
    data = list(set(data))
if args.subset_file:
    with open(args.subset_file, 'r') as f, open(args.output, 'w') as out:
        subset_headers = f.read().split("\n")
        for header,seq in data:
            if args.ambiguous == 'N':
                seq = ''.join(['N' if x in ambig else x for x in seq])
            elif args.ambiguous:
                seq = ''.join([ambig[x][random.randit(0, len(ambig[x]) - 1)] if x in ambig else x for x in seq])
            if not args.disjunction:
                if header in subset_headers:
                    out.write(">"+header+"\n"+seq+"\n")
            else:
                if header not in subset_headers:
                    out.write(">"+header+"\n"+seq+"\n")
else:
    with open(args.output, 'w') as out:
        for header,seq in data:
            if aa.search(seq):
                continue
            if args.ambiguous == 'N':
                seq = ''.join(['N' if x in ambig.keys() else x for x in seq])
            elif args.ambiguous:
                seq = ''.join([ambig[x][random.randint(0, len(ambig[x]) - 1)] if x in ambig.keys() else x for x in seq])
            out.write(">"+header+"\n"+seq+"\n")


