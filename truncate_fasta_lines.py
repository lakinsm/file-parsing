#!/usr/bin/env python

## Imports
import argparse

parser = argparse.ArgumentParser('python subsetFasta.py')
parser.add_argument('fasta_file', type=str, help="Input FASTA File")
parser.add_argument('output', type=str, help="Output FASTA File")
parser.add_argument('-s', '--subset_file', type=str, help="Optional file of headers for subset")
parser.add_argument('-q', '--fastq', action='store_true', help='Is the file fastq?')

args = parser.parse_args()


def fastq_parse(infile):
    """ Parses a fastq file in chunks of 4 lines, skipping lines that don't inform fasta creation.
    This script only accepts stdin, so use cat file | script.py for correct functionality.
    :return: generator yielding tuples of (read_name, seq)
    """
    with open(infile, 'r') as fastq_file:
        while True:
            line = fastq_file.readline()
            if line.startswith("@"):
                read_name = line.rstrip()[1:]
                seq = fastq_file.readline().rstrip()
                line3 = fastq_file.readline().rstrip('\n')
                line4 = fastq_file.readline().rstrip('\n')
            if not line:
                return  # stop iteration
            yield read_name, seq, line3, line4


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

if args.fastq:
    data = [x for x in fastq_parse(args.fasta_file)]
else:
    data = [x for x in fastaParse(args.fasta_file)]
if args.subset_file:
    with open(args.subset_file, 'r') as f, open(args.output, 'w') as out:
        subset_headers = f.read().split("\n")
        for header,seq in data:
            if header in subset_headers:
                out.write(">"+header+"\n"+seq+"\n")
else:
    with open(args.output, 'w') as out:
        if args.fastq:
            for header, seq, line3, line4 in data:
                if len(seq) > 124:
                    seq = seq[0:124]
                if len(line4) > 124:
                    line4 = line4[0:124]
                out.write("@"+header+'\n'+seq+'\n'+line3+'\n'+line4+'\n')
        else:
            for header,seq in data:
                if len(seq) > 124:
                    seq = seq[0:124]
                out.write(">"+header+"\n"+seq+"\n")


