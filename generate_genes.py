#!/usr/bin/env python


import random
import argparse


def fasta_parse(infile):
    """ Parses a fasta file in chunks of 2 lines.
    :param infile: path to the input fasta file
    :return: generator of (header, sequence) fasta tuples
    """
    with open(infile, 'r') as fasta_file:
        # Skip whitespace
        while True:
            line = fasta_file.readline()
            if line is "":
                return  # Empty file or premature end of file?
            if line[0] is ">":
                break
        while True:
            if line[0] is not ">":
                raise ValueError("Records in FASTA should begin with '>'")
            header = line[1:].rstrip()
            all_lines = []
            line = fasta_file.readline()
            while True:
                if not line:
                    break
                if line[0] is ">":
                    break
                all_lines.append(line.rstrip())
                line = fasta_file.readline()
            yield header, "".join(all_lines).replace(" ", "").replace("\r", "")
            if not line:
                return  # Stop Iteration
        assert False, "Should not reach this line"


def generate_genes(infile, outfile):
    candidates = []
    with open(outfile, 'aw') as g:
        # genome_name = infile.split('/')[-1].split('.')[0]
        iteration = 1
        for header, line in fasta_parse(infile):
            gene_name = header+'_'+str(iteration)
            slide = 0
            length = len(line)
            while True:
                window = interval = random.randint(300, 4000)
                if slide+window <= length:
                    subset = line[slide:slide+window]
                    slide += interval
                    iteration += 1
                    if any((x in subset for x in ['N', 'R', 'Y', 'S', 'W', 'K', 'M', 'B', 'D', 'H', 'V', 'N', '.', '-'])):
                        continue
                    candidates.append('>'+gene_name+'_'+str(iteration)+'\n'+subset+'\n')
                if slide+window > length:
                    subset = line[slide:slide+window]
                    if any((x in subset for x in ['N', 'R', 'Y', 'S', 'W', 'K', 'M', 'B', 'D', 'H', 'V', 'N', '.', '-'])):
                        break
                    candidates.append('>'+gene_name+'_'+str(iteration)+'\n'+subset+'\n')
                    iteration = 1
                    break
        if len(candidates) < 100:
            chosen = [candidates[x] for x in random.sample(xrange(0, len(candidates)), len(candidates))]
        else:
            chosen = [candidates[x] for x in random.sample(xrange(0, len(candidates)), 100)]
        for entry in chosen:
            g.write(entry)


## Main
random.seed(154)
parser = argparse.ArgumentParser('generate_genes.py')
parser.add_argument('-i', '--infiles', nargs="+", help='File path to fasta genome file')
parser.add_argument('db_file', type=str, help='File path to AMR database fasta file')
parser.add_argument('outfile', type=str, help='File path to output fasta file')


args = parser.parse_args()

db_fasta = [x for x in fasta_parse(args.db_file)]
subset = [db_fasta[x] for x in random.sample(xrange(0, len(db_fasta)), 400)]
for genome in args.infiles:
    generate_genes(genome, args.outfile)
with open(args.outfile, 'aw') as g:
    for entry in subset:
        g.write('>'+entry[0]+'\n'+entry[1]+'\n')


