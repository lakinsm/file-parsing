#!/usr/bin/env python

## This script filters the output results from BLASTing the HMM singletons against the nt BLAST database.
## Matches only pass filter if they meet the inclusion criteria of cd-hit (80% identity at 50% overlap of a given gene).

#############
## Imports ##
#############
import argparse


##########
## Vars ##
##########
pass_filter = []


#############
## Methods ##
#############
def filter_hits(infile, outfile):
    with open(infile, 'r') as hitfile, open(outfile, 'w') as out:
        data = hitfile.read().strip().split()
        for hit in data:
            info = hit.split('\t')
            gene_id, target_acc, target_id, gene_len, gene_start, gene_stop, target_len, _, evalue, pident, \
                target_title, seq = info
            if int(pident) > 99:
                continue
            if int(target_len) < (int(gene_len) / 2):
                continue






##############
## ArgParse ##
##############
parser = argparse.ArgumentParser('amr_skewness.py')
parser.add_argument('infile', type=str, help='File path to custom BLAST tabular output file')
parser.add_argument('outfile', type=str, help='File path to output file')
args = parser.parse_args()


##########
## Main ##
##########
if __name__ == '__main__':
    filter_hits(args.infile, args.outfile)
