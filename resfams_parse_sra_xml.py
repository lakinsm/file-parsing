#!/usr/bin/env python

## Parse the SRA XML metadata from the Resfams testsets, output a mapping from sample name to truth label

import argparse
from lxml import objectify


parser = argparse.ArgumentParser('hmmer_parse.py')
parser.add_argument('input', type=str, help='File path to SRA XML metadata file')
parser.add_argument('output', type=str, help='File path to output csv')
parser.add_argument('set', type=str, help='Which test set is this?')
args = parser.parse_args()

if args.set == 'pediatric':
    with open(args.input, 'r') as infile, open(args.output, 'w') as out:
        tree = objectify.parse(infile)
        root = tree.getroot()
        for i in root.getchildren():
            temp = str(i.getchildren()[0].getchildren()[1]).split(' ')
            library, truth = (temp[-1], temp[0])
            sample = i.getchildren()[6].getchildren()[0].attrib['accession']
            out.write(sample + ',' + library + ',' + truth + '\n')
