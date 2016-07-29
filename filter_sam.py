#!/usr/bin/env python3

import sys
import argparse


class SamParser(object):
    """
    Parse a SAM file and do something.  Leaving this class sparse for now,
    but it is meant to serve as a template for later.
    """

    def __init__(self, infile, len_threshold=None):
        self.sam_file = open(infile, 'r') if sys.stdin.isatty() else None
        self.len_threshold = len_threshold if len_threshold else None

    def parse_sam(self):
        if self.sam_file:
            line = self.sam_file.readline()
        else:
            line = sys.stdin.readline()
        while line:
            if line.startswith('@'):
                if self.sam_file:
                    line = self.sam_file.readline()
                else:
                    line = sys.stdin.readline()
                continue
            line = line.split('\t')
            if int(line[1]) & 4 != 0:
                if self.len_threshold and len(line[9]) >= self.len_threshold:
                    sys.stdout.write('>'+line[0]+'\n')
                    sys.stdout.write(line[9]+'\n')
            if self.sam_file:
                line = self.sam_file.readline()
            else:
                line = sys.stdin.readline()


parser = argparse.ArgumentParser()
parser.add_argument('sam_file', type=str, help='Input SAM file or - for stdin')
parser.add_argument('-l', '--len', type=int, default=None, help='Optional, sequences must be greater than this length')


if __name__ == '__main__':
    args = parser.parse_args()
    samp = SamParser(args.sam_file, args.len)
    try:
        samp.parse_sam()
    except IOError:
        pass





