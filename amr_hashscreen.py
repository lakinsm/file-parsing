#!/usr/bin/env python

""" Hash the 20-mers of a given database (AMR in our case) into a set.  Screen the reads in a fastq file against
the database hash table using each 20-mer in a sliding window along the read.  If any 20-mer hits the hash table,
then the read passes filter; forward and backward orientations are checked.  All hashes are implemented in bit form:
A, C, G, T = { 00, 01, 10, 11 }
"""


## Fastq format for reference:
# @HWI-ST916:307:C6039ACXX:5:1101:1242:2102 1:N:0:CGATGT
# TTGAAGCTTGGTGAGTTTGCCCCAACAAGAACTTTCCGTGGGCATGGAGGAGCAAAATCAACCGCTCCGGCACCGAAGAAGTAGGAGGGACATATGGCTA
# +
# BBBFFFFFFFFFFIIFIIIIIIIIFIIIIIIIIIIIIIFFIIIIIIIIIIIIIIIIIIIIFFFFFFFFFFFFFFFBFFFFFBBFFFFFFFFBFFFFFFFB



#############
## Imports ##
#############
import os.path
import argparse
import cPickle as pickle
import numpy as np
import sys


##########
## Vars ##
##########
db_hash = set()


#############
## Methods ##
#############
class HashScreen:
    """This object takes as input a fastq file against the and returns each header and sequence only.
    Only one line will be held in memory at a time using this method.
    """
    def __init__(self, filepath):
        """
        constructor
        @param filepath: filepath to the input fastq file.
        """
        if os.path.exists(filepath):  # if file is a file, read from the file
            self.fq_file = str(filepath)
            self.stdin = False
        elif not sys.stdin.isatty():  # else read from standard in
            self.stdin = True
        else:
            raise ValueError("Parameter filepath must be a fastq file")
        self.current_line = None
        self.read_name = None
        self.seq = None
        #self.bit_table = {'A': 0b00, 'C': 0b01, 'G': 0b10, 'T': 0b11}
        self.offset = 0

    def __iter__(self):
        return self

    @property
    def _iterate(self):
        while True:
            self.offset += 1
            if self.stdin:
                fq_line = sys.stdin.readline()  # read from stdin
            else:
                fq_line = self.fq_file.readline()  # read from file
            if not fq_line:
                return  # End of file
            else:
                if self.offset % 100000 == 0:  # update the counter on stderr every 100000 reads
                    sys.stderr.write("\rReads screened: {}".format(self.offset))
                    sys.stderr.flush()
                if fq_line.startswith("@"):
                    self.read_name = fq_line[1:].rstrip()
                else:
                    self.seq = fq_line.rstrip()
                    if self.stdin:
                        sys.stdin.readline()
                        sys.stdin.readline()
                    else:
                        self.fq_file.readline()
                        self.fq_file.readline()
                    return self.read_name, self.seq
        self.fq_file.close()  # catch all in case this line is reached
        assert False, "Should not reach this line"

    def next(self):
        if not self.stdin and type(self.fq_file) is str:  # only open file here if fastq is a str and not fileIO
            self.fq_file = open(self.fq_file, "r")
        value = self._iterate
        if not value:  # close file on EOF
            if not self.stdin:
                self.fq_file.close()
            raise StopIteration()
        else:
            return value


def fasta_parse(infile):
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


##############
## ArgParse ##
##############
parser = argparse.ArgumentParser('amr_skewness.py')
parser.add_argument('input', type=str, help='File path to AMR SAM file, or "-" for stdin')
parser.add_argument('database', type=str, help='File path to fasta database')
parser.add_argument('-p', '--pickle', nargs='?', default=None, help='Optional flag: database is a pickled hash table')
parser.add_argument('-s', '--save', nargs='?', default=None, help='Optional: save the hash table to a pickle file')


##########
## Main ##
##########
if __name__ == '__main__':
    ## Parse the arguments using ArgParse
    args = parser.parse_args()
    infile = args.input
    if args.pickle:
        db_hash = pickle.load(open(args.pickle, 'rb'))
    else:
        db_fasta = [x[1] for x in fasta_parse(args.database)]
    ## Construct a hash table for the fasta database
    if not args.pickle:
        for seq in db_fasta:
            window = 20
            for i in range(len(seq)-window):
                db_hash.add(seq[i:i+window])
    if args.save:
        pickle.dump(db_hash, open(args.save, 'wb'))
    ## Read in each fastq line and check for membership
    for name, seq in HashScreen(infile):
        window = 20
        if not seq:
            break
        for i in range(len(seq)-window):
            if seq[i:i+window] in db_hash:
                sys.stdout.write('>'+name+'\n'+seq+'\n')
                break
    sys.stderr.write('\nDone.\n')
