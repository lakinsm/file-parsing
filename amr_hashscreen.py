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
import argparse
import cPickle as pickle
import multiprocessing as mp
import sys
import logging


##########
## Vars ##
##########
db_hash = set()
chunksize = 50000
window = 20


#############
## Methods ##
#############
def fastq_parse():
    for x in range(chunksize):
        line = sys.stdin.readline()
        if line.startswith("@"):
            read_name = line.rstrip()[1:]
            seq = sys.stdin.readline().rstrip()
            sys.stdin.readline()
            sys.stdin.readline()
        if not line:
            return  # stop iteration
        yield read_name, seq


def worker(chunk):
    global db_hash
    for read_name, seq in chunk:
        global window
        for i in range(len(seq)-window):
            subseq = seq[i:i+window]
            if subseq in db_hash or subseq[::-1] in db_hash:
                logging.info('>'+read_name+'\n'+seq)
                break


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


def split(a, n):
    """
    :param a: list of arbitrary length
    :param n: number of groups to split into
    :return: generator of chunks
    """
    k, m = len(a) / n, len(a) % n
    return (a[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in xrange(n))


##############
## ArgParse ##
##############
parser = argparse.ArgumentParser('amr_skewness.py')
parser.add_argument('input', type=str, help='File path to AMR SAM file, or "-" for stdin')
parser.add_argument('database', type=str, help='File path to fasta database')
parser.add_argument('-p', '--pickle', nargs='?', default=None, help='Optional flag: database is a pickled hash table')
parser.add_argument('-s', '--save', nargs='?', default=None, help='Optional: save the hash table to a pickle file')
parser.add_argument('-n', '--num_process', type=int, default=1, help='Number of processes to run in parallel')


##########
## Main ##
##########
if __name__ == '__main__':
    ## Setup the logger
    root = logging.getLogger()
    root.setLevel(logging.DEBUG)
    handler = logging.StreamHandler(stream=sys.stdout)
    handler.setLevel(logging.DEBUG)
    root.addHandler(handler)
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
            for i in range(len(seq)-window):
                db_hash.add(seq[i:i+window])
    if args.save:
        pickle.dump(db_hash, open(args.save, 'wb'))
    ## Read in each fastq chunk and check for membership
    pool = mp.Pool()
    while True:
            chunks = [z for z in split([x for x in fastq_parse()], args.num_process)]
            check = sum([len(x) for x in chunks])
            if check is 0:
                break
            pool.map(worker, chunks)
            handler.flush()
            if check < chunksize:
                break

