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
import resource


##########
## Vars ##
##########
db_hash = set()  # object for the hash table
chunksize = 20000000  # limit memory consumption by reading in blocks
window = 20  # k-mer size
overall = 0  # counter for stderr writing


#############
## Methods ##
#############
def worker(chunk):
    """ This code block is used in a MapReduce manner: the chunk of code is executed many times across
    the data block (each worker receives a chunk of that block).  The reads that pass filter are written to
    the logging cache (this is because writing to stdout produces thrashing).  The logging cache is then flushed
    on every iteration of the outer loop.
    :param chunk: a chunk of reads divided amongst the pool of parallel workers
    :return: void
    """
    global db_hash
    for read_name, seq in chunk:
        global window
        for i in range(len(seq)-window):
            subseq = seq[i:i+window]
            if subseq in db_hash:
                logging.info('>'+read_name+'\n'+seq)
                break


def fastq_parse():
    """ Parses a fastq file in chunks of 4 lines, skipping lines that don't inform fasta creation.
    This script only accepts stdin, so use cat file | script.py for correct functionality.
    :return: generator yielding tuples of (read_name, seq)
    """
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


def split(a, n):
    """ Splits an input list into n equal chunks; this works even if modulo > 0.
    :param a: list of arbitrary length
    :param n: number of groups to split into
    :return: generator of chunks
    """
    k, m = len(a) / n, len(a) % n
    return (a[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in xrange(n))


def current_mem_usage():
    """
    :return: current memory usage in MB
    """
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024


##############
## ArgParse ##
##############
parser = argparse.ArgumentParser('amr_skewness.py')
parser.add_argument('-d', '--database', nargs='?', default=None, help='File path to fasta database if new hash')
parser.add_argument('-p', '--pickle', nargs='?', default=None, help='Optional flag: database is a pickled hash table')
parser.add_argument('-s', '--save', nargs='?', default=None, help='Optional: save the hash table to a pickle file')
parser.add_argument('-n', '--num_process', type=int, default=1, help='Number of processes to run in parallel')


##########
## Main ##
##########
if __name__ == '__main__':
    mp.freeze_support()
    ## Input must be on stdin; raise error if this is not the case
    if sys.stdin.isatty():
        raise IOError('Input must be on stdin.  Use stream redirect for correct functionality: cat file | script.py')

    ## Setup the logger for output of reads to stdout.  This is necessary because writing directly to stdout
    ## in parallel causes thrashing and variable results.  The logger caches reads passed to it on every loop
    ## and flushes to stdout after the reads have been tested for membership.
    root = logging.getLogger()
    root.setLevel(logging.DEBUG)
    handler = logging.StreamHandler(stream=sys.stdout)
    handler.setLevel(logging.DEBUG)
    root.addHandler(handler)

    ## Parse the arguments using ArgParse
    args = parser.parse_args()

    ## Either pickle or database must be defined
    if not args.pickle and not args.database:
        raise IOError('Either --database or --pickle must be defined')

    ## Unpickle if the database hash is in a pickle file, or read in the fasta if performing a new hash
    ## Construct a hash table for the fasta database for each unique k-mer.
    ## Forward and reverse k-mers are added to the table.
    if args.pickle:
        db_hash = pickle.load(open(args.pickle, 'rb'))
    else:
        db_fasta = [x[1] for x in fasta_parse(args.database)]
        for seq in db_fasta:
            for i in range(len(seq)-window):
                temp = seq[i:i+window]
                db_hash.add(temp)
                db_hash.add(temp[::-1])

    ## Pickle the file if save option set
    if args.save:
        pickle.dump(db_hash, open(args.save, 'wb'))

    ## Read in each fastq chunk and check for membership.  Chunk size should be set such that
    ## the block size doesn't overflow memory.  Keep in mind this block size has the potential to be doubled
    ## in the logging cache.  Chunksize is set to 20 million reads for Bovine.
    while True:
        pool = mp.Pool(processes=args.num_process)  # create pool of workers for parallel processing
        chunks = [z for z in split([x for x in fastq_parse()], args.num_process)]  # divide reads into chunks
        sys.stderr.write('\nMemory used: {}MB'.format(current_mem_usage()))
        check = sum([len(x) for x in chunks])  # this is the break condition for the while loop (count of reads)
        overall += check  # add to overall read count for reporting to stderr
        sys.stderr.write('\nTotal reads processed {}'.format(overall))
        if check is 0:
            pool.close()
            pool.join()
            pool.terminate()
            del pool
            break
        res = pool.map(worker, chunks)  # pool.map is MapReduce.  All workers must finish before proceeding.
        handler.flush()  # flush the logging cache to stdout
        sys.stderr.write('\nFinished block.  Loading next chunk\n')
        del chunks  # remove chunks from memory.  Otherwise memory usage will be doubled.
        if check < chunksize:
            pool.close()  # ask nicely
            pool.join()  # sigterm
            pool.terminate()  # sigkill
            del pool  # make sure pool is cleared
            break

