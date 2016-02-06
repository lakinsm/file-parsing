#!/usr/bin/env python3

""" Hash the k-mers of a given database (AMR in our case) into a set.  Screen the reads in a fastq file against
the database hash table using each k-mer in a sliding window along the read.  If any k-mer hits the hash table,
then the read passes filter; forward and backward orientations are hashed.

Credit for the compression algorithm code blocks goes to Cathal Garvey (https://github.com/cathalgarvey/dncode.git)
"""


#############
## Imports ##
#############
import argparse
import pickle
import multiprocessing as mp
import sys
import logging
import resource
import itertools
import collections


##########
## Vars ##
##########
db_hash = set()  # object for the hash table
uniq_hash = collections.Counter()  # store unique reads and counts
chunksize = 20000000  # limit memory consumption by reading in blocks
window = 20  # k-mer size
overall = 0  # counter for stderr writing
acgt_encoding_table = {}  # DNA encoding table
acgt_decoding_table = {}  # DNA decoding table
iupac = 'RYSWKMBDHVN'  # the DNA encoding algorithm only handles standard base codes


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
    barray = ()
    for read_name, seq in chunk:
        global window
        for i in range(len(seq) - window + 1):
            subseq = seq[i:i + window]
            if subseq in db_hash:
                if args.unique:
                    if any([x in seq for x in iupac]):
                        logging.info('>' + seq + '\n' + seq)
                        break
                    bits = encode(seq)
                    barray += (bits, )
                else:
                    logging.info('>' + seq + '\n' + seq)
                break
    return collections.Counter(barray)


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
    k, m = int(len(a) / n), len(a) % n
    return (a[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in range(n))


def current_mem_usage():
    """
    :return: current memory usage in MB
    """
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024


def get_encbyte(mask_len, rna=False, gzipped=False):
    """
    By Cathal Garvey (https://github.com/cathalgarvey/dncode.git)
    The final byte of an encoded sequence is of form (binary):
    |unused|unused|unused|unused|gzipped|rna|mask|mask|
    ..that is, the length of the final n-quad mask is obtained using only
    the final two bits of the mask, and the third-last bit is 1 if RNA, 0 else.
    The remaining 5 bits are unspecified at present.
    """
    if not 0 <= mask_len <= 3:
        raise ValueError("Mask length can only be 0-3.")
    mask  = 0
    mask |= mask_len
    if rna:
        mask |= 1<<2
    if gzipped:
        mask |= 1<<3
    return mask.to_bytes(1, 'big')


def encode(code, rna=False):
    """
    By Cathal Garvey (https://github.com/cathalgarvey/dncode.git)
    Outputs a bytestring consisting of the encoded DNA, plus a mask byte.
    The encoding is static 4:1 compression, with quadrets of DNA translating to
    single-byte values evenly.
    The last encoded byte may be padded to permit encoding. To account for this
    a "mask" byte is always appended to the sequence, which consists simply of
    a big-endian integer specifying which of the last 4n are padding to be
    discarded on decoding. There is extra "room" in this last byte for other
    formatting details, such as alphabet, as a future improvement.
    """
    # Iterate over the Dut DNA in four-letter chunks, maybe plus trailing.
    return b''.join(iter_encode(code, rna=rna))


def iter_encode(code_iter, table = acgt_encoding_table, rna=False):
    """
    By Cathal Garvey (https://github.com/cathalgarvey/dncode.git)
    Iterates over sequences, preventing large sequence files being loaded into RAM.
    Yields encoded sequence bytes plus final encoding byte.
    """
    # Create four iterators; one for each of four frames to iterate in steps of four.
    # Then zip these frames and use map to combine them into quadlets. Using
    # filter in the arguments to 'join' prevents a TypeError on last block if any
    # NoneType values are yielded from the zip_longest iterator.
    fr1, fr2, fr3, fr4 = [itertools.islice(code_iter, i, None, 4) for i in range(4)]
    framezip = itertools.zip_longest(fr1,fr2,fr3,fr4)
    zip_quads = map(lambda t: ''.join(filter(None,t)), framezip)
    seql = 0
    for quad in zip_quads:
        try:
            enc_dna = acgt_encoding_table[quad]
            seql += 4
        except KeyError as E:
            # Should only evaluate on last quad, making all preceding lookups faster.
            # TODO: Make at least a token effort to ensure it's not a real KeyError.
            if len(quad) < 4:
                enc_dna = acgt_encoding_table[ quad + ( "A" * (4 - len(quad)) ) ]
                seql += len(quad)
            else:
                raise E
        yield enc_dna
    else:
        # Generate and yield final byte.
        mask = (4 - len(quad)) % 4  # Can only be 0-3, leaving unused byte-space..
        yield get_encbyte(mask, rna)


def decode_encbyte(encbyte):
    """
    By Cathal Garvey (https://github.com/cathalgarvey/dncode.git)
    Decodes a byte as encoded by get_encbyte.
    """
    #encbyte = int.from_byte(encbyte, 'big') # Be explicit on endianness
    if isinstance(encbyte, bytes):
        if not len(encbyte) == 1:
            raise ValueError("encbyte MAY be bytes, but must be length 1!")
        else:
            encbyte = int.from_bytes(encbyte, "big")
    if (not isinstance(encbyte, int)) or encbyte > 255:
        raise TypeError("decode_encbyte only accepts a single byte or a byte-size int. "
                        "Value given was of type {} and value {}.".format(type(encbyte), encbyte))
    return dict(
        mask_len=encbyte & 3,
        rna=True if encbyte & 1 << 2 else False,
        gzipped=True if encbyte & 1 << 3 else False )


def iter_decode(enc_dna, encoding, table=acgt_decoding_table):
    """
    By Cathal Garvey (https://github.com/cathalgarvey/dncode.git)
    Returns an iterator to decode the input sequence or iterator over a sequence.
    """
    table = table.copy()
    if encoding['rna']:
        for k in table:
            table[k] = table[k].replace("T", "U")
    for b in enc_dna:
        yield table[b.to_bytes(1, 'big')]


def decode_all(enccode):
    """
    By Cathal Garvey (https://github.com/cathalgarvey/dncode.git)
    Straight decode of a sequence. Does not support iterators.
    """
    return ''.join(iter_decode(enccode[:-1], decode_encbyte(enccode[-1])))



##############
## ArgParse ##
##############
parser = argparse.ArgumentParser('amr_skewness.py')
parser.add_argument('-d', '--database', type=str, default=None, help='File path to fasta database if new hash')
parser.add_argument('-p', '--pickle', type=str, default=None, help='Optional flag: database is a pickled hash table')
parser.add_argument('-s', '--save', type=str, default=None, help='Optional: save the hash table to a pickle file')
parser.add_argument('-n', '--num_process', type=int, default=1, help='Number of processes to run in parallel')
parser.add_argument('-k', '--kmer', type=int, default=15, help='K-mer size')
parser.add_argument('-u', '--unique', type=str, default=None,
                    help='File to store hashes of unique read counts for HMMER (use for highly redundant fastq files)')


##########
## Main ##
##########
if __name__ == '__main__':
    mp.freeze_support()
    ## Parse the arguments using ArgParse
    args = parser.parse_args()
    window = args.kmer

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
            for i in range(len(seq) - window + 1):
                temp = seq[i:i + window]
                db_hash.add(temp)
                db_hash.add(temp[::-1])

    ## Pickle the file if save option set
    if args.save:
        pickle.dump(db_hash, open(args.save, 'wb'))

    ## If --unique is set, then create the DNA encoding table for compression of unique hashes
    ## See file header for credit to Cathal Garvey for the compression algorithm
    # For every four-letter combination of ACGT (256 possibilities):
    if args.unique:
        for n, n4 in enumerate((''.join(x) for x in itertools.product(*('ACGT',) * 4))):
            nb = n.to_bytes(1, 'big')
            acgt_encoding_table[n4] = nb
            acgt_decoding_table[nb] = n4

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
        if args.unique:
            uniq_hash += sum((x for x in res), collections.Counter())
        handler.flush()  # flush the logging cache to stdout
        sys.stderr.write('\nFinished block.  Loading next chunk\n')
        del res
        del chunks  # remove chunks from memory.  Otherwise memory usage will be doubled.
        if check < chunksize:
            pool.close()  # ask nicely
            pool.join()  # sigterm
            pool.terminate()  # sigkill
            del pool  # make sure pool is cleared
            break
    sys.stderr.write('\nTotal reads processed {}'.format(overall))

    ## Write to duplicate mapping file if --unique flag is set
    if args.unique:
        sys.stderr.write('\nWriting duplicate hash values\n')
        with open(args.unique, 'wb') as dupfile:
            for key, value in uniq_hash.items():
                if value > 1:
                    dupfile.write(key + b'\@@' + value.to_bytes(4, 'big') + b'\..')
                seq = decode_all(key)
                sys.stdout.write('>' + seq + '\n' + seq + '\n')

