#!/usr/bin/env python

## Note the following from the Kraken manual for parsing purposes:
##
# OUTPUT FORMAT
# =============
#
# Each sequence classified by Kraken results in a single line of
# output.  Output lines contain five tab-delimited fields; from
# left to right, they are:
#
#   1) "C"/"U": one letter code indicating that the sequence was
#      either classified or unclassified.
#   2) The sequence ID, obtained from the FASTA/FASTQ header.
#   3) The taxonomy ID Kraken used to label the sequence; this is
#      0 if the sequence is unclassified.
#   4) The length of the sequence in bp.
#   5) A space-delimited list indicating the LCA mapping of each k-mer
#      in the sequence.  For example, "562:13 561:4 A:31 0:1 562:3"
#      would indicate that:
#      - the first 13 k-mers mapped to taxonomy ID #562
#      - the next 4 k-mers mapped to taxonomy ID #561
#      - the next 31 k-mers contained an ambiguous nucleotide
#      - the next k-mer was not in the database
#      - the last 3 k-mers mapped to taxonomy ID #562

#############
## Imports ##
#############
import os.path
import argparse

##########
## Vars ##
##########


#############
## Methods ##
#############
class KrakenParser:
    """This object takes as input a Kraken raw output file path and constructs an iterable that outputs
    a hierarchically aggregated taxon count table.  Only one line will be held in memory at a time using this method.
    """
    def __init__(self, filepath):
        """
        constructor
        @param filepath: filepath to the input raw Kraken file.
        """
        if os.path.isfile(filepath):
            self.kraken_file = filepath
        else:
            raise ValueError("Parameter filepath must be a Kraken file")
        self.current_line = None
        self.classified = ""
        self.sequence_id = ""
        self.taxonomy_id = ""
        self.length = 0
        self.lca_mapping = []

    def __iter__(self):
        return self

    def _iterate(self):
        # Skip all leading whitespace
        while True:
            kraken_line = self.kraken_file.readline()
            if not kraken_line:
                return  # End of file
            if not kraken_line.isspace():
                break  # Break on finding something

        while True:
            self.classified, \
            self.sequence_id, \
            self.taxonomy_id, \
            self.length, \
            self.lca_mapping = kraken_line.rstrip().split('\t')
            self.lca_mapping = self.lca_mapping.split()
            ##
            ## Do something here to manipulate values before returning them
            ##
            return self.classified, self.sequence_id, self.taxonomy_id, self.length, self.lca_mapping
        self.kraken_file.close()
        assert False, "Should not reach this line"

    def next(self):
        if type(self.kraken_file) is str:
            self.kraken_file = open(self.kraken_file, "r")
        value = self._iterate()
        if not value:
            self.kraken_file.close()
            raise StopIteration()
        else:
            return value


##############
## Argparse ##
##############
parser = argparse.ArgumentParser('kraken_parse.py')
parser.add_argument('inputfile', type=str, help='File path to raw Kraken file')
parser.add_argument('outputfile', type=str, help='File path to desired output file')

##########
## Main ##
##########
if __name__ == '__main__':
    ## Parse the arguments using ArgParse
    args = parser.parse_args()
    infile = args.inputfile
    outfile = args.outputfile
    with open(outfile, 'w') as out:
        for classified, seq_id, tax_id, length, lca_mapping in KrakenParser(infile):
            ## Do something to each line of variables
