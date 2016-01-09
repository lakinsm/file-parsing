#!/usr/bin/env python

"""Parse a SAM file line by line, storing relevant information.  Output a vector for each gene that allows
for subsequent calculation of coverage and skewness of the read distribution.
"""

# From the SAM file format specifications, for SAM line reference:
# Col Field   Type    Regexp/Range    Brief description
# 1   QNAME   String  [!-?A-~]{1,254} Query template NAME
# 2   FLAG    Int     [0,2^16-1]      bitwise FLAG
# 3   RNAME   String  \*|[!-()+-<>-~][!-~]* Reference sequence NAME
# 4   POS     Int     [0,2^31-1]    1-based leftmost mapping POSition
# 5   MAPQ    Int     [0,2^8-1]     MAPping Quality
# 6   CIGAR   String  \*|([0-9]+[MIDNSHPX=])+   CIGAR string
# 7   RNEXT   String  \*|=|[!-()+-<>-~][!-~]*   Ref.  name of the mate/next read
# 8   PNEXT   Int     [0,2^31-1]    Position of the mate/next read
# 9   TLEN    Int     [-2^31+1,2^31-1]  observed Template LENgth
# 10  SEQ     String  \*|[A-Za-z=.]+    segment SEQuence
# 11  QUAL  String    [!-~]+        ASCII of Phred-scaled base QUALity+33
#





#############
## Imports ##
#############
import os.path
import argparse
import numpy as np
import sys
import matplotlib as mpl  # load mpl to set the output device, then load pyplot as plt
mpl.use('Agg')  # No X server running on Bovine, use Agg for png generation instead
import matplotlib.pyplot as plt


##########
## Vars ##
##########
total_reads_mapped = 0
plot_count = 0



#############
## Methods ##
#############
class SamParser:
    """This object takes as input a SAM file path and constructs an iterable that outputs
    hash-mapping of header to sequence information.  Only one line will be held in memory at a time using this method.
    """
    def __init__(self, filepath):
        """
        constructor
        @param filepath: filepath to the input raw SAM file.
        """
        if os.path.exists(filepath):  # if file is a file, read from the file
            self.sam_file = str(filepath)
            self.stdin = False
        elif not sys.stdin.isatty():  # else read from standard in
            self.stdin = True
        else:
            raise ValueError("Parameter filepath must be a SAM file")
        self.current_line = None
        self.reads_mapping = 0
        self.reads_total = 0
        # Allowable bitflags for SAM file -> reads with both mates mapping, regardless of other flags
        self.true_flags = (99, 147, 83, 163, 67, 131, 115, 179, 81, 161, 97, 145, 65, 129, 113, 177)

    def __iter__(self):
        return self

    @property
    def _iterate(self):
        # Skip all leading whitespace
        while True:
            if self.stdin:
                sam_line = sys.stdin.readline()  # read from stdin
            else:
                sam_line = self.sam_file.readline()  # read from file
            if not sam_line:
                return  # End of file
            if sam_line.startswith("@SQ"):  # these lines contain refseq length information
                temp = sam_line.split()
                return temp[1][3:], temp[2][3:]
            elif sam_line[0] != "@":  # these lines are the actual reads
                self.reads_total += 1
                if self.reads_total % 100000 == 0:  # update the counter on stdout every 100000 reads
                    sys.stdout.write("\rReads processed: {}".format(self.reads_total))
                    sys.stdout.flush()
                temp = sam_line.split()
                if int(temp[1]) in self.true_flags and temp[2] is not "*" and int(temp[3]) is not 0:
                    self.reads_mapping += 1
                    return temp[1], temp[2], temp[3], temp[9]
        self.sam_file.close()  # catch all in case this line is reached
        assert False, "Should not reach this line"

    def next(self):
        if not self.stdin and type(self.sam_file) is str:  # only open file here if sam_file is a str and not fileIO
            self.sam_file = open(self.sam_file, "r")
        value = self._iterate
        if not value:  # close file on EOF
            if not self.stdin:
                self.sam_file.close()
            sys.stdout.write("\n{:d} with both mates mapped out of {:d} total reads\n".format(self.reads_mapping, self.reads_total))
            sys.stdout.flush()
            raise StopIteration()
        else:
            return value



##############
## ArgParse ##
##############
parser = argparse.ArgumentParser('amr_skewness.py')
parser.add_argument('input', type=str, help='File path to AMR SAM file, or "-" for stdin')
parser.add_argument('outputfile', type=str, help='File path to desired output file (.csv format)')
parser.add_argument('graph_dir', type=str, help='Path to output directory for graphs')

##########
## Main ##
##########
if __name__ == '__main__':
    ## Parse the arguments using ArgParse
    args = parser.parse_args()
    infile = args.input
    outfile = args.outputfile
    graph_dir = args.graph_dir
    if not os.path.isdir(graph_dir):
        os.mkdir(graph_dir)

    ## Determine alignment positions and overlay the reads onto the gene vectors.
    ## Calculate coverage and measures of skewness.
    ## Output a master file describing each mapped vector/metric and a coverage graph for each mapped gene.
    with open(outfile, 'w') as out:
        vector_hash = {}  # This stores gene names that were referenced in the SAM file, along with their vectors
        for line in SamParser(infile):
            if len(line) is 2:
                vector_hash[line[0]] = np.zeros(int(line[1])).astype('int')  # new entry, initialize vector
            else:
                vector_hash[line[1]][(int(line[2])-1):(int(line[2])-1)+len(line[3])] += 1  # increment affected region
        out.write('Accession_Name,Accession_Length,Coverage,Shannon_Entropy,L2norm_Deviation,Vector\n')  # headers for outfile
        plt.ioff()  # no interactive mode
        plt.hold(False)  # don't keep plot
        for key, value in vector_hash.iteritems():
            if sum(value > 0):  # if the gene had reads aligned
                plot_count += 1
                sys.stdout.write("\rPlots generated: {}".format(plot_count))  # update counter for user benefit
                sys.stdout.flush()

                ## Calculate metrics
                coverage = float(sum(value > 0)) / len(value)  # what percentage of the gene has a read aligned?
                norm_vec = value**float(1) / sum(value)  # normalize so the vector sums to 1 (frequency)
                shannon = np.negative(sum([x*np.log2(x) for x in norm_vec if x > 0]))  # Shannon entropy
                l2norm = ((np.sqrt(sum(norm_vec*norm_vec))*np.sqrt(len(norm_vec))) - 1) / (np.sqrt(len(norm_vec)) - 1)  # Deviation from the L2 norm unit sphere
                out.write(",".join([key, str(len(value)), str(coverage), str(shannon), str(l2norm), " ".join([str(x) for x in value])])+'\n')

                ## Plot figure
                plt.plot(value)
                plt.xlabel('Nucleotide Position')
                plt.ylabel('Observed Count')
                plt.title('Coverage Plot for {}'.format(key))
                plt.savefig(graph_dir+'/'+'{}.png'.format(key))
                plt.close()  # make sure plot is closed
