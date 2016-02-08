#!/usr/bin/env python3

"""Parse a HMMer tblout file line by line, storing relevant information. Optional output of coverage vectors for each
HMM with the skewness pipeline.
Note: Reads in the original fasta/fastq file MUST have a unique header for this to work properly.
"""


#############
## Imports ##
#############
import os.path
import argparse
import numpy as np
import sys
import re
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
class TbloutParser:
    """This object takes as input a HMMer tblout file path and constructs an iterable that outputs
    hash-mapping of header to sequence information.  Only one line will be held in memory at a time using this method.
    """

    def __init__(self, filepath, outpath, filename, length, evalue=10, multi=False, clstr_file=None,
                 annot_file=None, dupfile=None):
        """
        Constructor; can't touch this.  This is a hellish nightmare of an __init__ function.
        All of sound mind, turn back now.
        :param filepath: filepath to input hmmer tblout file
        :param outpath: path to output directory
        :param filename: basename for the output file (excluding extension)
        :param length: filepath to the file containing the hmm lengths
        :param evalue: Evalue threshold below which to keep a hit
        :param multi: optional flag to correct for multiple reads, maintaining a 1 to 1 ratio of reads to counts
        :param clstr_file: optional filepath to the .clstr file from cd-hit from hmm generation step
        :param annot_file: optional filepath to the hmmer_annotation file for each gene
        :param kmer: optional prefix to append to file basename if analyzing kmer outputs
        :param dupfile: optional filepath to the table of duplicate counts if used in previous pipeline steps
        :return: name, len, start, stop, str of filepath
        """
        if os.path.exists(filepath):  # if file is a file, read from the file
            self.hmmer_file = str(filepath)
            self.stdin = False
        elif not sys.stdin.isatty():  # else read from standard in
            self.stdin = True
        else:
            raise ValueError("Parameter filepath must be a HMMer tblout file")
        self.outpath = outpath
        self.filename = filename
        self.reads_mapping = 0
        self.reads_total = 0
        self.ethreshold = float(evalue)
        self.multicorrect = multi
        self.duplicate = dupfile
        self.duplicate_table = {}  # Counts of duplicate reads from kmer screen
        self.hmm_observed = {}  # Mapping of HMMs to hits observed
        self.class_observed = {}  # Mapping of Classes to hits observed
        self.mech_observed = {}  # Mapping of Mechanisms to hits observed
        self.group_observed = {}  # Mapping of Groups to hits observed
        self.clstr_members = {}  # Mapping of hmm # -> genes used to build that HMM
        self.hmm_lengths = {}  # Length of each hmm, hmm # -> length
        self.gene_multihits = {}  # Will store HMM hits for each gene
        with open(length, 'r') as hmm_length:
            data = hmm_length.read().split('\n')[1:]
            for line in data:
                line = line.split()
                if line:
                    self.hmm_lengths.setdefault(line[0], int(line[1]))
        with open(clstr_file, 'r') as f:
            ## Map HMM to its gene members, initialize the HMM two-by-two
            header_reg = re.compile(r'>(.+?)\.\.\.')
            line = f.readline()
            cluster_num = -1  # For zero indexing offset
            clstr = []
            while line:
                if line[0] is ">":
                    cluster_num += 1
                    if cluster_num > 0:  # Is this the first entry? If so, skip it
                        if len(clstr) is 1:  # No singletons
                            clstr = []
                        else:
                            self.clstr_members.setdefault(str(cluster_num), clstr)
                            clstr = []
                else:
                    clstr.append(header_reg.findall(line)[0])  # Begin appending headers before pushing to dict
                line = f.readline()
        self.gene_annots = {}  # Mapping from gene name -> annotations
        with open(annot_file, 'r') as annot:
            ## Map the gene names to their annotations, initialize the hierarchy two-by-twos
            data = annot.read().split('\n')
            for line in data:
                temp = line.split(',')
                if temp[0]:
                    self.gene_annots.setdefault(temp[0], temp[1:])
        self.hmm_annots = {}  # HMM annotation mapping from key (hmm #) -> annotations
        for key, values in self.clstr_members.items():
            ## Calculate the annotations of each HMM, combining around a pipe '|' if multiple
            classes, mechs, groups = zip(*[self.gene_annots[x] for x in values])
            classes = [x for x in classes if x]
            mechs = [x for x in mechs if x]
            groups = [x for x in groups if x]
            self.hmm_annots.setdefault(key, ['|'.join(set(classes)), '|'.join(set(mechs)), '|'.join(set(groups))])
        if dupfile:
            with open(dupfile, 'r') as indup:
                ## Put each gene header into a key and its counts into the value.  Initialize the obs count dict
                data = indup.read().split('\n')
                for line in data:
                    if line:
                        temp = line.split()
                        if temp and int(temp[0]) > 1:
                            self.duplicate_table.setdefault(temp[1], int(temp[0]))

    def __iter__(self):
        return self

    @property
    def _iterate(self):
        # Skip all leading whitespace
        while True:
            if self.stdin:
                hmmer_line = sys.stdin.readline()  # read from stdin
            else:
                hmmer_line = self.hmmer_file.readline()  # read from file
            if not hmmer_line:
                return  # End of file
            if hmmer_line.startswith("#"):  # these lines contain header information
                continue
            else:  # all other lines are hits
                self.reads_total += 1
                if self.reads_total % 100000 == 0:  # update the counter on stdout every 100000 reads
                    sys.stdout.write("\rHits processed: {}".format(self.reads_total))
                    sys.stdout.flush()
                temp = hmmer_line.split()
                read_name = temp[0]
                if float(temp[12]) <= float(self.ethreshold):
                    self.reads_mapping += 1
                    ## Basic observation increment rules
                    if self.duplicate:
                        try:
                            multiplier = self.duplicate_table[read_name]
                        except KeyError:
                            multiplier = 1
                        try:
                            self.hmm_observed[temp[2]] += multiplier
                        except KeyError:
                            self.hmm_observed.setdefault(temp[2], multiplier)
                        try:
                            self.gene_multihits[read_name][temp[2]] += multiplier
                        except KeyError:
                            try:
                                self.gene_multihits[read_name].setdefault(temp[2], multiplier)
                            except KeyError:
                                self.gene_multihits.setdefault(read_name, {temp[2]: multiplier})
                    else:
                        try:
                            self.hmm_observed[temp[2]] += 1
                        except KeyError:
                            self.hmm_observed.setdefault(temp[2], 1)
                        try:
                            self.gene_multihits[read_name][temp[2]] += 1
                        except KeyError:
                            try:
                                self.gene_multihits[read_name].setdefault(temp[2], 1)
                            except KeyError:
                                self.gene_multihits.setdefault(read_name, {temp[2]: 1})
                    return temp[2], self.hmm_lengths[temp[2]], temp[4], temp[5], temp[11]  # name, len, start, stop, str
        self.hmmer_file.close()  # catch all in case this line is reached
        assert False, 'Should not reach this line'

    def correct_multihit(self):
        """
        For our purposes, reads need to maintain a one-to-one mapping with hits counted.  Therefore, we both need to
        correct for hits across multiple HMMs due to homology and also for different motifs within the same read
        hitting multiple times within the same HMM. Note that for this to work, reads must have a unique header in the original fastq/fasta file
        that is at the end of the read, separated by a pipe '|'.

        Correct the observed counts for each hierarchy and the HMMs, since we can't
        a priori know what is true and what is not.  THIS CORRECTION SHOULD NOT BE USED FOR ASSEMBLED CONTIGS, since
        we can't know from assembled contigs whether they should truly multi-map or not.  It is advised to only use
        this feature for raw reads.
        :return: void
        """
        if self.multicorrect:
            for key, subdict in self.gene_multihits.items():
                if subdict:
                    ## Correct for HMM counts
                    if len(subdict) > 1:
                        for nkey, nvalue in subdict.items():
                            self.hmm_observed[nkey] -= nvalue
                            self.hmm_observed[nkey] += float(nvalue) / sum(subdict.values())
                    if len(subdict) == 1:
                        for nkey, nvalue in subdict.items():
                            self.hmm_observed[nkey] -= nvalue
                            self.hmm_observed[nkey] += 1

    def aggregate_hierarchy(self):
        """
        Create the actual aggregated hierarchy counts based on the annotation file.  All counts are fundamentally
        aggregated from the HMM assignments.  We are interested in the accuracy of HMM assignment and classification,
        not necessarily in the absolute annotation of the reads themselves.  By default, reads are split evenly if more
        than one annotation is present in the HMM when aggregating up.  The user is left to decide whether to round the
        output values or leave them as floats.
        :return: void
        """
        for key, value in self.hmm_observed.items():
            annots = self.hmm_annots[key]
            class_annot = annots[0].split('|')
            mech_annot = annots[1].split('|')
            group_annot = annots[2].split('|')
            for entry in class_annot:
                try:
                    self.class_observed[entry] += value / float(len(class_annot))
                except KeyError:
                    self.class_observed.setdefault(entry, value / float(len(class_annot)))
            for entry in mech_annot:
                try:
                    self.mech_observed[entry] += value / float(len(mech_annot))
                except KeyError:
                    self.mech_observed.setdefault(entry, value / float(len(mech_annot)))
            for entry in group_annot:
                try:
                    self.group_observed[entry] += value / float(len(group_annot))
                except KeyError:
                    self.group_observed.setdefault(entry, value / float(len(group_annot)))

    def write_observed(self):
        """
        Output the counts (or corrected counts) to the output file based on the hierarchy of annotations.
        :return: void
        """
        pathname = self.outpath + '/' + self.filename + '.csv'
        with open(pathname, 'w') as observed_file:
            observed_file.write('Hierarchy,Name,Abundance,E-value\n')
            for key, values in self.hmm_observed.items():
                observed_file.write('HMM,' + key + ',' + str(values) + ',' + str(self.ethreshold) + '\n')
            for key, values in self.class_observed.items():
                observed_file.write('Class,' + key + ',' + str(values) + ',' + str(self.ethreshold) + '\n')
            for key, values in self.mech_observed.items():
                observed_file.write('Mechanism,' + key + ',' + str(values) + ',' + str(self.ethreshold) + '\n')
            for key, values in self.group_observed.items():
                observed_file.write('Group,' + key + ',' + str(values) + ',' + str(self.ethreshold) + '\n')

    def __next__(self):
        if not self.stdin and type(self.hmmer_file) is str:  # only open file here if hmmer_file is a str and not fileIO
            self.hmmer_file = open(self.hmmer_file, 'r')
        value = self._iterate
        if not value:  # close file on EOF
            if not self.stdin:
                self.hmmer_file.close()
            if self.multicorrect:
                self.correct_multihit()
            self.aggregate_hierarchy()
            self.write_observed()
            raise StopIteration()
        else:
            return value


##############
## ArgParse ##
##############
parser = argparse.ArgumentParser('hmmer_parse.py')
parser.add_argument('input', type=str, help='File path to HMMer tblout file, or "-" for stdin')
parser.add_argument('output', type=str, help='Base directory path for desired outputs')
parser.add_argument('filename', type=str, help='File name for this particular file (.csv format')
parser.add_argument('hmm_len', type=str, help='Path to file containing HMM lengths')
parser.add_argument('clstr', type=str, help='Path to file containing clstr generation info')
parser.add_argument('annots', type=str, help='Path to annotation file')
parser.add_argument('-d', '--dupfile', type=str, default=None,
                    help='Path to duplicate count file if used in the kmer screen')
parser.add_argument('-e', '--evalue', type=float, default=1e-14, help='Evalue under which to keep hits')
parser.add_argument('-m', '--multicorrect', action='store_true', default=False,
                    help='If set, reads have a one to one mapping with reported hits')
parser.add_argument('-s', '--skewfile', type=str, default=None, help='Optional output file for HMM skewness metrics')
parser.add_argument('-g', '--graphs', type=str, help='Path to output directory for graphs')


##########
## Main ##
##########
if __name__ == '__main__':
    ## Parse the arguments using ArgParse
    args = parser.parse_args()
    infile = args.input
    outpath = args.output
    if args.skewfile:
        if args.graphs:
            graph_dir = args.graphs
            if not os.path.isdir(graph_dir):
                os.mkdir(graph_dir)
        else:
            sys.stdout.write('\nGraph directory not defined, using output directory by default\n')
            graph_dir = outpath
    if not os.path.isdir(outpath):
        os.mkdir(outpath)

    ## Determine alignment positions and overlay the reads onto the gene vectors.
    ## Calculate coverage and measures of skewness.
    ## Output a master file describing each mapped vector/metric and a coverage graph for each mapped gene.
    vector_hash = {}  # This stores gene names that were referenced in the SAM file, along with their vectors
    vector_counts = {}  # This stores gene names as in the other dictionary, but stores read counts instead
    for line in TbloutParser(infile, outpath, args.filename, args.hmm_len, args.evalue, args.multicorrect,
                             args.clstr, args.annots, args.dupfile):
        if args.skewfile:
            if not line:
                continue
            if int(line[2]) < int(line[3]):
                start = line[2]
                stop = line[3]
            else:
                start = line[3]
                stop = line[2]
            try:
                vector_hash[line[0]][(int(start) - 1):(int(stop))] += 1  # increment affected region
            except KeyError:
                vector_hash[line[0]] = np.zeros(int(line[1])).astype('int')  # new entry, initialize vector
            ## If skewness metrics are to be written, calculate the entropy measures and output coverage plots.
            with open(args.skewfile, 'w') as out:
                out.write('Accession_Name,Accession_Length,Coverage,Shannon_Entropy,L2norm_Deviation,Vector\n')  # headers for outfile
                plt.ioff()  # no interactive mode
                plt.hold(False)  # don't keep plot
                for key, value in vector_hash.items():
                    if sum(value > 0):  # if the gene had reads aligned
                        plot_count += 1
                        sys.stdout.write("\rPlots generated: {}".format(plot_count))  # update counter for user benefit
                        sys.stdout.flush()

                        ## Calculate metrics
                        coverage = float(sum(value > 0)) / len(value)  # what percentage of the gene has a read aligned?
                        norm_vec = value ** float(1) / sum(value)  # normalize so the vector sums to 1 (frequency)
                        max_entropy = np.ones(len(norm_vec)) / len(norm_vec)  # this is used to calculate the maximum shannon entropy for a given vector
                        shannon = np.negative(sum([x * np.log2(x) for x in norm_vec if x > 0])) / np.negative(
                                sum([x * np.log2(x) for x in max_entropy]))  # % of max possible Shannon entropy
                        l2norm = 1 - ((np.sqrt(sum(norm_vec * norm_vec)) * np.sqrt(len(norm_vec))) - 1) / (
                            np.sqrt(len(norm_vec)) - 1)  # % of max possible conformity to the L2 norm unit sphere
                        out.write(",".join([key, str(len(value)), str(coverage), str(shannon), str(l2norm),
                                            " ".join([str(x) for x in value])]) + '\n')

                        ## Plot figure
                        plt.plot(value)
                        plt.xlabel('Nucleotide Position')
                        plt.ylabel('Observed Count')
                        plt.title('Coverage Plot for {}'.format(key))
                        plt.savefig(graph_dir + '/' + '{}.png'.format(key))
                        plt.close()  # make sure plot is closed
