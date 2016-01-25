#!/usr/bin/env python

"""Parse a HMMer tblout file line by line, storing relevant information.  Output a vector for each gene that allows
for subsequent calculation of coverage and skewness of the read distribution for each HMM.  If a truthset is provided,
calculate two-by-two accuracies. Optionally output a ROC graph and statistics flat file.
Note: Reads in the original fasta/fastq file MUST have a unique header for this to work properly.
"""

# Header values for reference:
# 1 target name
# 2 accession
# 3 query name
# 4 accession
# 5 hmmfrom
# 6 hmm to
# 7 alifrom
# 8 ali to
# 9 envfrom
# 10 env to
# 11 sq len
# 12 strand
# 13 E-value
# 14 score
# 15 bias
# 16 description of target (multi-value)


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
class HmmerTime:
    """This object takes as input a HMMer tblout file path and constructs an iterable that outputs
    hash-mapping of header to sequence information.  Only one line will be held in memory at a time using this method.
    The object walks along the file and, if a truthset is provided, outputs two-by-two values for accuracy.
    """

    def __init__(self, filepath, outpath, filename, length, evalue=10, multi=False, truthset=None, clstr_file=None,
                 annot_file=None, kmer=None):
        """
        Constructor; can't touch this.  This is a hellish nightmare of an __init__ function.
        All of sound mind, turn back now.
        :param filepath: filepath to input hmmer tblout file
        :param evalue: Evalue threshold below which to keep a hit
        :param truthset: optional filepath to counts from each unique truth accession (format is output of uniq -c)
        :param clstr_file: optional filepath to the .clstr file from cd-hit from hmm generation step
        :param annot_file: optional filepath to the hmmer_annotation file for each gene
        :param hmm_annot: optional filepath mapping each hmm_query name (this is a number) to its annotations
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
        self.kmer = kmer
        self.reads_mapping = 0
        self.reads_total = 0
        self.ethreshold = float(evalue)
        self.multicorrect = multi
        self.hmm_observed = {}  # Mapping of HMMs to hits observed
        self.class_observed = {}  # Mapping of Classes to hits observed
        self.mech_observed = {}  # Mapping of Mechanisms to hits observed
        self.group_observed = {}  # Mapping of Groups to hits observed
        self.clstr_members = {}  # Mapping of hmm # -> genes used to build that HMM
        self.truthset = False
        self.hmm_lengths = {}  # Length of each hmm, hmm # -> length
        self.gene_multihits = {}  # Will store HMM hits for each gene
        self.hmm_twobytwo = {}  # HMM matrix for ROC generation, only used for truthset
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
                            self.hmm_twobytwo.setdefault(str(cluster_num), [0, 0, 0, 0])  # TP, FP, FN, TN
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
        for key, values in self.clstr_members.iteritems():
            ## Calculate the annotations of each HMM, combining around a pipe '|' if multiple
            classes, mechs, groups = zip(*[self.gene_annots[x] for x in values])
            classes = [x for x in classes if x]
            mechs = [x for x in mechs if x]
            groups = [x for x in groups if x]
            self.hmm_annots.setdefault(key, ['|'.join(set(classes)), '|'.join(set(mechs)), '|'.join(set(groups))])
        if truthset:
            ## This section is complicated.  To calculate two-by-twos, we need to keep track of several things:
            ## 1. The true number of reads that map to a given gene from the test set
            ## 2. The number of reads (gene fragments) mapping to each HMM and which ones they are (their gene header)
            ## 3. The hierarchy of annotations for each gene (Class -> Mechanism -> Group -> HMM)
            ## 4. The cluster membership of each HMM (which genes truly belong to each HMM?)
            ## 5. The HMM membership of each gene (to which HMMs did the gene's hits map?)
            ## 6. If there are more than 1 in #5, then we need to split the reads or choose the best candidate
            ##
            ## Additionally, we need to find the truth values for the hierarchy of annotations and then aggregate the
            ## counts for each HMM into Group -> Mechanism -> Class.  The 1-to-1 mapping of reads is essential here
            ## so that the value of the false negative counts is not a negative number.
            self.truthset = True
            self.truthset_counts = {}  # True counts for each gene
            with open(truthset, 'r') as truth:
                ## Put each gene header into a key and its counts into the value.  Initialize the obs count dict
                data = truth.read().split('\n')
                for line in data:
                    if line:
                        temp = line.split()
                        temp = [temp[0], " ".join(temp[1:])]
                        if temp:
                            self.truthset_counts.setdefault(temp[1], int(temp[0]))
            self.group_twobytwo = {}  # Group matrix for ROC generation
            self.mech_twobytwo = {}  # Mechanism matrix for ROC generation
            self.class_twobytwo = {}  # Class matrix for ROC generation
            self.hmm_truth = {}  # True counts for each HMM
            self.class_truth = {}  # True counts for each class
            self.mech_truth = {}  # True counts for each mechanism
            self.group_truth = {}  # True counts for each group
            for key, values in self.clstr_members.iteritems():
                ## Generate true counts for each HMM
                total = sum([int(self.truthset_counts[x]) for x in values if x in self.truthset_counts])
                self.hmm_truth.setdefault(key, total)
            for key, values in self.hmm_truth.iteritems():
                ## Generate aggregated hierarchy truth counts from the HMM truth counts
                annots = self.hmm_annots[key]
                class_annot = annots[0].split('|')
                mech_annot = annots[1].split('|')
                group_annot = annots[2].split('|')
                for entry in class_annot:
                    try:
                        list1 = np.array(self.class_truth[entry])
                        list2 = np.array(values)
                        self.class_truth[entry] = list1 + (list2 / float(len(class_annot)))
                    except KeyError:
                        self.class_truth.setdefault(entry, np.array(values) / float(len(class_annot)))
                for entry in mech_annot:
                    try:
                        list1 = np.array(self.mech_truth[entry])
                        list2 = np.array(values)
                        self.mech_truth[entry] = list1 + (list2 / float(len(mech_annot)))
                    except KeyError:
                        self.mech_truth.setdefault(entry, np.array(values) / float(len(mech_annot)))
                for entry in group_annot:
                    try:
                        list1 = np.array(self.group_truth[entry])
                        list2 = np.array(values)
                        self.group_truth[entry] = list1 + (list2 / float(len(group_annot)))
                    except KeyError:
                        self.group_truth.setdefault(entry, np.array(values) / float(len(group_annot)))

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
                temp[0] = '|'.join(temp[0].split('|')[:-1])  # This is a predetermined format correction for testset
                if float(temp[12]) < float(self.ethreshold):
                    self.reads_mapping += 1
                    ## Basic observation increment rules
                    if self.truthset:
                        ## TRUTHSET ENABLED: What category of HMM-level two-by-two does this hit fall under?
                        ## Can only calculate TP/FP here; the others are done at the end
                        if temp[0] in self.clstr_members[temp[2]]:
                            self.hmm_twobytwo[temp[2]][0] += 1
                        else:
                            self.hmm_twobytwo[temp[2]][1] += 1
                        ## TRUTHSET ENABLED: Keep track of whether a gene hits across multiple HMMs
                        try:
                            self.gene_multihits[read_name][temp[2]] += 1
                        except KeyError:
                            try:
                                self.gene_multihits[read_name].setdefault(temp[2], 1)
                            except KeyError:
                                self.gene_multihits.setdefault(read_name, {temp[2]: 1})
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
        hitting multiple times within the same HMM.  This correction will be optional in the production version of
        the script.  Note that for this to work, reads must have a unique header in the original fastq/fasta file
        that is at the end of the read, separated by a pipe '|'.

        For non-truthset corrections, just correct the observed counts for each hierarchy and the HMMs, since we can't
        a priori know what is true and what is not.  THIS CORRECTION SHOULD NOT BE USED FOR ASSEMBLED CONTIGS, since
        we can't know from assembled contigs whether they should truly multi-map or not.  It is advised to only use
        this feature for raw reads.
        :return: void
        """
        if self.multicorrect and self.truthset:
            for key, subdict in self.gene_multihits.iteritems():
                key = '|'.join(key.split('|')[0:-1])
                if subdict:
                    ## Correct for HMM counts
                    if len(subdict) > 1:
                        for nkey, nvalue in subdict.iteritems():
                            if key in self.clstr_members[nkey]:
                                self.hmm_twobytwo[nkey][0] -= nvalue
                                self.hmm_twobytwo[nkey][0] += float(nvalue) / sum(subdict.values())
                            else:
                                self.hmm_twobytwo[nkey][1] -= nvalue
                                self.hmm_twobytwo[nkey][1] += float(nvalue) / sum(subdict.values())
                    if len(subdict) == 1:
                        for nkey, nvalue in subdict.iteritems():
                            if key in self.clstr_members[nkey]:
                                self.hmm_twobytwo[nkey][0] -= nvalue
                                self.hmm_twobytwo[nkey][0] += 1
                            else:
                                self.hmm_twobytwo[nkey][1] -= nvalue
                                self.hmm_twobytwo[nkey][1] += 1
        elif self.multicorrect and not self.truthset:
            for key, subdict in self.gene_multihits.iteritems():
                if subdict:
                    ## Correct for HMM counts
                    if len(subdict) > 1:
                        for nkey, nvalue in subdict.iteritems():
                            self.hmm_observed[nkey] -= nvalue
                            self.hmm_observed[nkey] += float(nvalue) / sum(subdict.values())
                    if len(subdict) == 1:
                        for nkey, nvalue in subdict.iteritems():
                            self.hmm_observed[nkey] -= nvalue
                            self.hmm_observed[nkey] += 1

    def aggregate_hierarchy(self):
        """
        Create the actual aggregated hierarchy counts based on the annotation file.  All counts are fundamentally
        aggregated from the HMM assignments, since this is comparing apples to apples with the truthset.  We are
        interested in the accuracy of HMM assignment and classification, not necessarily in the absolute annotation of
        the reads themselves.  By default, reads are split evenly if more than one annotation is present in the HMM
        when aggregating up.  The user is left to decide whether to round the output values or leave them as floats.
        :return: void
        """
        if self.truthset:
            for key, values in self.hmm_twobytwo.iteritems():
                annots = self.hmm_annots[key]
                class_annot = annots[0].split('|')
                mech_annot = annots[1].split('|')
                group_annot = annots[2].split('|')
                for entry in class_annot:
                    try:
                        list1 = np.array(self.class_twobytwo[entry])
                        list2 = np.array(values)
                        self.class_twobytwo[entry] = list1 + (list2 / float(len(class_annot)))
                    except KeyError:
                        self.class_twobytwo.setdefault(entry, np.array(values) / float(len(class_annot)))
                for entry in mech_annot:
                    try:
                        list1 = np.array(self.mech_twobytwo[entry])
                        list2 = np.array(values)
                        self.mech_twobytwo[entry] = list1 + (list2 / float(len(mech_annot)))
                    except KeyError:
                        self.mech_twobytwo.setdefault(entry, np.array(values) / float(len(mech_annot)))
                for entry in group_annot:
                    try:
                        list1 = np.array(self.group_twobytwo[entry])
                        list2 = np.array(values)
                        self.group_twobytwo[entry] = list1 + (list2 / float(len(group_annot)))
                    except KeyError:
                        self.group_twobytwo.setdefault(entry, np.array(values) / float(len(group_annot)))
        else:
            for key, value in self.hmm_observed.iteritems():
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

    def calculate_false(self):
        """
        This function is only utilized when a truthset is provided; it calculates the True/False negatives
        for each HMM and hierarchy two-by-two.  Note that this can only be done after all reads have been taken into
        account.
        :return: void
        """
        for key, values in self.hmm_twobytwo.iteritems():
            if key in self.hmm_truth:
                values[2] = int(self.hmm_truth[key]) - values[0]
                values[3] = sum(self.truthset_counts.itervalues()) - sum(values)
            else:
                values[2] = 0 - values[0]
        for key, values in self.class_twobytwo.iteritems():
            if key in self.class_truth:
                values[2] = int(self.class_truth[key]) - values[0]
                values[3] = sum(self.truthset_counts.itervalues()) - sum(values)
        for key, values in self.mech_twobytwo.iteritems():
            if key in self.mech_truth:
                values[2] = int(self.mech_truth[key]) - values[0]
                values[3] = sum(self.truthset_counts.itervalues()) - sum(values)
        for key, values in self.group_twobytwo.iteritems():
            if key in self.group_truth:
                values[2] = int(self.group_truth[key]) - values[0]
                values[3] = sum(self.truthset_counts.itervalues()) - sum(values)

    def write_twobytwos(self):
        if self.kmer:
            pathname = self.outpath + '/' + self.kmer + '_' + self.filename + '_%.0e.csv' % self.ethreshold
        else:
            pathname = self.outpath + '/' + self.filename + '.csv'
        with open(pathname, 'w') as twobytwo_file:
            twobytwo_file.write('Hierarchy,Name,True_Positive,False_Positive,False_Negative,True_Negative\n')
            for key, values in self.hmm_twobytwo.iteritems():
                twobytwo_file.write('HMM,' + key + ',' + ','.join([str(x) for x in values]) + '\n')
            for key, values in self.class_twobytwo.iteritems():
                twobytwo_file.write('Class,' + key + ',' + ','.join([str(x) for x in values]) + '\n')
            for key, values in self.mech_twobytwo.iteritems():
                twobytwo_file.write('Mechanism,' + key + ',' + ','.join([str(x) for x in values]) + '\n')
            for key, values in self.group_twobytwo.iteritems():
                twobytwo_file.write('Group,' + key + ',' + ','.join([str(x) for x in values]) + '\n')

    def write_observed(self):
        if self.kmer:
            pathname = self.outpath + '/' + self.kmer + '_' + self.filename + '_%.0e.csv' % self.ethreshold
        else:
            pathname = self.outpath + '/' + self.filename + '.csv'
        with open(pathname, 'w') as observed_file:
            observed_file.write('Hierarchy,Name,True_Positive,False_Positive,False_Negative,True_Negative\n')
            for key, values in self.hmm_observed.iteritems():
                observed_file.write('HMM,' + key + ',' + str(values) + '\n')
            for key, values in self.class_observed.iteritems():
                observed_file.write('Class,' + key + ',' + str(values) + '\n')
            for key, values in self.mech_observed.iteritems():
                observed_file.write('Mechanism,' + key + ',' + str(values) + '\n')
            for key, values in self.group_observed.iteritems():
                observed_file.write('Group,' + key + ',' + str(values) + '\n')

    def next(self):
        if not self.stdin and type(self.hmmer_file) is str:  # only open file here if hmmer_file is a str and not fileIO
            self.hmmer_file = open(self.hmmer_file, 'r')
        value = self._iterate
        if not value:  # close file on EOF
            if not self.stdin:
                self.hmmer_file.close()
            if self.multicorrect:
                self.correct_multihit()
            self.aggregate_hierarchy()
            if self.truthset:
                self.calculate_false()
                self.write_twobytwos()
            else:
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
parser.add_argument('graph_dir', type=str, help='Path to output directory for graphs')
parser.add_argument('clstr', type=str, help='Path to file containing clstr generation info')
parser.add_argument('annots', type=str, help='Path to annotation file')
parser.add_argument('-e', '--evalue', type=float, default=1e-14, help='Evalue under which to keep hits')
parser.add_argument('-t', '--truthset', nargs='?', default=None,
                    help='Path to file containing uniq -c style truth set counts')
parser.add_argument('-m', '--multicorrect', action='store_true', default=False,
                    help='If set, reads have a one to one mapping with reported hits')
parser.add_argument('-s', '--skewfile', nargs='?', default=None, help='Optional output file for HMM skewness metrics')
parser.add_argument('-k', '--kmer', nargs='?', default=None,
                    help='Optional k-mer length to append to filename for debug or testing purposes')

##########
## Main ##
##########
if __name__ == '__main__':
    ## Parse the arguments using ArgParse
    args = parser.parse_args()
    infile = args.input
    outpath = args.output
    graph_dir = args.graph_dir
    if not os.path.isdir(graph_dir):
        os.mkdir(graph_dir)
    if not os.path.isdir(outpath):
        os.mkdir(outpath)

    ## Determine alignment positions and overlay the reads onto the gene vectors.
    ## Calculate coverage and measures of skewness.
    ## Output a master file describing each mapped vector/metric and a coverage graph for each mapped gene.
    vector_hash = {}  # This stores gene names that were referenced in the SAM file, along with their vectors
    vector_counts = {}  # This stores gene names as in the other dictionary, but stores read counts instead
    for line in HmmerTime(infile, outpath, args.filename, args.hmm_len, args.evalue, args.multicorrect, args.truthset,
                          args.clstr, args.annots, args.kmer):
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
        if args.skewfile:
            with open(args.skewfile, 'w') as out:
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
                        norm_vec = value ** float(1) / sum(value)  # normalize so the vector sums to 1 (frequency)
                        max_entropy = np.ones(len(norm_vec)) / len(norm_vec)  # this is used to calculate the maximum shannon entropy for a given vector
                        shannon = np.negative(sum([x * np.log2(x) for x in norm_vec if x > 0])) / np.negative(
                                sum([x * np.log2(x) for x in max_entropy]))  # Shannon entropy
                        l2norm = 1 - ((np.sqrt(sum(norm_vec * norm_vec)) * np.sqrt(len(norm_vec))) - 1) / (
                            np.sqrt(len(norm_vec)) - 1)  # Deviation from the L2 norm unit sphere
                        out.write(",".join([key, str(len(value)), str(coverage), str(shannon), str(l2norm),
                                            " ".join([str(x) for x in value])]) + '\n')

                        ## Plot figure
                        plt.plot(value)
                        plt.xlabel('Nucleotide Position')
                        plt.ylabel('Observed Count')
                        plt.title('Coverage Plot for {}'.format(key))
                        plt.savefig(graph_dir + '/' + '{}.png'.format(key))
                        plt.close()  # make sure plot is closed
