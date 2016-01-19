#!/usr/bin/env python

"""Parse a HMMer tblout file line by line, storing relevant information.  Output a vector for each gene that allows
for subsequent calculation of coverage and skewness of the read distribution for each HMM.  If a truthset is provided,
calculate two-by-two accuracies. Optionally output a ROC graph and statistics flat file.
"""

# Header values for reference:
#1 target name
#2 accession
#3 query name
#4 accession
#5 hmmfrom
#6 hmm to
#7 alifrom
#8 ali to
#9 envfrom
#10 env to
#11 sq len
#12 strand
#13 E-value
#14 score
#15 bias
#16 description of target (multi-value)


## To do list:
# splitting of multi-hits (need to check region?)
# Class, mech, gene level two-by-twos

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
class HmmerWalk:
    """This object takes as input a HMMer tblout file path and constructs an iterable that outputs
    hash-mapping of header to sequence information.  Only one line will be held in memory at a time using this method.
    The object walks along the file and, if a truthset is provided, outputs two-by-two values for accuracy.
    """
    def __init__(self, filepath, length, evalue=10, diff=0, truthset=None, clstr_file=None, annot_file=None,):
        """
        constructor
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
        self.reads_mapping = 0
        self.reads_total = 0
        self.ethreshold = float(evalue)
        self.observed_counts = {}
        self.observed_group = {}
        self.observed_mech = {}
        self.observed_class = {}
        self.clstr_members = {}
        self.truthset = False
        self.hmm_lengths = {}
        with open(length, 'r') as hmm_length:
            data = hmm_length.read().split('\n')[1:]
            for line in data:
                line = line.split()
                if line:
                    self.hmm_lengths.setdefault(line[0], int(line[1]))
        if truthset:
            if not clstr_file or not annot_file:
                raise ValueError("If truthset is provided, all non-default parameters must be defined")
            self.truthset = True
            self.truthset_counts = {}
            self.off_target_hits = {}
            with open(truthset, 'r') as truth:
                data = truth.read().split('\n')
                for line in data:
                    temp = line.split()
                    if temp:
                        self.truthset_counts.setdefault(temp[1], int(temp[0]))
                        self.observed_counts.setdefault(temp[1], 0)
            self.hmm_twobytwo = {}
            with open(clstr_file, 'r') as f:
                header_reg = re.compile(r'>(.+?)\.\.\.')
                line = f.readline()
                cluster_num = -2
                clstr = []
                while line:
                    if line[0] is ">":
                        cluster_num += 1
                        if cluster_num > 0:
                            if len(clstr) is 1:
                                clstr = []
                            else:
                                self.clstr_members.setdefault(str(cluster_num), clstr)
                                self.hmm_twobytwo.setdefault(str(cluster_num), [0, 0, 0, 0])  # TP, FP, FN, TN
                                clstr = []
                    else:
                        clstr.append(header_reg.findall(line)[0])
                    line = f.readline()
            self.group_twobytwo = {}
            self.mech_twobytwo = {}
            self.class_twobytwo = {}
            self.gene_annots = {}
            with open(annot_file, 'r') as annot:
                data = annot.read().split('\n')
                for line in data:
                    temp = line.split(',')
                    if temp[0]:
                        self.gene_annots.setdefault(temp[0], temp[1:])
                        self.class_twobytwo.setdefault(temp[1], [0, 0, 0, 0])
                        if temp[2]:
                            self.mech_twobytwo.setdefault(temp[2], [0, 0, 0, 0])
                        if temp[3]:
                            self.group_twobytwo.setdefault(temp[3], [0, 0, 0, 0])
            self.hmm_annots = {}
            for key, values in self.clstr_members.iteritems():
                classes, mechs, groups = zip(*[self.gene_annots[x] for x in values])
                classes = [x for x in classes if x]
                mechs = [x for x in mechs if x]
                groups = [x for x in groups if x]
                self.hmm_annots.setdefault(key, ['|'.join(set(classes)), '|'.join(set(mechs)), '|'.join(set(groups))])
            self.class_counts = {}
            self.mech_counts = {}
            self.group_counts = {}
            for key, value in self.truthset_counts.iteritems():
                if key in self.gene_annots:
                    try:
                        gene_class = self.gene_annots[key][0]
                        self.class_counts[gene_class] += value
                    except KeyError:
                        self.class_counts.setdefault(gene_class, value)
                    try:
                        gene_mech = self.gene_annots[key][1]
                        self.mech_counts[gene_mech] += value
                    except KeyError:
                        self.mech_counts.setdefault(gene_mech, value)
                    try:
                        gene_group = self.gene_annots[key][2]
                        self.group_counts[gene_group] += value
                    except KeyError:
                        self.group_counts.setdefault(gene_group, value)

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
                if float(temp[12]) < float(self.ethreshold):
                    self.reads_mapping += 1
                    ## Basic observation increment rules
                    if self.truthset:
                        try:
                            self.observed_counts[temp[0]] += 1
                        except KeyError:
                            try:
                                self.off_target_hits[temp[9]] += 1
                            except KeyError:
                                self.off_target_hits.setdefault(temp[0], 1)
                        ## TRUTHSET ENABLED: What category of HMM-level two-by-two does this hit fall under?
                        ## Can only calculate TP/FP here; the others are done at the end
                        if temp[0] in self.clstr_members[temp[2]]:
                            self.hmm_twobytwo[temp[2]][0] += 1
                        else:
                            self.hmm_twobytwo[temp[2]][1] += 1
                        ## TRUTHSET ENABLED: What category of Group-level two-by-two does this hit fall under?
                        ## Can only calculate TP/FP here; the others are done at the end
                        hmm_annot = self.hmm_annots[temp[2]]
                        hmm_groups = hmm_annot[2].split('|')
                        if temp[0] in self.gene_annots:
                            gene_group = self.gene_annots[temp[0]][2]
                            if gene_group and gene_group in hmm_groups:
                                self.group_twobytwo[gene_group][0] += 1
                            elif gene_group and gene_group not in hmm_groups:
                                self.group_twobytwo[gene_group][1] += 1
                        ## TRUTHSET ENABLED: What category of Mechanism-level two-by-two does this hit fall under?
                        ## Can only calculate TP/FP here; the others are done at the end
                        hmm_mechs = hmm_annot[1].split('|')
                        if temp[0] in self.gene_annots:
                            gene_mech = self.gene_annots[temp[0]][1]
                            if gene_mech and gene_mech in hmm_mechs:
                                self.mech_twobytwo[gene_mech][0] += 1
                            elif gene_mech and gene_mech not in hmm_mechs:
                                self.mech_twobytwo[gene_mech][1] += 1
                        ## TRUTHSET ENABLED: What category of Class-level two-by-two does this hit fall under?
                        ## Can only calculate TP/FP here; the others are done at the end
                        hmm_classes = hmm_annot[0].split('|')
                        if temp[0] in self.gene_annots:
                            gene_class = self.gene_annots[temp[0]][0]
                            if gene_class and gene_class in hmm_classes:
                                self.class_twobytwo[gene_class][0] += 1
                            elif gene_class and gene_class not in hmm_classes:
                                self.class_twobytwo[gene_class][1] += 1
                        ## TRUTHSET ENABLED: Increment the hits on each gene
                        hmm_classes = hmm_annot[0].split('|')
                        if temp[0] in self.gene_annots:
                            gene_class = self.gene_annots[temp[0]][0]
                            if gene_class and gene_class in hmm_classes:
                                self.class_twobytwo[gene_class][0] += 1
                            elif gene_class and gene_class not in hmm_classes:
                                self.class_twobytwo[gene_class][1] += 1
                        ## TRUTHSET ENABLED: Add each hit on HMM to gene hash
                        hmm_classes = hmm_annot[0].split('|')
                        if temp[0] in self.gene_annots:
                            gene_class = self.gene_annots[temp[0]][0]
                            if gene_class and gene_class in hmm_classes:
                                self.class_twobytwo[gene_class][0] += 1
                            elif gene_class and gene_class not in hmm_classes:
                                self.class_twobytwo[gene_class][1] += 1
                    else:
                        try:
                            self.observed_counts[temp[0]] += 1
                        except KeyError:
                            self.observed_counts.setdefault(temp[0], 1)
                    return temp[2], self.hmm_lengths[temp[2]], temp[4], temp[5], temp[11]  # name, len, start, stop, str
        self.hmmer_file.close()  # catch all in case this line is reached
        assert False, 'Should not reach this line'

    #def calculate_false(self):
        ## Bookmark

    def next(self):
        if not self.stdin and type(self.hmmer_file) is str:  # only open file here if hmmer_file is a str and not fileIO
            self.hmmer_file = open(self.hmmer_file, 'r')
        value = self._iterate
        if not value:  # close file on EOF
            if not self.stdin:
                self.hmmer_file.close()
            sum = 0
            for key, value in self.truthset_counts.iteritems():
                if key in self.gene_annots:
                    print key, value
                    sum += value
            print sum

            #print([x for x in self.observed_counts.itervalues() if x])
            #self.calculate_false()
            #self.write_stats()  # Write the calculated dictionaries to the appropriate files (WIP)
            ## Remember to calculate the true/false negatives here
            ## Also to calculate observed aggregated values
            raise StopIteration()
        else:
            return value



##############
## ArgParse ##
##############
parser = argparse.ArgumentParser('amr_skewness.py')
parser.add_argument('input', type=str, help='File path to HMMer tblout file, or "-" for stdin')
parser.add_argument('hmm_len', type=str, help='Path to file containing HMM lengths')
parser.add_argument('outputfile', type=str, help='File path to desired output file (.csv format)')
parser.add_argument('graph_dir', type=str, help='Path to output directory for graphs')
parser.add_argument('--evalue', type=float, default=10, help='Evalue under which to keep hits')
parser.add_argument('--truthset', nargs='?', default=None, help='Path to file containing uniq -c style truth set counts')
parser.add_argument('--clstr', nargs='?', default=None, help='Path to file containing clstr generation info')
parser.add_argument('--annots', nargs='?', default=None, help='Path to annotation file')
parser.add_argument('--diff', nargs ='?', default=0, help='Difference needed to declare an Evalue truth over another')


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
        vector_counts = {}  # This stores gene names as in the other dictionary, but stores read counts instead
        for line in HmmerWalk(infile, args.hmm_len, args.evalue, args.diff, args.truthset, args.clstr, args.annots):
            if int(line[2]) < int(line[3]):
                start = line[2]
                stop = line[3]
            else:
                start = line[3]
                stop = line[2]
            try:
                vector_hash[line[0]][(int(start)-1):(int(stop))] += 1  # increment affected region
            except KeyError:
                vector_hash[line[0]] = np.zeros(int(line[1])).astype('int')  # new entry, initialize vector

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
                max_entropy = np.ones(len(norm_vec)) / len(norm_vec)  # this is used to calculate the maximum shannon entropy for a given vector
                shannon = np.negative(sum([x*np.log2(x) for x in norm_vec if x > 0])) / np.negative(sum([x*np.log2(x) for x in max_entropy])) # Shannon entropy
                l2norm = 1 - ((np.sqrt(sum(norm_vec*norm_vec))*np.sqrt(len(norm_vec))) - 1) / (np.sqrt(len(norm_vec)) - 1)  # Deviation from the L2 norm unit sphere
                out.write(",".join([key, str(len(value)), str(coverage), str(shannon), str(l2norm), " ".join([str(x) for x in value])])+'\n')

                ## Plot figure
                plt.plot(value)
                plt.xlabel('Nucleotide Position')
                plt.ylabel('Observed Count')
                plt.title('Coverage Plot for {}'.format(key))
                plt.savefig(graph_dir+'/'+'{}.png'.format(key))
                plt.close()  # make sure plot is closed
