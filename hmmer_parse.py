#!/usr/bin/env python

"""Parse a HMMer tblout file line by line, storing relevant information.  Output a vector for each gene that allows
for subsequent calculation of coverage and skewness of the read distribution for each HMM.  If a truthset is provided,
calculate two-by-two accuracies. Optionally output a ROC graph and statistics flat file.
Note: Reads in the original fasta/fastq file MUST have a unique header for this to work properly.
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
    def __init__(self, filepath, length, evalue=10, diff=0, truthset=None, clstr_file=None, annot_file=None,):
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
        self.reads_mapping = 0
        self.reads_total = 0
        self.ethreshold = float(evalue)
        self.observed_counts = {}  # Observed counts for each gene
        self.observed_group = {}  # Observed counts for each group, aggregated values
        self.observed_mech = {}  # Observed counts for each mechanism, aggregated values
        self.observed_class = {}  # Observed counts for each class, aggregated values
        self.clstr_members = {}  # Mapping of hmm # -> genes used to build that HMM
        self.truthset = False
        self.hmm_lengths = {}  # Length of each hmm, hmm # -> length
        self.gene_multihits = {}  # Will store HMM hits for each gene
        self.gene_multihit_evalues = {}  # Will store evalues for each hit
        with open(length, 'r') as hmm_length:
            data = hmm_length.read().split('\n')[1:]
            for line in data:
                line = line.split()
                if line:
                    self.hmm_lengths.setdefault(line[0], int(line[1]))
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
            if not clstr_file or not annot_file:
                raise ValueError("If truthset is provided, all non-default parameters must be defined")
            self.truthset = True
            self.truthset_counts = {}  # True counts for each gene
            self.off_target_hits = {}  # If a hit doesn't map to its proper (truthset) or any (normal) HMM, put it here
            with open(truthset, 'r') as truth:
                ## Put each gene header into a key and its counts into the value.  Initialize the obs count dict
                data = truth.read().split('\n')
                for line in data:
                    temp = line.split()
                    if temp:
                        self.truthset_counts.setdefault(temp[1], int(temp[0]))
                        self.observed_counts.setdefault(temp[1], 0)
            self.hmm_twobytwo = {}  # HMM matrix for ROC generation
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
            self.group_twobytwo = {}  # Group matrix for ROC generation
            self.mech_twobytwo = {}  # Mechanism matrix for ROC generation
            self.class_twobytwo = {}  # Class matrix for ROC generation
            self.gene_annots = {}  # Mapping from gene name -> annotations
            with open(annot_file, 'r') as annot:
                ## Map the gene names to their annotations, initialize the hierarchy two-by-twos
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
            self.hmm_annots = {}  # HMM annotation mapping from key (hmm #) -> annotations
            for key, values in self.clstr_members.iteritems():
                ## Calculate the annotations of each HMM, combining around a pipe '|' if multiple
                classes, mechs, groups = zip(*[self.gene_annots[x] for x in values])
                classes = [x for x in classes if x]
                mechs = [x for x in mechs if x]
                groups = [x for x in groups if x]
                self.hmm_annots.setdefault(key, ['|'.join(set(classes)), '|'.join(set(mechs)), '|'.join(set(groups))])
            self.hmm_truth = {}  # True counts for each HMM
            self.class_truth = {}  # True counts for each class
            self.mech_truth = {}  # True counts for each mechanism
            self.group_truth = {}  # True counts for each group
            for key, values in self.clstr_members.iteritems():
                ## Generate true counts for each HMM
                total = sum([int(self.truthset_counts[x]) for x in values if x in self.truthset_counts])
                self.hmm_truth.setdefault(key, total)
            for key, value in self.truthset_counts.iteritems():
                if key in self.gene_annots:
                    ## Generate true counts for the hierarchy by passing gene name through the annotations dict
                    if self.gene_annots[key][0]:
                        class_entry = self.gene_annots[key][0]
                        try:
                            self.class_truth[class_entry] += value
                        except KeyError:
                            self.class_truth.setdefault(class_entry, value)
                    if self.gene_annots[key][1]:
                        mech_entry = self.gene_annots[key][1]
                        try:
                            self.mech_truth[mech_entry] += value
                        except KeyError:
                            self.mech_truth.setdefault(mech_entry, value)
                    if self.gene_annots[key][2]:
                        group_entry = self.gene_annots[key][2]
                        try:
                            self.group_truth[group_entry] += value
                        except KeyError:
                            self.group_truth.setdefault(group_entry, value)

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
                            if temp[2] == '989':
                                print temp[0]
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
                        ## TRUTHSET ENABLED: Keep track of whether a gene hits across multiple HMMs
                        try:
                            self.gene_multihits[temp[0]][temp[2]] += 1
                        except (TypeError, KeyError):
                            self.gene_multihits.setdefault(temp[0], {temp[2]: 1})
                    else:
                        try:
                            self.observed_counts[temp[0]] += 1
                        except KeyError:
                            self.observed_counts.setdefault(temp[0], 1)
                        try:
                            self.gene_multihits[temp[0]][temp[2]] += 1
                        except (TypeError, KeyError):
                            self.gene_multihits.setdefault(temp[0], {temp[2]: 1})
                        try:
                            self.gene_multihit_evalues[temp[0]][temp[2]] += (temp[12], )
                        except (TypeError, KeyError):
                            self.gene_multihit_evalues.setdefault(temp[0], {temp[2]: (temp[12], )})
                    return temp[2], self.hmm_lengths[temp[2]], temp[4], temp[5], temp[11]  # name, len, start, stop, str
        self.hmmer_file.close()  # catch all in case this line is reached
        assert False, 'Should not reach this line'

    def correct_multihit(self):
        """
        For our purposes, reads need to maintain a one-to-one mapping with hits counted.  Therefore, we both need to
        correct for hits across multiple HMMs due to homology and also for different motifs within the same read
        hitting multiple times within the same HMM.  This correction will be optional in the production version of
        the script.  Note that for this to work, reads must have a unique header in the original fastq/fasta file.
        :return: void
        """
        for key, subdict in self.gene_multihits.iteritems():
            if subdict:
                multiclass = {self.hmm_annots[k][0]:v for k, v in subdict.iteritems()}
                multimech = {self.hmm_annots[k][1]:v for k, v in subdict.iteritems()}
                multigroup = {self.hmm_annots[k][2]:v for k, v in subdict.iteritems()}
                ## Correct for HMM counts
                if len(subdict) > 1:
                    for nkey, nvalue in subdict.iteritems():
                        if key in self.clstr_members[nkey]:
                            self.hmm_twobytwo[nkey][0] -= nvalue
                            self.hmm_twobytwo[nkey][0] += float(nvalue) / sum(subdict.values())
                        else:
                            self.hmm_twobytwo[nkey][1] -= nvalue
                            self.hmm_twobytwo[nkey][1] += float(nvalue) / sum(subdict.values())
                ## Correct for Class counts
                if len(multiclass) > 1:
                    for nkey, nvalue in multiclass.iteritems():
                        if self.gene_annots[key][0]:
                            if nkey is self.gene_annots[key][0]:
                                self.class_twobytwo[nkey][0] -= nvalue
                                self.class_twobytwo[nkey][0] += float(nvalue) / sum(multiclass.values())
                            else:
                                self.class_twobytwo[nkey][1] -= nvalue
                                self.class_twobytwo[nkey][1] += float(nvalue) / sum(multiclass.values())
                ## Correct for Mech counts
                if len(multimech) > 1:
                    for nkey, nvalue in multimech.iteritems():
                        if self.gene_annots[key][1]:
                            if nkey is self.gene_annots[key][1]:
                                self.mech_twobytwo[nkey][0] -= nvalue
                                self.mech_twobytwo[nkey][0] += float(nvalue) / sum(multimech.values())
                            else:
                                self.mech_twobytwo[nkey][1] -= nvalue
                                self.mech_twobytwo[nkey][1] += float(nvalue) / sum(multimech.values())
                ## Correct for Group counts
                if len(multigroup) > 1:
                    for nkey, nvalue in multigroup.iteritems():
                        if self.gene_annots[key][2]:
                            if nkey is self.gene_annots[key][2]:
                                self.group_twobytwo[nkey][0] -= nvalue
                                self.group_twobytwo[nkey][0] += float(nvalue) / sum(multigroup.values())
                            else:
                                self.group_twobytwo[nkey][1] -= nvalue
                                self.group_twobytwo[nkey][1] += float(nvalue) / sum(multigroup.values())
                ## Now correct for within-HMM multihits by reducing any multiple read hits to a single hit
                elif len(subdict) == 1:
                    for nkey, nvalue in subdict.iteritems():
                        if key in self.clstr_members[nkey]:
                            self.hmm_twobytwo[nkey][0] -= nvalue
                            self.hmm_twobytwo[nkey][0] += 1
                        else:
                            self.hmm_twobytwo[nkey][1] -= nvalue
                            self.hmm_twobytwo[nkey][1] += 1
                ## Correct for Class counts
                if len(multiclass) == 1:
                    for nkey, nvalue in multiclass.iteritems():
                        if self.gene_annots[key][0]:
                            if nkey is self.gene_annots[key][0]:
                                self.class_twobytwo[nkey][0] -= nvalue
                                self.class_twobytwo[nkey][0] += 1
                            else:
                                self.class_twobytwo[nkey][1] -= nvalue
                                self.class_twobytwo[nkey][1] += 1
                ## Correct for Mech counts
                if len(multimech) == 1:
                    for nkey, nvalue in multimech.iteritems():
                        if self.gene_annots[key][1]:
                            if nkey is self.gene_annots[key][1]:
                                self.mech_twobytwo[nkey][0] -= nvalue
                                self.mech_twobytwo[nkey][0] += 1
                            else:
                                self.mech_twobytwo[nkey][1] -= nvalue
                                self.mech_twobytwo[nkey][1] += 1
                ## Correct for Group counts
                if len(multigroup) == 1:
                    for nkey, nvalue in multigroup.iteritems():
                        if self.gene_annots[key][2]:
                            if nkey is self.gene_annots[key][2]:
                                self.group_twobytwo[nkey][0] -= nvalue
                                self.group_twobytwo[nkey][0] += 1
                            else:
                                self.group_twobytwo[nkey][1] -= nvalue
                                self.group_twobytwo[nkey][1] += 1

    def calculate_false(self):
        for key, values in self.hmm_twobytwo.iteritems():
            if key in self.hmm_truth:
                values[2] = int(self.hmm_truth[key]) - values[0]
                values[3] = sum(self.observed_counts.itervalues()) - sum(values)
        for key, values in self.class_twobytwo.iteritems():
            if key in self.class_truth:
                values[2] = int(self.class_truth[key]) - values[0]
                values[3] = sum(self.observed_counts.itervalues()) - sum(values)
        for key, values in self.mech_twobytwo.iteritems():
            if key in self.mech_truth:
                values[2] = int(self.mech_truth[key]) - values[0]
                values[3] = sum(self.observed_counts.itervalues()) - sum(values)
        for key, values in self.group_twobytwo.iteritems():
            if key in self.group_truth:
                values[2] = int(self.group_truth[key]) - values[0]
                values[3] = sum(self.observed_counts.itervalues()) - sum(values)

    def next(self):
        if not self.stdin and type(self.hmmer_file) is str:  # only open file here if hmmer_file is a str and not fileIO
            self.hmmer_file = open(self.hmmer_file, 'r')
        value = self._iterate
        if not value:  # close file on EOF
            if not self.stdin:
                self.hmmer_file.close()
            self.correct_multihit()
            self.calculate_false()
            for key, value in self.class_twobytwo:
                print key, value
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
        for line in HmmerTime(infile, args.hmm_len, args.evalue, args.diff, args.truthset, args.clstr, args.annots):
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

        # out.write('Accession_Name,Accession_Length,Coverage,Shannon_Entropy,L2norm_Deviation,Vector\n')  # headers for outfile
        # plt.ioff()  # no interactive mode
        # plt.hold(False)  # don't keep plot
        # for key, value in vector_hash.iteritems():
        #     if sum(value > 0):  # if the gene had reads aligned
        #         plot_count += 1
        #         sys.stdout.write("\rPlots generated: {}".format(plot_count))  # update counter for user benefit
        #         sys.stdout.flush()
        #
        #         ## Calculate metrics
        #         coverage = float(sum(value > 0)) / len(value)  # what percentage of the gene has a read aligned?
        #         norm_vec = value**float(1) / sum(value)  # normalize so the vector sums to 1 (frequency)
        #         max_entropy = np.ones(len(norm_vec)) / len(norm_vec)  # this is used to calculate the maximum shannon entropy for a given vector
        #         shannon = np.negative(sum([x*np.log2(x) for x in norm_vec if x > 0])) / np.negative(sum([x*np.log2(x) for x in max_entropy])) # Shannon entropy
        #         l2norm = 1 - ((np.sqrt(sum(norm_vec*norm_vec))*np.sqrt(len(norm_vec))) - 1) / (np.sqrt(len(norm_vec)) - 1)  # Deviation from the L2 norm unit sphere
        #         out.write(",".join([key, str(len(value)), str(coverage), str(shannon), str(l2norm), " ".join([str(x) for x in value])])+'\n')
        #
        #         ## Plot figure
        #         plt.plot(value)
        #         plt.xlabel('Nucleotide Position')
        #         plt.ylabel('Observed Count')
        #         plt.title('Coverage Plot for {}'.format(key))
        #         plt.savefig(graph_dir+'/'+'{}.png'.format(key))
        #         plt.close()  # make sure plot is closed
