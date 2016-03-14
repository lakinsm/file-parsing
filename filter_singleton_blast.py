#!/usr/bin/env python2

## This script filters the output results from BLASTing the HMM singletons against the nt BLAST database.
## Matches only pass filter if they meet the inclusion criteria of cd-hit (80% identity at 50% overlap of a given gene).

#############
## Imports ##
#############
import argparse
import numpy as np


##########
## Vars ##
##########
pass_filter = {}
identical = {}
gene_hash = set()
annodict = {}


#############
## Methods ##
#############
def index_max(values):
    return max(xrange(len(values)), key=values.__getitem__)


def filter_hits(infile, outfile, annotfile, tranfile):
    with open(infile, 'r') as hitfile, open(outfile, 'w') as out, open(tranfile, 'w') as tranout, open(annotfile, 'r') as annot:
        annotations = annot.read().strip().split('\n')
        [annodict.setdefault(x.split(',')[0], []).append(x.split(',')[1:4]) for x in annotations]
        data = hitfile.read().strip().split('\n')
        for hit in data:
            info = hit.split('\t')
            gene_id, target_acc, target_id, gene_len, gene_start, gene_stop, target_len, _, _, evalue, pident, \
                _, _, target_title, _, _, _, seq = info
            name = str(target_id)+str(target_acc)
            gene_hash.add(gene_id)
            if (int(gene_stop) - int(gene_start)) < (int(gene_len) / 2):
                continue
            if int(float(pident)) is 100:
                identical.setdefault(gene_id, []).append((name, int(gene_stop) - int(gene_start), seq))
                continue
            pass_filter.setdefault(gene_id, []).append((name, seq))
        for gene in gene_hash:
            if gene not in pass_filter.keys():
                if gene not in identical.keys():
                    print(gene)
                    continue
                zipped = zip(*identical[gene])
                max_indices = np.where(np.array(zipped[1]) == max(zipped[1]))[0].tolist()
                [pass_filter.setdefault(gene, []).append((y[0], y[2])) for y in [identical[gene][x] for x in max_indices]]
            if gene not in pass_filter.keys():
                raise Exception('Should not reach this line')
        for key, entry in pass_filter.iteritems():
            for hit in entry:
                tranout.write(hit[0]+'|singleton_addition,'+','.join(annodict[key][0])+'\n')
                out.write('>'+hit[0]+'|singleton_addition\n'+hit[1].replace('-','')+'\n')


##############
## ArgParse ##
##############
parser = argparse.ArgumentParser('amr_skewness.py')
parser.add_argument('infile', type=str, help='File path to custom BLAST tabular output file')
parser.add_argument('annotation_file', type=str, help='File path to annotation file for translation')
parser.add_argument('outfile', type=str, help='File path to output file')
parser.add_argument('new_annotations_file', type=str, help='File path to output new annotations')
args = parser.parse_args()


##########
## Main ##
##########
if __name__ == '__main__':
    filter_hits(args.infile, args.outfile, args.annotation_file, args.new_annotations_file)

