#!/usr/bin/env python3

## Locate the HMM model numbers for the MEG HMMs where the HMM was constructed from an augmented singleton cluster

#############
## Imports ##
#############
import re
import argparse


##########
## Vars ##
##########


#############
## Methods ##
#############
def find_singletons(clstr_file):
    clstr_members = {}
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
                        clstr_members.setdefault(str(cluster_num), clstr)
                        clstr = []
            else:
                clstr.append(header_reg.findall(line)[0])  # Begin appending headers before pushing to dict
            line = f.readline()
        singleton_reg = re.compile(r'singleton_addition')
        singleton_members = {}
        for key, value in clstr_members.items():
            if any([singleton_reg.findall(x) for x in value]):
                singleton_members.setdefault(key, value)
    return singleton_members


##############
## ArgParse ##
##############
parser = argparse.ArgumentParser('locate_singleton_hmms.py')
parser.add_argument('input', type=str, help='File path to HMM cluster file (output from CD-HIT')
parser.add_argument('output', type=str, help='File path to singleton cluster output (.txt)')


##########
## Main ##
##########
if __name__ == '__main__':
    ## Parse the arguments using ArgParse
    args = parser.parse_args()
    infile = args.input
    data = find_singletons(infile)
    with open(args.output, 'w') as out:
        for key, value in data.items():
            out.write(key+'\n')
