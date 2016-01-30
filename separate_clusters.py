#!/usr/bin/env python

## This script generates fasta files that will be fed into muscle for MPA

import re

singletons = '/s/bovine/index/databases/resistance_databases/analyze/lakinsm/hmm/remaining_singletons.fasta'
infile = '/s/bovine/index/databases/resistance_databases/analyze/lakinsm/hmm/clstr/no_singletons.fasta.clstr'
outdir = '/s/bovine/index/databases/resistance_databases/analyze/lakinsm/hmm/fasta/'
database = '/s/bovine/index/databases/resistance_databases/analyze/lakinsm/hmm/final_hmmer_templates.fasta'
lite = '/s/bovine/index/databases/resistance_databases/analyze/lakinsm/hmm/low_freq_hmm_names.txt'
litedir = '/s/bovine/index/databases/resistance_databases/analyze/lakinsm/hmm/lite_fasta/'


def fastaParse(infile):
    with open(infile, 'r') as fastaFile:
        # Skip whitespace
        while True:
            line = fastaFile.readline()
            if line is "":
                return  # Empty file or premature end of file?
            if line[0] is ">":
                break
        while True:
            if line[0] is not ">":
                raise ValueError("Records in FASTA should begin with '>'")
            header = line[1:].rstrip()
            allLines = []
            line = fastaFile.readline()
            while True:
                if not line:
                    break
                if line[0] is ">":
                    break
                allLines.append(line.rstrip())
                line = fastaFile.readline()
            yield header, "".join(allLines).replace(" ", "").replace("\r", "")
            if not line:
                return  # Stop Iteration
        assert False, "Should not reach this line"


def separate_clusters():
    cluster_num = 0
    clstr = []
    header_reg = re.compile(r'>(.+?)\.\.\.')
    with open(infile, 'r') as f, open(singletons, 'w') as s, open(lite, 'r') as hmmlite:
        dbdata = {header: seq for header, seq in fastaParse(database)}
        lite_hmms = [int(x) for x in hmmlite.read().split('\n') if x]
        print(len(lite_hmms))
        line = f.readline()
        while line:
            if line[0] is ">":
                cluster_num += 1
                if cluster_num > 1:
                    if len(clstr) is 1:
                        for header in clstr:
                            s.write('>' + header + '\n' + dbdata[header] + '\n')
                        clstr = []
                    else:
                        with open(outdir + str(cluster_num - 1) + '.fasta', 'w') as out:
                            for header in clstr:
                                out.write('>' + header + '\n' + dbdata[header] + '\n')
                            if cluster_num - 1 not in lite_hmms:
                                with open(litedir + str(cluster_num - 1) + '.fasta', 'w') as outlite:
                                    for header in clstr:
                                        outlite.write('>' + header + '\n' + dbdata[header] + '\n')
                            clstr = []
            else:
                clstr.append(header_reg.findall(line)[0])
            line = f.readline()
        print(cluster_num)


if __name__ == '__main__':
    separate_clusters()
