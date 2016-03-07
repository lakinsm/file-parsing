#!/usr/bin/env python

## Download all gene sequences from a few model organisms for use in the simulated truth set for the HMM testing


#############
## Imports ##
#############
from Bio import Entrez
import itertools
import time
from urllib2 import HTTPError

##########
## Vars ##
##########
iteration = 0
#'txid9913', 'txid7227', 'txid10090',


##########
## Main ##
##########
Entrez.email = 'Steven.Lakin@colostate.edu'
search = ['txid11103', 'txid4932']
with open('/home/lakinsm/Documents/phd/phdenv/benchmarking_hmm/simulated_truthset/outfile.fasta', 'aw') as outfile:
    for term in search:
        print term
        handle = Entrez.read(Entrez.esearch(db="sequences", term=term, retmax=2000))
        if len(handle['IdList']) > 200:
            gis = [x for x in itertools.izip_longest(*(iter(handle['IdList']),) * 200)]
        else:
            gis = [handle['IdList']]
        for group in gis:
            clean = [x for x in group if x]
            iteration += 1
            print "Group {}".format(iteration)
            request = Entrez.epost('gene', id=",".join(map(str, clean)))
            result = Entrez.read(request)
            webEnv = result['WebEnv']
            queryKey = result['QueryKey']
            attempt = 1
            try:
                handle = Entrez.efetch(db='nucleotide', rettype='fasta', retmode='text', webenv=webEnv, query_key=queryKey)
            except HTTPError as err:
                if 500 <= err.code <= 599:
                    print 'Attempt {} of 3'.format(attempt)
                    attempt +=1
                    time.sleep(15)
                    raise
                else:
                    raise
            data = handle.read()
            handle.close()
            outfile.write(data)
