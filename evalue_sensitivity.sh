#!/bin/bash

## This script analyses the sensitivity of hmmer output to the choice
## of evalue in hmmer_parse.py

## Author: Steven Lakin
## Email: Steven.Lakin@colostate.edu

##########
## Vars ##
##########
skip=""
threads=1
j=1
k=-25
mu=0
num=1

shopt -s extglob


#############
## Methods ##
#############
display_help() {
    echo "
        Usage: kmer_sensitivity.sh [options]

            -h | --help		Display this menu

        Input Options

            -s | --skip don't perform the analysis (debugging)
            -t | --threads INT	Number of threads to use [1]

    "
}


##########
## Main ##
##########
while [[ "${1+defined}" ]]; do
    case "$1" in
        -i | --input )
            infiles=("$2")
            shift 2
            ;;
        -h | --help )
            display_help
            exit 0
            ;;
        -t | --threads )
            threads=$2
            shift 2
            ;;
        --)  # End of options
            shift 1
            break
            ;;
        -*)
            echo "Error: Unknown option $1" >&2
            exit 1
            ;;
        *)  # No more options
            break
            ;;
    esac
done


## Create dirs if not present
if [[ ! -d outputs ]]; then
    mkdir outputs
fi


## Generate tblout files for each kmer choice
for i in ${infiles[@]}; do
    j=1
    k=-25
    mu=0
    num=1
	while [[ $(echo "if (${mu} < 10) 1 else 0" | bc) -eq 1 ]]; do
		mu=$(echo "scale=25;$j*10^$k" | bc)
		handle=$( basename "$i" | sed 's/.tblout.scan//' )
		kmer=$( echo "$handle" | grep -o -P "^\d{2}" )
		if [ -e $kmer ]; then
		    kmer=0
		fi
		#echo "$num" "$mu" >> evalues.txt
		hmmer_parse.py -m -e "$mu" -t small_testset_counts.txt -k "$kmer" "$i" outputs/full_analysis truthset_screened hmm_lengths.tsv outputs/graphs no_singletons.fasta.clstr final_hmmer_annotations.csv
		[[ j -eq 1 ]] && j=4
		[[ j -eq 5 ]] && j=0 && ((k++))
		((j++))
		((num++))
	done
done


exit 0

