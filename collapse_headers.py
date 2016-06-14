#!/usr/bin/env python3

import sys

annot_table = {}


def load_annots(infile):
    global annot_table
    with open(infile, 'r') as afile:
        data = afile.read().strip().split('\n')
        for entry in data[1:]:
            entry = entry.split(',')
            annot_table.setdefault('|'.join(entry[0].split('|')[:-3]), entry[1:])


def fasta_parse(infile):
    """ Parses a fasta file in chunks of 2 lines.
    :param infile: path to the input fasta file
    :return: generator of (header, sequence) fasta tuples
    """
    with open(infile, 'r') as fasta_file:
        # Skip whitespace
        while True:
            line = fasta_file.readline()
            if line is "":
                return  # Empty file or premature end of file?
            if line[0] is ">":
                break
        while True:
            if line[0] is not ">":
                raise ValueError("Records in FASTA should begin with '>'")
            header = line[1:].rstrip()
            all_lines = []
            line = fasta_file.readline()
            while True:
                if not line:
                    break
                if line[0] is ">":
                    break
                all_lines.append(line.rstrip())
                line = fasta_file.readline()
            yield header, "".join(all_lines).replace(" ", "").replace("\r", "")
            if not line:
                return  # Stop Iteration
        assert False, "Should not reach this line"


if __name__ == '__main__':
    load_annots(sys.argv[1])
    with open('production_resistance_annotations.csv', 'w') as out:
        for header, seq in fasta_parse(sys.argv[2]):
            header = '|'.join(header.split('|')[:-3])
            if header not in annot_table:
                sys.stderr.write('{}\n'.format(header))
            else:
                header2 = header + '|'+'|'.join(annot_table[header]).replace(' ', '_')
                sys.stdout.write('>{}\n{}\n'.format(header2, seq))
                annots = annot_table[header]
                out.write('{},{},{},{}\n'.format(header2, annots[0], annots[1], annots[2]))

