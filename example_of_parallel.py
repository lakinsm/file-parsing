#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# multiproc_sums.py
"""A program that reads integer values from a CSV file and writes out their
sums to another CSV file, using multiple processes if desired.
"""

import csv
import multiprocessing
import optparse
import sys

NUM_PROCS = multiprocessing.cpu_count()

def make_cli_parser():
    """Make the command line interface parser."""
    usage = "\n\n".join(["python %prog INPUT_CSV OUTPUT_CSV",
            __doc__,
            """
ARGUMENTS:
    INPUT_CSV: an input CSV file with rows of numbers
    OUTPUT_CSV: an output file that will contain the sums\
"""])
    cli_parser = optparse.OptionParser(usage)
    cli_parser.add_option('-n', '--numprocs', type='int',
            default=NUM_PROCS,
            help="Number of processes to launch [DEFAULT: %default]")
    return cli_parser

class CSVWorker(object):
    def __init__(self, numprocs, infile, outfile):
        self.numprocs = numprocs
        self.infile = infile
        self.outfile = outfile
        self.in_csvfile = csv.reader(self.infile)
        self.inq = multiprocessing.Queue()
        self.outq = multiprocessing.Queue()

        self.pin = multiprocessing.Process(target=self.parse_input_csv, args=())
        self.pout = multiprocessing.Process(target=self.write_output_csv, args=())
        self.ps = [ multiprocessing.Process(target=self.sum_row, args=())
                        for i in range(self.numprocs)]

        self.pin.start()
        self.pout.start()
        for p in self.ps:
            p.start()

        self.pin.join()
        i = 0
        for p in self.ps:
            p.join()
            print "Done", i
            i += 1

        self.pout.join()
        self.infile.close()

    def parse_input_csv(self):
            """Parses the input CSV and yields tuples with the index of the row
            as the first element, and the integers of the row as the second
            element.

            The index is zero-index based.

            The data is then sent over inqueue for the workers to do their
            thing.  At the end the input process sends a 'STOP' message for each
            worker.
            """
            for i, row in enumerate(self.in_csvfile):
                row = [ int(entry) for entry in row ]
                self.inq.put( (i, row) )

            for i in range(self.numprocs):
                self.inq.put("STOP")

    def sum_row(self):
        """
        Workers. Consume inq and produce answers on outq
        """
        tot = 0
        for i, row in iter(self.inq.get, "STOP"):
                self.outq.put( (i, sum(row)) )
        self.outq.put("STOP")

    def write_output_csv(self):
        """
        Open outgoing csv file then start reading outq for answers
        Since I chose to make sure output was synchronized to the input there
        is some extra goodies to do that.

        Obviously your input has the original row number so this is not
        required.
        """
        cur = 0
        stop = 0
        buffer = {}
        # For some reason csv.writer works badly across processes so open/close
        # and use it all in the same process or else you'll have the last
        # several rows missing
        outfile = open(self.outfile, "w")
        self.out_csvfile = csv.writer(outfile)

        #Keep running until we see numprocs STOP messages
        for works in range(self.numprocs):
            for i, val in iter(self.outq.get, "STOP"):
                # verify rows are in order, if not save in buffer
                if i != cur:
                    buffer[i] = val
                else:
                    #if yes are write it out and make sure no waiting rows exist
                    self.out_csvfile.writerow( [i, val] )
                    cur += 1
                    while cur in buffer:
                        self.out_csvfile.writerow([ cur, buffer[cur] ])
                        del buffer[cur]
                        cur += 1

        outfile.close()

def main(argv):
    cli_parser = make_cli_parser()
    opts, args = cli_parser.parse_args(argv)
    if len(args) != 2:
        cli_parser.error("Please provide an input file and output file.")

    c = CSVWorker(opts.numprocs, args[0], args[1])

if __name__ == '__main__':
    main(sys.argv[1:])