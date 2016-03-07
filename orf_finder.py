__author__ = 'lakinsm'

import os.path


class FastaParser:
    """This object takes as input a FASTA filepath and constructs an iterable that outputs
    (header, sequence) tuples.  Only one header/sequence will be held in memory at a time using this method.
    """

    def __init__(self, filepath):
        """
        constructor
        @param filepath: filepath to the input FASTA file.
        """
        if os.path.isfile(filepath):
            self.fasta_file = filepath
        else:
            raise ValueError("Parameter filepath must be a FASTA file")
        self.current_line = None
        self.current_header = None

    def __iter__(self):
        return self

    def _iterate(self):
        # Skip all non-fasta entry lines
        while True:
            fasta_line = self.fasta_file.readline()
            if fasta_line == "":
                return
            if fasta_line.startswith(">"):
                break

        while True:
            if fasta_line[0] != ">":
                self.fasta_file.close()
                raise ValueError("FASTA records should begin with '>'")
            self.current_header = fasta_line[1:].rstrip()
            fasta_lines = []
            fasta_line = self.fasta_file.readline()
            iline = self.fasta_file.tell() # Seek to here if new entry found in the next section
            while True:
                if not fasta_line:
                    break
                if fasta_line.startswith(">"):
                    break
                fasta_lines.append(fasta_line.rstrip())
                fasta_line = self.fasta_file.readline()
            self.fasta_file.seek(iline)
            return self.current_header, "".join(fasta_lines).replace(" ", "").replace("\r", "")
        self.fasta_file.close()
        assert False, "Should not reach this line"

    def next(self):
        if type(self.fasta_file) is str:
            self.fasta_file = open(self.fasta_file, "r")
        value = self._iterate()
        if not value:
            self.fasta_file.close()
            raise StopIteration()
        else:
            return value


class OrfFinder:
    def __init__(self, entry, transl_table=1):
        self.start = 0
        self.stop = 0
        self.loc = None

a = FastaParser('/home/lakinsm/Downloads/Lakin_Homework2.latex')


