import re, string, math
import TranscriptionFactor

def matches_accession_number(line):
    """Return an re.match object for the accession number pattern """
    return re.search(r'^AC\s+([A-Z0-9]+)', line)

def matches_column(line):
    """Return an re.match object for the column pattern"""
    regexp = """
^[0-9]+     #Begins with column number
\s+         #Followed by some whitespace
            #and now numbers corresponding to the base count at that column:
([.0-9]+)\s+ #number of As
([.0-9]+)\s+ #Cs
([.0-9]+)\s+ #Gs
([.0-9]+)    #Ts
"""
    return re.search(regexp, line, re.VERBOSE)

class TransfacTable(object):
    """Represents matrix.dat"""
    def __init__(self, filename):
        """Accept a filename (matrix.dat) and return a
        dictionary representation thereof"""
        lines = open(filename).readlines()
        accession_chunks = split_on(lines, matches_accession_number)[1:]
        self.entries = map(self.parse_accession_chunk, accession_chunks)

    def parse_accession_chunk(self, chunk):
        """Accept an accession chunk represented as a list of entries
        corresponding to lines of the input file and return a PSSM and
        a dictionary associating the two letter prefixes of the file
        with the line contents"""
        return TranscriptionFactor(chunk)
    
