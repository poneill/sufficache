import string
from utils import *

class TranscriptionFactor(object):
    """Represents the entries in matrix.dat"""
    def __init__(self, chunk):
        """Accepts a chunk of lines from matrix.dat and returns an
        object representation"""
        column_lines, other_lines = separate(matches_column, chunk)
        matrix = self.matrix_from_lines(column_lines)
        usable_lines = filter(self.usable_line,other_lines)
        self.attr_dict = {}
        for line in usable_lines:
            tlid = line[:2]
            attr = line[2:].strip()
            if not tlid in self.attr_dict:
                self.attr_dict[tlid] = [attr]
            else:
                self.attr_dict[tlid].append(attr)
        self.attr_dict["pssm"] = self.pssm_from_matrix(matrix)
        for key in self.attr_dict:
            val = self.attr_dict[key]
            exec("""storable_val = val if len(val) > 1 else val[0];self.{0} = storable_val""".format(key,val))
        
    def pssm_from_matrix(self, matrix,background_probs = (0.25,)*4):
        """Accept count matrix (as nested list) and return pssm (as nested list)"""
        def convert_column(col):
            return [safe_log2(c/p) for (c,p) in zip(normalize(col),background_probs)]
        return [convert_column(col) for col in matrix]

    def usable_line(self,line):
        bad_tlids = ['XX', #separates entries in chunk
                     '//', #separates chunks in file
                     'P0'] #signals beginning of PSWM
        return not any(line.startswith(tlid) for tlid in bad_tlids)

    def matrix_from_lines(self,lines):
        """Convert raw column lines into matrix.  Assumes all lines are
        column lines"""
        return [map(float,matches_column(line).groups()) for line in lines]

