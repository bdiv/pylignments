#!/usr/bin/env python

#imports
import argparse    # for argument parsing
import re          # for regex

#classes

class aligner:
    """ Alignes two sequences using Needleman-Wunsch, Smith&Waterman or Endfree alignment """
    # attributes
    def __init__(self, conf):
        self.config = conf        # store config
        self.m = len(conf.seq1)   # store length of sequence1
        self.n = len(conf.seq2)   # store length of sequence2

class config:
    """Generates the config out of a single file input"""
    # attributes
    path = ""      # path of config file
    algorithm = "" # algorithm to use
    match = 0      # matching score
    mismatch = 0   # mismatch score
    indel = 0      # indel score
    seq1 = ""      # sequence 1
    seq2 = ""      # sequence 2
    # constructor - already reads the file and fills everything
    def __init__(self, path):
        file = open(path, 'r')  # open config file in read mode
        text = file.read()      # read everything for convenience
        algo = re.search("alg:\s(sw|nw|ef|SW|NW|EF)", text)           # search the algorithm
        match = re.search("match:\s([\+|-][\d]+)", text)       # search the matching score
        mismatch = re.search("mismatch:\s([\+|-][\d]+)", text) # search the mismatching score
        indel = re.search("indel:\s([\+|-][\d]+)", text)       # search the indel score
        seq1 = re.search("seq1: ([A-Z]{1,200})", text)       # search sequence 1
        seq2 = re.search("seq2: ([A-Z]{1,200})", text)       # search sequence 2
    
        # check all the fields and safe them if they're ok    
        if(algo):
            self.algorithm = algo.group(1)
        else:
            raise Exception("Algorithm not defined in file " + path + ". Please insert 'alg: XX' statement, where XX is either sw, nw or ef.")
        
        if(match):
            self.match = int(match.group(1))
        else:
            raise Exception("Match score not defined in file " + path + ". Please insert 'match: yX', where y is either + or - and X is a number.")
        
        if(mismatch):
            self.mismatch = int(mismatch.group(1))
        else:
            raise Exception("Mismatch score not defined in file " + path + ". Please insert 'mismatch: yX', where y is either + or - and X is a number.")      
        
        if(indel):
            self.indel = int(indel.group(1))
        else:
            raise Exception("Indel score not defined in file " + path + ". Please insert 'indel: yX', where y is either + or - and X is a number.")
        
        if(seq1):
            self.seq1 = seq1.group(1)
        else:
            raise Exception("Sequence one not found in file " + path + ". Please insert 'seq1: S', where S can be a sequence with a length between one and 200 characters.")
        
        if(seq2):
            self.seq2 = seq2.group(1)
        else:
            raise Exception("Sequence one not found in file " + path + ". Please insert 'seq2: S', where S can be a sequence with a length between one and 200 characters.")
    
    # methods


###############################
#main##########################
###############################

# commandline argument configuration
parser = argparse.ArgumentParser()
parser.add_argument("file", help="config file, which specifies the sequences, scores and the algorithm to use")
args = parser.parse_args()
conf = config(args.file)



