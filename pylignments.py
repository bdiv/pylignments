#!/usr/bin/env python

#imports
import argparse    # for argument parsing
import re          # for regex
import numpy as np # for numpy arrays

#classes

class aligner:
    """ Alignes two sequences using Needleman-Wunsch, Smith&Waterman or Endfree alignment """
    # attributes
    def __init__(self, conf):
        self.config = conf        # store config
        self.m = len(conf.seq1)   # store length of sequence1
        self.n = len(conf.seq2)   # store length of sequence2
        

    def initNeedlemanWunsch(self):
        n = self.n
        m = self.m
        value = 0
        self.path = np.array([[[[0]*n]*m]*n]*m)
        self.matrix = np.array([[[0,0]]*n+1]*m+1)
        for x in range(0,m+1):
	    self.matrix[0][x] = [value,1]
            value = value + self.config.indel
	value = 0
	for x in range(0,n+1):
	    self.matrix[x][0] = [value,1]
            value = value + self.config.indel
    
    def needlemanWunsch(self):
        
    def needlemanWunsch_iterative(self):
        
        
    def needlemanWunsch_recursive(self, m, n):
        if(self.matrix[m-1][n-1][1] == 0):
            needlemanWunsch_recursive(m-1,n-1)   # calculate diagonal if not set
        if(self.matrix[m-1][n][1] == 0):
            needlemanWunsch_recursive(m-1,n)     # calculate left neighbor
        if(self.matrix[m][n-1][1] == 0):
            needlemanWunsch_recursive(m,n-1)     # calculate upper neighbor
	self.needlemanWunsch_calcScore(m,n)

     def needlemanWunsch_calcScore(m,n):
        # calculate all possible incoming scores
        pDiagonalScore = 0
        if(self.config.seq1[m-1] == self.config.seq2[n-1]):
            pDiagonalScore = self.matrix[m-1][n-1][0] + self.config.match
        else
            pDiagonalScore = self.matrix[m-1][n-1][0] + self.config.mismatch
        pUpperScore = self.matrix[m][n-1][0] + self.config.indel
        pLeftScore  = self.matrix[m-1][n][0] + self.config.indel
        #
        # now check scores and set the current score at (m,n) to the highest. (+ raise the set flag and save paths)
        if(pDiagonalScore >= pUpperScore and pDiagonalScore >= pLeftScore):
            self.matrix[m][n] = [pDiagonalScore,1]
            self.path[m][n][m-1][n-1] = 1
            self.path[m-1][n-1][m][n] = 1
        if(pUpperScore >= pDiagonalScore and pUpperScore >= pLeftScore):
            self.matrix[m][n] = [pUpperScore,1]
            self.path[m][n][m][n-1] = 1
            self.path[m][n-1][m][n] = 1
        if(pLeftScore >= pDiagonalScore and pLeftScore >= pUpperScore):
            self.matrix[m][n] = [pLeftScore,1]
            self.path[m][n][m-1][n] = 1
            self.path[m-1][n][m][n] = 1
        

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



