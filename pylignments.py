#!/usr/bin/env python

#imports
import argparse    # for argument parsing
import re          # for regex
import numpy as np # for numpy arrays

#classes

class aligner:
    """ Alignes two sequences using Needleman-Wunsch, Smith&Waterman or Endfree alignment """
    # constructor
    def __init__(self, conf):
        self.config = conf          # store config
        self.row = len(conf.seq1)   # store length of sequence1
        self.col = len(conf.seq2)   # store length of sequence2
        self.sequences = []         # will contain the aligned sequences
        # run the respective routine depending on the configured algorithm
        if(conf.algorithm == "nw" or conf.algorithm == "NW"):
            self.needlemanWunsch()
        elif(conf.algorithm == "ef" or conf.algorithm == "EF"):
            self.endFree()
        elif(conf.algorithm == "sw" or conf.algorithm == "SW"):
            self.smithWaterman()
        else:
            raise Exception("Algorithm \"" + conf.algorithm  + "\" not implemented")

    # needlemanWunsch part
    def initNeedlemanWunsch(self):
        col = self.col     # set local n 
        row = self.row     # set local m
        col1 = col+1       # stupid stuff bc of typeError during concatenation
        row1 = row+1       #  ...
        value = 0
        self.path = np.array([[[[0]*col1]*row1]*col1]*row1)    # create 4 dimmensional array to safe baths 
        self.matrix = np.array([[[0,0]]*col1]*row1)  # create our matrix
        for x in range(0,row1):                      # initialize first row
            self.matrix[x][0] = [value,1]            # 1 is the "set"-flag
            try:
                self.path[x][0][x+1][0] = 1
            except IndexError:
                1 == 1 # do nothing
            value = value + self.config.indel       # alter value for the next field
        value = 0                  
        for x in range(0,col1):                     # initialize first column
            self.matrix[0][x] = [value,1] 
            try:
                self.path[0][x][0][x+1] = 1
            except IndexError:
                1 == 1 # do nothing
            value = value + self.config.indel
    
    def needlemanWunsch(self):
        self.initNeedlemanWunsch()              # initialize marices (scores & path)
        self.recursive_calc(self.row,self.col)  # calculate actual scores
        self.backtrace()                        # backtrace, store aligned sequences in self.sequences and store best score
    
    def needlemanWunsch_iterative(self):
        print("not implemented")
        
    def recursive_calc(self, row, col):
        if(self.matrix[row-1][col-1][1] == 0):
            self.recursive_calc(row-1,col-1)   # calculate diagonal if not set
        if(self.matrix[row-1][col][1] == 0):
            self.recursive_calc(row-1,col)     # calculate left neighbor
        if(self.matrix[row][col-1][1] == 0):
            self.recursive_calc(row,col-1)     # calculate upper neighbor
        self.calcScore(row,col)

    def calcScore(self,row,col):
        # calculate all possible incoming scores
        pDiagonalScore = 0
        if(self.config.seq1[row-1] == self.config.seq2[col-1]):
            pDiagonalScore = self.matrix[row-1][col-1][0] + self.config.match
        else:
            pDiagonalScore = self.matrix[row-1][col-1][0] + self.config.mismatch
        pUpperScore = self.matrix[row][col-1][0] + self.config.indel
        pLeftScore  = self.matrix[row-1][col][0] + self.config.indel
        #
        # now check scores and set the current score at (m,n) to the highest. (+ raise the set flag and save paths)
        if(pDiagonalScore >= pUpperScore and pDiagonalScore >= pLeftScore):
            # if the algorithm is smith waterman and the score is below 0
            if((self.config.algorithm == "sw" or self.config.algorithm == "SW") and pDiagonalScore < 0):
                self.matrix[row][col] = [0,1]              # set score to 0, no path is saved
            else:
                self.matrix[row][col] = [pDiagonalScore,1] # else save score and path
                self.path[row][col][row-1][col-1] = 1
                self.path[row-1][col-1][row][col] = 1
        if(pUpperScore >= pDiagonalScore and pUpperScore >= pLeftScore):
            # if the algorithm is smith waterman and the score is below 0
            if((self.config.algorithm == "sw" or self.config.algorithm == "SW") and pUpperScore < 0):
                self.matrix[row][col] = [0,1]              # set score to 0, no path is saved
            else:
                self.matrix[row][col] = [pUpperScore,1]    # else save score and path
                self.path[row][col][row][col-1] = 1
                self.path[row][col-1][row][col] = 1
        if(pLeftScore >= pDiagonalScore and pLeftScore >= pUpperScore):
            if((self.config.algorithm == "sw" or self.config.algorithm == "SW") and pLeftScore < 0):
                self.matrix[row][col] = [0,1]
            else:
                self.matrix[row][col] = [pLeftScore,1]
                self.path[row][col][row-1][col] = 1
                self.path[row-1][col][row][col] = 1
    
    # end free
    def initEndFree(self):
        col = self.col     # set local n
        row = self.row     # set local m
        col1 = col+1       # stupid stuff bc of typeError during concatenation while initializing the matrices below
        row1 = row+1       #  ...
        self.path = np.array([[[[0]*col1]*row1]*col1]*row1)    # create 4 dimmensional array to safe baths
        self.matrix = np.array([[[0,0]]*col1]*row1)  # create our matrix
        for x in range(0,row1):                      # initialize first row
            self.matrix[x][0] = [0,1]                # 1 is the "set"-flag
        for x in range(0,col1):                      # initialize first column
            self.matrix[0][x] = [0,1]
    
    def endFree(self):
        self.initEndFree()                       # initialize score matrix and path matrix
        self.recursive_calc(self.row, self.col)  # calculate scores in the matrix, this uses the self.config.algorithm setting to determine the MO
        self.backtrace()                         # backtrace uses the algorithm value to determine it's MO
  
    # smith waterman
    def initSmithWaterman(self):
        self.initEndFree()    # uses the same initialization as end free does
        
    def smithWaterman(self):
        self.initSmithWaterman()                 # initialize score and path matrices
        self.recursive_calc(self.row, self.col)  # calculate scores recursively
        self.backtrace()                         # backtrace optimal paths, store alignments and score
  
    # general methods        
    # print the matrices with predefined style characters from config
    def printMatrix(self, backtrace):
        # get the configured characters for visulization
        right = self.config.charPathRight
        down = self.config.charPathDown
        diag = self.config.charPathDiag
        con = self.config.charSpacerCon
        hSpacer = self.config.charHSpacer
        vSpacer = self.config.charVSpacer
        head = "  "+"{:>6s}"*len(self.config.seq2) # head string
        print(head.format(*self.config.seq2))      # print head string
        x = len(self.config.seq2)+1                
        print(" " + con + (hSpacer*5 + con )*x)                     # print spacer, x is the number of fields
        for row in range(0,len(self.config.seq1)+1): # for every line   
            line1 = " " + vSpacer                  # line1 contains the scores and path flags and indicators
            try: 
                line2 = self.config.seq1[row] + con  # line2 is the spacer and contains path indicators
            except IndexError:                     # if we're in the last line
                line2 = " " + con
            for col in range(0, len(self.config.seq2)+1):                  # for every column in the line
                if(backtrace == 1 or self.bpath[row][col] == self.config.charBacktraceEnd):
                    line1 = line1 + "{:4d}".format(self.matrix[row][col][0]) + self.bpath[row][col] # score field + backtrace path
                else:
                    line1 = line1 + "{:4d} ".format(self.matrix[row][col][0])    # score field
                try:                                                     
                    if(self.path[row][col][row][col+1]==1):              # if path exist
                        line1 = line1 + right
                    else:
                        line1 = line1 + vSpacer
                except IndexError:
                    line1 = line1 + vSpacer                              # if we're at the last column...
                line2 = line2 + hSpacer*2                                # first part of the spacer
                try: 
                    if(self.path[row][col][row+1][col]==1):                      # middle part of the spacer, the path indicator
                        line2 = line2 + down
                    else:
                        line2 = line2 + hSpacer
                except IndexError:
                    line2 = line2 + hSpacer                              # if path.[row][col][row+1][col] is outside of our matrix
                line2 = line2 + hSpacer*2                                # third part of the spacer
                try:                                                     # last part, it's a path indicator again. Prints a + if we're out of bounds
                    if(self.path[row][col][row+1][col+1]==1):
                        line2 = line2 + diag
                    else: 
                        line2 = line2 + con
                except IndexError:
                    line2 = line2 + con
            print(line1)
            print(line2)
    
    # print the spread sheet version shown in the example
    def printSpreadSheet(self):
        print("Forward part:")
        print("")
        self.printMatrix(0)
        print("")
        print("")
        print("Backward part:")
        print("")
        align.printMatrix(1)
        print("")
        print("")
        print("Alignments (Score: " + str(align.score) + "):")
        for ali in align.sequences:
            print("")
            print(ali[1])
            print(ali[0])
        print("")
    
    def backtrace(self):
        # base method for backtracing, will call a recursive routine
        col1 = self.col + 1                            # for convenience
        row1 = self.row + 1                            # for convenience
        self.bpath = np.array([[" "]*col1]*row1)       # stores the backtrace path as characters
        if(self.config.algorithm == "nw" or self.config.algorithm == "NW"):
            # if we use needleman wunsch
            self.backtrace_recursive(self.row, self.col, ("",""))            # start recursive routine
            self.bpath[self.row][self.col] = self.config.charBacktraceEnd    # mark the end of the path
            self.score = self.matrix[self.row][self.col][0]                  # save the score
        elif(self.config.algorithm == "sw" or self.config.algorithm == "SW"):
            # if we use smith waterman
            s = (0, [(0,0)])                                         # will contain the scores and and array containing the corresponding points
            # loop through the matrix and find the highest scores
            for row in range(0, self.row + 1):
                for col in range(0, self.col + 1):
                    if(self.matrix[row][col][0] > s[0]):             # if we found a higher score
                        s = (self.matrix[row][col][0],[(row, col)])  # overwrite s with new values
                    elif(self.matrix[row][col][0] == s[0]):          # if the score is equal to our stored one
                        s[1].append((row, col))                      # append to the list of points
            self.score = s[0]                                        # store score
            for point in s[1]:                                       # loop through points
                self.backtrace_recursive(point[0], point[1], ("",""))# start backtrace for every point with the highest score
            for point in s[1]:
                self.bpath[point[0]][point[1]] = self.config.charBacktraceEnd # mark the end of every path
        elif(self.config.algorithm == "ef" or self.config.algorithm == "EF"):
            # if we use end free
            s = (0, [(0,0)])                                             # used to store the highest score and the respective points
            for row in range(0, self.row):                               # loop through all points of the right boundary 
                if(self.matrix[row][self.col][0] > s[0]):                # if we found a higher score
                    s = (self.matrix[row][self.col][0],[(row, self.col)])# overwrite the current s
                elif(self.matrix[row][self.col][0] == s[0]):             # if we found an equal score
                    s[1].append((row, self.col))                         # append the point
            for col in range(0, self.col+1):                             # loop through the points of the lower boundary
                if(self.matrix[self.row][col][0] > s[0]):                # if we found a higher score
                    s = (self.matrix[self.row][col][0], [(self.row, col)])# overwrite s
                elif(self.matrix[self.row][col][0] == s[0]):             # if we found an equal score
                    s[1].append((self.row, col))                         # append point to s
            self.score = s[0]                                            # store score
            for point in s[1]:                                           # loop through all points we found
                self.backtrace_recursive(point[0],point[1], ("",""))     # backtrace them
            for point in s[1]:                                           # loop through points again
                self.bpath[point[0]][point[1]] = self.config.charBacktraceEnd# mark the ends
        else:
            raise Exception("Cannot backtrace because of undefined algorithm")   # this should never happen bc we check the algorithm while parsing the config file

    def backtrace_recursive(self,row,col, seq):
        # recursive backtrace routine
        try:
            seq1 = seq[0] # for convenience
            seq2 = seq[1] # for convenience
            self.bpath[row][col] = self.config.charBacktracePart # mark the current path field
            upPath = self.path[row][col][row-1][col]             # contains 0 or 1, depending on the existence of a path
            leftPath = self.path[row][col][row][col-1]           # contains 0 or 1, depending on the existence of a path
            diagPath = self.path[row][col][row-1][col-1]         # contains 0 or 1, depending on the existence of a path
            # print(str(seq) + " row: " + str(row) + " col: " + str(col) + " up: " + str(upPath) + " left: " + str(leftPath) + " diag: " + str(diagPath))
            # we have different conditions for exit
            if((self.config.algorithm == "nw" or self.config.algorithm == "NW") and row == 0 and col == 0): # if we reach (0,0) in needlemanWunsch
                self.sequences.append(seq)
                return            
            if((self.config.algorithm == "ef" or self.config.algorithm == "EF") and (row == 0 or col == 0)):# if we reach the upper or left boundary in end free
                self.sequences.append(seq)
                return
            if((self.config.algorithm == "sw" or self.config.algorithm == "SW") and (upPath == 0 and leftPath == 0 and diagPath == 0 )):
                # if we have nowhere to go in smith waterman 
                # note: the scoring algorithm (calcscore_recursive())marks the path. If a score gets below 0 it is set to 0 and the path ends.
                # we could've also just checked the score instead
                self.sequences.append(seq)
                return
            if(leftPath == 1):                                  # if we can go left
                # print("left")
                s1 = "-" + seq1                                 # append "-"
                s2 = self.config.seq2[col-1] + seq2             # append next character to the sequence
                self.backtrace_recursive(row, col-1, (s1,s2))   # go to the left field
            if(diagPath == 1):                                  # if we can go diagonal
                # print("diag")
                s1 = self.config.seq1[row-1] + seq1             # align char of seq1 with...
                s2 = self.config.seq2[col-1] + seq2             # ...char of seq2
                self.backtrace_recursive(row-1, col-1,(s1,s2))  # go diagonal
            if(upPath == 1):                                    # if we can go up
                # print("up")
                s1 = self.config.seq1[row-1] + seq1             # append next character to the sequence
                s2 = "-" + seq2                                 # appen "-"
                self.backtrace_recursive(row-1,col, (s1,s2))    # go up
            # up = self.matrix[row-1][col][0]
            # left = self.matrix[row][col-1][0]
            # diag = self.matrix[row-1][col-1][0]
            '''if(upPath == 1 and leftPath == 1 and diagPath == 1):
                if(up >= left and up >= diag):
                    self.backtrace_recursive(row-1,col)
                if(left >= up and left >= diag):
                    self.backtrace_recursive(row, col-1)
                if(diag >= up and diag >= left):
                    self.backtrace_recursive(row-1, col-1)
            elif(upPath == 1 and leftPath == 1 and diagPath == 0):
                if(up >= left):
                    self.backtrace_recursive(row-1,col)
                if(left >= up):
                    self.backtrace_recursive(row, col-1)
            elif(upPath == 1 and leftPath == 0 and diagPath == 0):
                self.backtrace_recursive(row-1,col)
            elif(upPath == 0 and leftPath == 1 and diagPath == 1):
                if(left >= diag):
                    self.backtrace_recursive(row, col-1)
                if(diag >= left):
                    self.backtrace_recursive(row-1, col-1)
            elif(upPath == 0 and leftPath == 0 and diagPath == 1):
                self.backtrace_recursive(row-1, col-1)
            elif(upPath == 0 and leftPath == 1 and diagPath == 0):
                self.backtrace_recursive(row, col-1)
            elif(upPath == 1 and leftPath == 0 and diagPath == 1):
                if(up >= diag):
                    self.backtrace_recursive(row-1,col)
                if(diag >= up):
                    self.backtrace_recursive(row-1, col-1)
            else: 
                print("fuck my life")'''
            
        except IndexError:
            self.sequences.append(seq)
            return

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
    charPathDown = "|"      # for the matrix display
    charPathRight = "-"     # for the matrix display
    charPathDiag = "\\"     # for the matrix display
    charBacktracePart = "*" # for the matrix display
    charBacktraceEnd = "#"  # for the matrix display
    charVSpacer = "|"       # for the matrix display
    charHSpacer = "-"       # for the matrix display
    charSpacerCon = "+"
    # constructor - already reads the file and fills everything
    def __init__(self, path):
        file = open(path, 'r')                                 # open config file in read mode
        text = file.read()                                     # read everything for convenience
        file.close()                                           # close file
        # using regex because regex are cool
        algo = re.search("alg:\s(sw|nw|ef|SW|NW|EF)", text)    # search the algorithm
        match = re.search("match:\s([\+|-][\d]+)", text)       # search the matching score
        mismatch = re.search("mismatch:\s([\+|-][\d]+)", text) # search the mismatching score
        indel = re.search("indel:\s([\+|-][\d]+)", text)       # search the indel score
        seq2 = re.search("seq1: ([A-Z]{1,200})", text)         # search sequence 1
        seq1 = re.search("seq2: ([A-Z]{1,200})", text)         # search sequence 2
        #
        # check all the fields and save them if they're ok    
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
conf = config(args.file)             # create config object (parses config file)
align = aligner(conf)                # create aligner. this will start the alignment based on the provided configs
#print(align.matrix)
align.printSpreadSheet()             # print out neat little spread sheet 
