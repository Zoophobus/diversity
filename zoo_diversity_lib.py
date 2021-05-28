#!/usr/bin/env python3.7

import sys,os,re#,glob
#import argparse


class Base(object):
    def __init__(self,n,nextBase=None):
        self.data = n
        self.nxt = nextBase
        self.prv = None
    def get_data(self):
        return self.data
    def set_data(self,newData):
        self.data = newData
    def next(self):
        return self.nxt
    def prev(self):
        return self.prv
    def set_next(self,new_base):
        self.nxt = new_base
    def set_prev(self,new_base):
        self.prv = new_base

class Reads(object):
    ''' was formerly the class Sequences '''
    def __init__(self,data=None):
        self.sequence = ""
        self.begins = 0
        self.finishes = 1
        if data is not None:
            self.finishes = len(data)
            self.add(data)

## TODO the following has been adjusted so it can be used for larger sequence alignments
## it may be required to change this at a later date for shorter reads again
    def add(self,seq): # this adds the sequence to the Sequence object, it also provides
        self.sequence = seq # an integer that is set to the actual beginning of the 
        self.begins = 0
        self.finishes = len(self.sequence)

    def get(self,seq):
        if self.begins > seq or self.finishes < seq: # this condition might need
            return "-"     # to be changed at some point  
        return self.sequence[seq]
    def start(self):
        return self.begins
    def end(self):
        return self.finishes
    def __str__(self):
        return self.sequence
'''
class Sequence(object):
    def __init__(self,data=None):
        self.head = None
        self.tail = None
        self.max_length = window_length
        self.count = 0
        self.seq_ref = None
        if data is not None:
            self.insert(data)
    def insert(self,seq): # setting up the window from the very first character
        self.seq_ref = seq
        new_base = Base(self.seq_ref.get(0)) 
        new_base.set_next(self.head)
        self.head = new_base
        self.tail = new_base
        self.count += 1
        for i in range(1,self.max_length,1): # setup the initial stretch of the window 
            new_base = Base(self.seq_ref.get(i))
            new_base.set_next(self.head)
            self.head.set_prev(new_base)
            self.head = new_base
            self.count += 1
    def size(self):
        return self.max_length
    def add_base(self): # adds the next base to the window, while simultaneously
        discard_base = self.tail # removing the first base
        self.tail = self.tail.prev()
        self.tail.set_next(None)
        new_base = Base(self.seq_ref.get(self.count))
        new_base.set_next(self.head)
        self.head.set_prev(new_base)
        self.head = new_base
        self.count += 1
        return discard_base
    def last(self):
        return self.head
    def next(self,base=None): # iterates through the doubly linked list
        if base == None:
            return self.head
        return base.next()
    def prev(self,base=None):
        if base == None:
            return self.tail
        return base.prev()
    def start_from(self):
        return self.seq_ref.start()
    def end_from(self):
        return self.seq_ref.end()
    def position(self):
        return (self.count - (self.max_length / 2))
    def __str__(self):
        bIter = self.prev(None)
        currentStr = ""
        while bIter != None:
            currentStr += bIter.get_data()
            bIter = self.prev(bIter)
        return currentSt
'''
class Sequence(object):
    def __init__(self,seq_name=None,data=None,length=66):
        self.head = None
        self.tail = None
        self.max_length = length
        self.count = 0
        self.seq_ref = None
        self.name = None
        if seq_name is not None:
            self.name = seq_name
        if data is not None:
            self.insert(data)
    def insert(self,seq): # setting up the window from the very first character
        self.seq_ref = seq
        new_base = Base(self.seq_ref.get(0)) 
        new_base.set_next(self.head)
        self.head = new_base
        self.tail = new_base
        self.count += 1
        for i in range(1,self.max_length,1): # setup the initial stretch of the window 
            new_base = Base(self.seq_ref.get(i))
            new_base.set_next(self.head)
            self.head.set_prev(new_base)
            self.head = new_base
            self.count += 1
    def size(self):
        return self.max_length
    def addAndRetrieve(self,item):
        if item == None:
            new_base = Base(self.seq_ref.get(self.count))
        else:
            new_base = Base(item)
        discard_base = self.tail # removing the first base
        self.tail = self.tail.prev()
        self.tail.set_next(None)
        self.head.set_prev(new_base)
        self.head = new_base
        self.count += 1
        return discard_base

    def add_base(self,item=None): # adds the next base to the window, while simultaneously
        if item != None and self.count < self.max_length:
            new_base = Base(item)
            new_base.set_next(self.head)
            self.head.set_prev(new_base)
            self.head = new_base
            self.counter += 1
            return None
        discard_base = self.addAndRetrieve(item)
        return discard_base
    def last(self):
        return self.head
    def next(self,base=None): # iterates through the doubly linked list
        if base == None:
            return self.head
        return base.next()
    def prev(self,base=None):
        if base == None:
            return self.tail
        return base.prev()
    def start_from(self):
        return self.seq_ref.start()
    def end_from(self):
        return self.seq_ref.end()
    def position(self):
        return (self.count - (self.max_length / 2))
    def __str__(self):
        bIter = self.prev(None)
        currentStr = self.name + "\n"
        while bIter != None:
            currentStr += bIter.get_data()
            bIter = self.prev(bIter)
        return currentStr

class Windows(object):
    def __init__(self,window_length):
        self.windows_start_at = 0
        self.windows_end_at = 0
        self.sub_windows = []
        self.no_windows = 0
        self.wattersons = True
        self.diversity = 0
        self.seq_diversity = []
        self.consensus = []
        self.location = 0
        self.thisMismatch = 0 
        self.divergence = 0
        self.window_length = window_length
        self.denominator = 0.0
    def size(self):
        return self.no_windows
    def insert(self,name,seq):
        self.sub_windows.append(Sequence(name,seq))
        self.no_windows += 1
    def theta(self,version):
        self.wattersons = version
    def set_denominator(self):
        if self.wattersons:
            for i in range(1,self.no_windows):
                self.denominator += 1.0/i
        else:
            self.denominator = (self.no_windows * (self.no_windows -1) /2)
    def valid_data(self,aChar):
        if aChar == '-' or aChar == "?" or aChar == "N" or aChar == "n":
            return False 
        return True

    def set_windows(self,other=None,shorah=False,end=1):
        if not shorah:
            self.windows_start_at = 0
            self.windows_end_at = end
            return
        ''' first we find the largest region that all the sequences cover '''
        k = self.sub_windows[0].start_from()
        e = self.sub_windows[0].end_from()
        for i in range(1,self.no_windows):
            if k < self.sub_windows[i].start_from():
                k = self.sub_windows[i].start_from()
            if e > self.sub_windows[i].end_from():
                e = self.sub_windows[i].end_from()
        ''' Then if we also need to compare this to another group we include them in the above search '''
        if not other is None:
            for i in range(other.no_windows):
                if k < other.sub_windows[i].start_from():
                    k = other.sub_windows[i].start_from()
                if e > other.sub_windows[i].end_from():
                    e = other.sub_windows[i].end_from()
            other.windows_start_at = k
            other.windows_end_at = e
            for i in range(other.no_windows):
                l = 0
                while l < other.windows_start_at:
                    other.sub_windows[i].add_base()
                    l += 1
        self.windows_start_at = k
        self.windows_end_at = e
        ''' finally we get this group of sequences to move up to the first base that
            is covered by all of the available sequences '''
        for i in range(self.no_windows):
            l = 0
            while l < self.windows_start_at:
                self.sub_windows[i].add_base()
                l += 1

    def window_stats(self,thisPosition,otherPosition,I=None,other=None): 
#        thisMismatch = 0  ## NOW THESE FOLLOWING TWO VARS ARE CLASS
#        divergence = 0 ## VARIABLES
        seqDivergence = 0.0
        ####################################################
        ### changing method structure with self.sub_windows.size() loop
        ####################################################
        thisBaseFreq = {}
        localMismatch = 0
        missingData = re.compile(r'-\?')
        if other is None:   # if we just have to do the one!
        # should set this section to another method then we can use it three times!!
            for k in range(self.no_windows):
                thisPosition[k] = self.sub_windows[k].prev(thisPosition[k])
                if not thisPosition[k].get_data() in thisBaseFreq:
                    thisBaseFreq[thisPosition[k].get_data()] = 1
                else:
                    thisBaseFreq[thisPosition[k].get_data()] += 1
                if k > 0: # skip the first element as this is comparing the same position 
                    for l in range(k):
                        if thisPosition[k].get_data() != thisPosition[l].get_data() and self.valid_data(thisPosition[k].get_data()) and self.valid_data(thisPosition[l].get_data()):
                            if self.wattersons:
                                localMismatch = 1
                                break
                            else:
                                localMismatch += 1
            self.consensus.append(self.get_consensus(thisBaseFreq))
            self.thisMismatch += localMismatch
#            if I >= (self.window_length /2):
            self.location += 1
            if self.thisMismatch == 0:
                self.seq_diversity.append(0)
                self.diversity = 0
            else:
                self.seq_diversity.append((self.thisMismatch /(self.denominator * (I+1))))
                self.diversity = self.thisMismatch / self.denominator
#            else:
#                self.seq_diversity.append((self.thisMismatch /(self.denominator * (I+1))))
#                self.location += 1
#                self.diversity = self.thisMismatch / self.denominator
            return self.diversity
        ''' Otherwise we have to do both'''
        otherBaseFreq = {}
        for k in range(self.no_windows):
            otherPosition[k] = other.sub_windows[k].prev(otherPosition[k])
            thisPosition[k] = self.sub_windows[k].prev(thisPosition[k])
            if not thisPosition[k].get_data() in thisBaseFreq:
                thisBaseFreq[thisPosition[k].get_data()] = 1
            else:
                thisBaseFreq[thisPosition[k].get_data()] += 1
            if not otherPosition[k].get_data() in otherBaseFreq:
                otherBaseFreq[otherPosition[k].get_data()] = 1
            else:
                otherBaseFreq[otherPosition[k].get_data()] += 1
            if k > 0: # skip the first element as this is comparing the same position 
                for l in range(k):
                    if thisPosition[k].get_data() != thisPosition[l].get_data() and self.valid_data(thisPosition[k].get_data()) and self.valid_data(thisPosition[l].get_data()):
                        if self.wattersons:
                            self.thisMismatch = 1
                        else:
                            self.thisMismatch += 1
                    if otherPosition[k].get_data() != otherPosition[l].get_data() and self.valid_data(thisPosition[k].get_data()) and self.valid_data(thisPosition[l].get_data()):
                        if self.wattersons:
                            other.thisMismatch = 1
                        else:
                            other.thisMismatch += 1
        self.consensus.append(self.get_consensus(thisBaseFreq))
        other.consensus.append(other.get_consensus(otherBaseFreq))
        if self.consensus[I] != other.consensus[I]:
            self.divergence += 1
        if I >= (self.window_length / 2):
            self.seq_diversity.append((self.thisMismatch /(self.denominator * I)))
            other.seq_diversity.append((other.thisMismatch /(self.denominator * I)))
            seqDivergence = (self.divergence / i)
            self.location += 1
            other.location += 1
            self.diversity = self.thisMismatch / self.denominator
            other.diversity = other.thisMismatch / self.denominator
        return seqDivergence
    def slide(self,other=None):
        divergence = 0
#       if not other is None:
#           discard = other.slide()
        if self.location >= (self.windows_end_at + 1 - self.windows_start_at): # statistics
            return None
        position = []
        latest = []
        for k in range(self.no_windows):
            position.append(None)
            latest.append(None)
        mismatch = 0
        by = 1
        while by > 0:
            if not other is None:
                discard = other.slide()
            baseFreq = {}
            for i in range(self.no_windows):
                position[i] = self.sub_windows[i].add_base() # the base that is removed
                latest[i] = self.sub_windows[i].last() # the new base
                if not latest[i].get_data() in baseFreq:
                    baseFreq[latest[i].get_data()] = 1
                else:
                    baseFreq[latest[i].get_data()] += 1
                if i > 0:
                    for j in range(i):
                        if position[i].get_data() != position[j].get_data() and self.valid_data(position[i].get_data()) and self.valid_data(position[j].get_data()):
                            if self.wattersons:
                                mismatch = -1
                            else:
                                mismatch -= 1
                        if (self.location + (self.window_length / 2)) < self.windows_end_at - self.windows_start_at and latest[i].get_data() != latest[j].get_data() and self.valid_data(latest[i].get_data()) and self.valid_data(latest[j].get_data()):
                        # if we have reached the end of the sequence we don't want to add any more mismatches
                            if self.wattersons:
                                mismatch = 1
                            else:
                                mismatch += 1
            by -= 1
            self.consensus.append(self.get_consensus(baseFreq)) # TODO change this
            self.diversity = self.diversity + (mismatch / self.denominator)
            divisor = self.window_length
            if self.location + (self.window_length/2) > self.windows_end_at:
                divisor = ((self.window_length/2) + (self.windows_end_at - self.location))
            self.seq_diversity.append(self.diversity/divisor)
            self.location += 1
            divergentLoss = int(self.location - (self.window_length/2))
            if not other is None:
                if self.consensus[divergentLoss] != other.consensus[divergentLoss]:
                    divergence -= 1 
                if self.consensus[self.location] != other.consensus[other.location]:
                    divergence += 1 
        if other is None:
            return 0
        return divergence / divisor # returning the consensus base
    def get_consensus_string(self): 
        return self.consensus
    def get_consensus(self,consensusMap):
        score = 0
        baseChar = ""
        for key,val in consensusMap.items():
            if val > score:
                score = val
                baseChar = key
        return baseChar
    def get_diversity(self):
        return self.seq_diversity
    def get_average_diversity(self):
        mid = self.windows_end_at - self.windows_start_at
        return self.seq_diversity[mid]

    def __str__(self):
        theStr = self.consensus()
        theStr += "\n"
        return theStr

class Alignment(object):
    # TODO need to modify this for the other method argument to be None
    def __init__(self,input_file,window_length=66):
        self.input_name = input_file
        self.no_sequences = 0
        self.sequences = []
        self.windows = Windows(window_length)
        self.position = 0
        self.max_length = window_length
        self.divergence = []
        self.required_coverage = 0
        self.required_posterior_prob = 0.95
        self.shorah = False

    def add(self,name,seq): # adds a sequence to the alignment, including to the array of
        exp = re.compile(r'.*?posterior=(\S+)')
        prob = re.match(exp,name)
        ave_cov = re.compile(r'.*?ave_reads=(\S+)')
        cov_g = re.match(ave_cov,name)
        if not self.shorah or float(prob.group(1)) > self.required_posterior_prob and float(cov_g.group(1)) > self.required_coverage:
            self.sequences.append(Reads(seq)) # Sequences, as well as the
            self.windows.insert(name,self.sequences[self.no_sequences])
            self.no_sequences += 1 # to the windows

    def set_shorah_input(self,shorah=False):
        self.shorah = shorah

    def get(self,seq):
        assert (seq < self.no_sequences), "Attempting to retrieve unavailable sequences"
        return self.sequences[seq]

    def set_Wattersons(self,value=True):
        self.windows.theta(value)

    def set_minimum_coverage(self,value=0):
        self.required_coverage = value

    def set_required_posterior(self,value=0):
        self.required_posterior_prob = value

    def set_required_no_haplotypes(self,value=2):
        self.minimum_haplotypes = value

    def start(self,other=None):
        if self.no_sequences < self.minimum_haplotypes:
            return
        tmpSeq = []
        groupOne = []
        groupTwo = []
        otherWindows = None
        if not other is None:
            otherWindows = other.windows
        self.windows.set_windows(otherWindows,self.shorah,self.sequences[0].end())
        self.windows.set_denominator()
        for i in range(self.windows.size()):
            groupOne.append(None)
            groupTwo.append(None)
        for i in range(self.max_length): ## TODO checking the code for different usage here
            value = self.windows.window_stats(groupOne,groupTwo,i,otherWindows)
            if i >= int(self.max_length/2):
#                tmpSeq.append(self.windows.window_stats(groupOne,groupTwo,i,otherWindows))
                tmpSeq.append(value)
        for i in range(int(self.max_length/2)):
            self.divergence.append(tmpSeq[i])
        divergenceIncr = self.windows.slide(otherWindows)
        if divergenceIncr is None: # this means more than half of the sequence of one
            return                  # haplotype is a gap
        self.divergence.append(self.divergence[-1]+divergenceIncr)
        while divergenceIncr != None:
            divergenceIncr = self.windows.slide(otherWindows)
            if not divergenceIncr is None:
                self.divergence.append(self.divergence[-1]+divergenceIncr)
    def output(self,fileOutstream,other=None):
#        self.print_header(fileOutstream,other)
        if self.no_sequences < self.minimum_haplotypes:
            if other is None:
                outputRow = self.input_name + ", " + "NA" + "\n"
                fileOutstream.write(outputRow)
                return
        consensusA = self.get_consensus()
        diversityA = self.get_diversity()
        if not other is None:
            consensusB = other.get_consensus()
            diversityB = other.get_diversity()
            for i in range(len(self.divergence)):
                outputRow = consensusA[i] + ", " + "{0:.4f}".format(diversityA[i]) + ", " + consensusB[i] + ", " + "{0:.4f}".format(diversityB[i]) + ", " + "{0:.4f}".format(self.divergence[i]) + "\n"
                fileOutstream.write(outputRow)
#        for i in range(len(self.divergence)):
#           outputRow = consensusA[i] + ", " + "{0:.4f}".format(diversityA[i]) + "\n"
        else:
            diversityA = self.windows.get_average_diversity()
            outputRow = self.input_name + ", " + "{0:.4f}".format(diversityA) + "\n"
            fileOutstream.write(outputRow)

    def __str__(self):
        theSeq = self.windows.get_consensus() + "\n" + self.windows.get_diversity()
        return theSeq
    def get_consensus(self):
        return self.windows.get_consensus_string()
    def get_diversity(self):
        return self.windows.get_diversity()
    def compare(self,other,fileOutStream):
        if not isinstance(other,Alignment):
            print("wrong object type provided to Alignment.compare(self,other)")
        fileOutStream.writelines(self.windows.compare(other))
    def get_divergence(self):
        return self.divergence
    def print_header(self,fileOutstream,other=None):
        if other is None:
            fileOutstream.write("Consensus, diversity\n")
            return
        fileOutstream.write("Consensus-A, diversity-A, Consensus-B, diversity-B, divergence\n") 





#seqA1 = "--CGTACAGTTACAAAGTTTGTTACTTGTCGTTACGGGCCG"
#seqA2 = "---GTACAGTTACAAAGTATGTTACTTGTCGTTACGGGCCG"
#seqA3 = "--GGTACAGATACAAAGTTTGTAACTTGTCGTTACGCGCCG"


#seqB1 = "ACCGGGCCGAGACAAAGTTTGTAACTTGTCGTTACGGGC--"
#seqB2 = "--CGGGCCGAGTCAATGTTTGTTACTTGTCGTTACCGGCC-"
#seqB3 = "-CCGGGCCGAGACAAAGTTGGTAACTTGTCGTTACGGGCC-"



#alignmentA.start()
#alignmentB.start()
#alignmentA.output()
#alignmentB.output()



def loadSeq(fileIn):
    refSeq = {}
    name = ''
    seq = ''
    for line in fileIn:
        if re.match(r'^[@|>]',line):
            if len(seq) > 2:
                refSeq[name] = seq
                seq = ""
            name = str(line)
        else:
            seq+=str(line.strip(' \t\n\r'))
    refSeq[name] = seq
    return refSeq

# TODO make a function that will just take the sequences of interest from the above map
# TODO so these sequences can be separated out and compared separately and then against each
# TODO other. To produce the full range of statistics.



# TODO get the main function to look like this
#diversity_analyses.main(arguments.shorah_input,arguments.output_file,arguments.reference_file,arguments.input_is_directory=False,arguments.contiguous_filter=10,arguments.identity=0.85,arguments.MEGA=None)

