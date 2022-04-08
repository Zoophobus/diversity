#!/usr/bin/env python3.8

import sys,os,re
from numpy import mean
######
### UPDATES TO MAKE:
### 1) Update the SequenceIO class to address the possibility of 
###     file input
### 2) Integrate within the Django project to allow the statistics
###     to be presented in a figure
### 3) Correct the usage of a reference, enabling the module to work
###     without any provided reference. Presently the reference
###     indicates the overall length of the sequence from which the
###     records are obtained.
### 4) Improve the consensus sequence array to avoid the "last" empty
###     character.
######

class Base(object):
    def __init__(self,n,prvBase=None):
        self.data = n
        self.nxt = None
        self.prv = prvBase 
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

class Bases:
    def __init__(self,data = None):
        self.head = None
        self.tail = None
        if data != None and type(data) == str:
            bsItr = None
            for c in data:
                if bsItr == None:
                    self.head = Base(c)
                    bsItr = self.head
                else:
                    tmpItr = bsItr
                    newBase = Base(c,bsItr)
                    bsItr.set_next(newBase)
                    bsItr = bsItr.next()
            self.tail = bsItr
                
    def set(self,data):
        if type(data) == str:
            bsItr = None
            for c in data:
                if bsItr == None:
                    self.head = Base(c)
                    bsItr = self.head
                else:
                    tmpItr = bsItr
                    bsItr.set_next(Base(c))
                    bsItr.next()
                    bsItr.set_prev(tmpItr)
            self.tail = bsItr

    def next(self,bs=None):
        if bs == None:
            return self.head
        return bs.next()

    def prev(self,bs=None):
        if bs == None:
            return self.tail
        return bs.prev()

    def __str__(self):
        __str = ""
        itr = self.next(None)
        while itr != None:
            __str += itr.get_data()
            itr = self.next(itr)
        return __str


class Read(Bases):
    ''' 
        A class to define the reads that will are aligned to a 
        reference sequence and can be used to generate measures
        of variation in bases.

        Has a doubly linked list structure defined for the class 
        variable sequence
    '''
    def __init__(self,data=None,start = 0):
        self.begins = start
        self.finishes = self.begins + len(data)
        self.head = None
        if data is not None:
            super().__init__(data)

    def start(self, start = None):
        if start != None:
            self.begins = start
        else:
            return self.begins

    def end(self, end = None):
        if end != None:
            self.finishes = end
        else:
            return self.finishes

    def __str__(self):
        return self.sequence

class Sequence(object):
    """
    This class Sequence is a slight misnomer. 

    This class is used to define the sequence objects within
    a sliding window.

    After iterating for the length arguments the head of the
    "list" is to be updated, to maintain the length, until
    the end of the data.
    """
    def __init__(self,seq_name=None,data=None,length=66,start_offset=0):
        """
        """
        self.name = seq_name
        self.length = length
        self.start_offset = start_offset
        self.offset = start_offset
        self.data = Bases(data)
        self.head = self.data.next()
        self.tail = None
        self.base = None
        self.iteration = 0
        self.currentStr = ""
        self.title = ""
            
    def __str__(self):
        bIter = None
        currentStr = ">" + self.name + ":\t"
        for i in range(0,self.offset):
            currentStr += '-'
        while (bIter := self.next(bIter)) != None:
            currentStr += bIter.get_data()
        return currentStr

    def toString(self,base=-1):
        """
        Method:
            Produces a string for the current iteration of the
            window, with the sequence name included. THIS IS 
            PRIMARILY FOR USE WITH DEBUGGING!
        Method Input:
            this takes a parameter argument base. Upon first use of
            this method this should be set to None. To continue to
            iterate through the sequence, this can then be ignored.
        """

        # First we need to set up the initial string, including
        # the name and initialize some variables for iterating 
        if base == None:
            self.iteration = 0
            self.offset = self.start_offset
            self.head = self.data.head
            base = self.head
            seqHdChck = self.data.head
            self.currentStr = ""
            self.title = ">" + self.name + ":\t"
        else:
        # With the following while our "currentStr" is ten bases
        # long the actual string is longer, so we mask it 
            self.currentStr = re.subn(r'[a-zA-Z]','-',self.currentStr,1)[0]
            base = self.base
            seqHdChck = self.head

        while base != None:
            if seqHdChck != self.head:
                self.base = base
                return self.title + self.currentStr + "\n"
            if self.offset > 0:
                self.offset -= 1
                self.currentStr += '-'
            else:
                self.currentStr += base.get_data()
                base = self.next(base)
        self.base = base
        return self.title + self.currentStr  + "\n"

    def iterate(self,base=-1):
        """
        Method Input:
            this takes a parameter argument base. Upon first use of
            this method this should be set to None. To continue to
            iterate through the sequence, this can then be ignored.
        """

        # First we need to set up the initial string, including
        # the name and initialize some variables for iterating 
        if base == None:
            self.iteration = 0
            self.offset = self.start_offset
            self.head = self.data.head
            self.base = self.head
            return self.base
        self.base = self.next(self.base)
        return self.base

    def currentBase(self):
        """ 
        Method output:
            The handler that indicates where in iterating through
            the sequence "we" are.
        """
        return self.base

    def set(self,data,name="",length=66):
        """
        Method input:
            The sequence data, the name is optional, and length are
            optional.
        """
        self.data = data    
        self.name = name
        self.length = length
        self.head = self.data.next()
        self.tail= self.head

    def next(self,bs=None):
        """
        Method input:
            bs an object of type Base
        Method purpose:
            iterates through the reference to the Bases, or thereof
            derived Object. Contains bounds checks and limits on
            the tail, and head pointers for the reference.
        """
        self.iteration += 1
        if bs == None:
            return self.head
        if self.iteration >= self.length:
            self.control()
        return self.data.next(bs)

    def preview(self,bs=-1):
        if bs == None:
            return self.head
        return self.data.next(self.base)

    def lastBase(self,bs=None):
        """
        Method input:
            bs an object of type Base that is used to traverse the
            list. This methods does not have bounds checking and can
            overrun the intended list. NON-FATAL.
        """
        if bs == None:
            return self.head
        return self.data.prev(bs)

    def control(self):
        """
        Method output: this returns the first Base of the Sequence
            and resets the head variable to the second element on the
            Sequence
        """
        tmpObj = self.head
        self.head = self.data.next(tmpObj)
        return tmpObj

    def dropped(self):
        """
        Method output:
            Returns the base if there has been a character dropped,
            otherwise should return a None
        """
        return self.head.prev()

class SequenceObject(Sequence):
    def __init__(self,name,data,window_length,start_offset):
        super().__init__(name,data,window_length,start_offset)
        self.forward = None
        self.backward = None
    def nextSeq(self):
        return self.forward
    def prevSeq(self):
        return self.backward
    def setBack(self,sq):
        self.backward = sq
    def setForward(self,sq):
        self.forward = sq
    def __str__(self):
        return super().__str__()

class SequenceList(object):
    """
    CLASS PURPOSE:
    1) PROVIDE A LINKED LIST CONTAINER TO ITERATE
        THROUGH Sequence INSTANCES
    2) PROVIDE AN EASY FACILITY TO DROP SEQUENCES
        THAT HAVE ALREADY REACHED THE TERMINAL BASE
    3) ALLOW FOR NEW SEQUENCES TO BE APPENDED WHEN
        NEEDED
    """
    def __init__(self):
        self.first = None
        self.last = None

    def append(self,sq):
        """
        Method input:
            sq parameter, this parameter refers to an
            instance of the SequenceObject that is a wrapper
            to provide a linkedlist interface for a Sequence.
            This is performed through inheritance.
        """
        if self.first == None:
            self.first = sq
            self.last = sq
        else:
            self.first.setForward(sq)
            sq.setBack(self.first)
            self.first = sq

    def next(self,sqObj=None):
        """
        Method input:
            A handler to a SequenceObject instance, as this
            works via a linkedlist interface initial use requires
            no argument. Subsequent calls require the output to
            continue iterating.
        Method output:
            The corresponding handler to a SequenceObject instance.
            Used to gain the Sequence
        """

        if sqObj == None:
            return self.last
        return sqObj.nextSeq()
    
    def debugSequence(self,seqStrt = -1):
        """
        Method: To run through the Sequence instances printing
            them out with the windows in the order they are
            appended to the linked list. THIS IS PRIMARILY FOR
            USE IN DEBUGGING!
        Method input: seqStrt this is a value that indicates
            whether the call to debugSequence is the first such
            call. On the first call an argument with None type
            is required, subsequent calls require no input value.
        """
        seqObject = None
        algnStrngs = ""
        while (seqObject := self.next(seqObject)) != None:
            algnStrngs += seqObject.toString(seqStrt)

        return algnStrngs

    def drop(self):
        """
        Method output:
            This returns an integer indicating the number of sequences
            that have been iterated through. These should be subsequently
            dropped. With a sorted list these are the first elements in
            the linked list.
        """
        sq = None
        lost = 0
        old = None
        while (sq := self.next(sq)) != None:
            if sq.preview() == None: 
                # to help with debugging and for clarity we save a
                # sequence we will remove
                old = sq
                if sq.prevSeq() != None and sq.nextSeq() != None:
                    # if the sequence is in the middle of the linkedlist
                    sq.prevSeq().setForward(sq.nextSeq())
                    sq.nextSeq().setBack(sq.prevSeq())
                elif sq.prevSeq() != None:
                    # if the sequence is at the end of the list
                    sq.prevSeq().setForward(None)
                    self.first = sq.prevSeq()
                elif sq.nextSeq() != None:
                    # if the sequence is at the beginning of the list
                    sq.nextSeq().setBack(None)
                    self.last = sq.nextSeq()
                lost += 1
        if old != None:
            return old.name
        return None



class Windows(object):
    """
    CLASS PURPOSE:
    1) TO PROVIDE THE VARIABLES AND CODE FOR 
        CALCULATING SEQUENCE STATISTICS
    2) PROVIDE THE MECHANISM FOR MAKING A SINGLE
        ITERATION THROUGH THE CONTAINED SEQUENCES
    """
    def __init__(self,window_length,reference_length):
        self.windows_start_at = 0
        self.windows_end_at = reference_length
        self.sequences = SequenceList()
        self.diversity = 0
        self.seq_diversity = []
        self.consensus = []
        self.location = 0
        self.thisMismatch = 0 
        self.window_length = window_length
    def size(self):
        return self.no_windows
    def insert(self,seq):
        """
        Method input: seq is of type Bases, or derived thereof
        Method purpose: as an interface for adding sequence/read
            references to the Window object
        """
        self.sequences.append(SequenceObject(seq.name,seq.sequence,self.window_length,seq.start_offset))
    def theta(self,version):
        self.wattersons = version

    def valid_data(self,aChar):
        if aChar == '-' or aChar == "?" or aChar == "N" or aChar == "n":
            return False 
        return True

    def restart(self):
        self.sequences = SequenceList()
        self.diversity = 0
        self.seq_diversity = []
        self.consensus = []
        self.location = 0
        self.thisMismatch = 0 


    def setWindow(self,seqs,other=None,end=1):
        """
        Method input: 
            seqs an array that should contain references to Objects 
                derived from the Bases class

            ref a reference to a Bases object that contains the reference

            other a secondary window object, for the moment this is not
                implemented
        METHOD DEFUNCT
        """
        self.windows_start_at = 0
        self.windows_end_at = end

        for seq in seqs:
            self.insert(seq)

    def slide(self,startSequences=-1,sortedReads=None):
        """ 
        Method input:
            sortedReads a reference to the collection of reads, this should
                be ordered by starting position, so that these can be added 
                in a simple iterative style.
            startSequences a parameter that should be set to None for the
                first occurrence of the methods use, triggers the 
                initialization for the sliding window

        Method output:
            an integer indicating the current position that the window is 
                at within the sequences, this position should be relative
                to the reference
        """
        """ ENFORCING A SLIDE OF ONLY A SINGLE BASE """
        def updateMap(m,x):
            if not x.get_data() in m and self.valid_data(x.get_data()):
                m[x.get_data()] = 1
            else:
                m[x.get_data()] += 1
        def mismatches(m,s):
            mm = 0
            k = len(m)
            for frq in m.values():
                k -= 1
                if k > 0:
                    mm += frq * (s - frq)
            if mm == 0:
                return 0
            return mm / (s *(s-1)/2)

        mismatch = 0
        newBase = {}
        lostBase = {}
        latestSum = 0
        headBase = None
        headSum = 0

        seqObject = None
        while (seqObject := self.sequences.next(seqObject)) != None:
            # for each Sequence instance available iterate ONCE!
            latestBase = (seqObject.iterate(seqObject.currentBase()))

            # Check the newest base in the window
            if latestBase != None:
                latestSum += 1 
                updateMap(newBase,latestBase)
            # Find the bases that are dropped, or None type
            headBase = seqObject.dropped() 
            # Check the dropped bases
            if headBase != None:
                headSum += 1
                updateMap(lostBase,headBase)
        # adjust the number of sequences, if we need to drop some
        finished_seq = self.sequences.drop()
        # calculate the change in the numbers of mismatches per comparison
        mismatch = mismatches(newBase,latestSum) - mismatches(lostBase,headSum)
        # record the consensus
        self.consensus.append(self.get_consensus(newBase)) 
        # adjust the running score for the number of mismatches per
        # comparison
        self.diversity +=  mismatch
        divisor = self.window_length

        # then find the window length (towards the beginning and end
        # of the reference this won't equal the value set) For the 
        # per base value
        if self.location >= self.window_length and self.location <= (self.windows_end_at - self.window_length/2):
            # numerically meaningless calculations to be avoided
            if self.diversity == 0:
                self.seq_diversity.append(0)
            else:
                self.seq_diversity.append(self.diversity/divisor)
        self.location += 1
        return self.location


    def get_consensus_string(self): 
        return self.consensus
    def get_consensus(self,consensusMap):
        """
        Method input: a map that contains the bases that occur at a single
            site, the base that has the highest frequency will be returned
        Method output: the base from the input that has the highest 
            frequency
        """
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
        return mean(self.seq_diversity)

    def __str__(self):
        theStr = self.consensus()
        theStr += "\n"
        return theStr

class SequenceIO(object):
    """
    CLASS PURPOSE:
    1) PROVIDE A COMMON INTERFACE FOR SEQUENCE INPUT TO BE
        READ AND UTILIZED
        (CURRENTLY ONLY DEVELOPED TO USE A LIST AS INPUT)
    2) ALLOW THE SORTING OF UNSORTED DATA (REQUIRES PARTICULAR
        INPUT STYLE)
    """
    def __init__(self,inputType,inputRef):
        self.type = inputType
        self.input = inputRef
        self.position = 0

    def sort(self):
        if self.type == 'file':
            return self.input
        elif self.type == 'list':
            self.input = sorted(self.input,key = lambda x: x.start_offset) 
        else:
            raise Exception("Poor input type for sorting, only 'file' or 'list' may be specified")

    def next(self):
        self.position += 1
        return self.input[self.position-1]
    def where(self):
        if self.position >= len(self.input):
            return -1
        return self.input[self.position].start_offset

class Alignment(object):
    """
    CLASS PURPOSE:
    1) THROUGH INCLUDING THE SOURCE OF INPUT AND THE
        WINDOW ALLOW FOR THE DATA TO BE USED MORE 
        EFFICIENTLY, LOADING SEQUENCE DATA WHEN THEY
        NEED TO BE LOADED
    2) PROVIDE THE FACILITY FOR EXTRACTING STATISTICS
    """
    def __init__(self,sequence_input,window_length=66,reference=None):
        """
        Method input:
            sequence_input a necessary dictionary object of the format
                { 'file' : <file-name>, 'list' <list-reference> }
                either, or one of these needs to be included. If both are
                included errors should be expected. The list should be a
                reference to a list of educational.Sequence type objects.
        """
        if 'file' in sequence_input and 'list' in sequence_input:
            raise Exception("Too many parameters passed to Alignment.__init__")
        elif 'file' in sequence_input:
            self.input = SequenceIO('file',sequence_input['file'])
        else:
            self.input = SequenceIO('list',sequence_input['list'])
        self.sort()
        self.sequences = []
        self.reference = reference
        if self.reference != None:
            self.reference_length = len(self.reference)
        else:
            self.reference_length = 0
        self.windows = Windows(window_length,self.reference_length)

    def addReference(self,reference):
        self.reference = reference
        self.reference_length = len(self.reference)

    def sort(self):
        return self.input.sort()

    def run(self,other=None):
        position = 0
        while position <= self.reference_length:
            while self.input.where() == position:
                self.windows.insert(self.input.next())
            self.windows.slide()
            position += 1
 
    def output(self,fileOutstream,other=None):
        pass

    def __str__(self):
        theSeq = self.windows.get_consensus() + "\n" + self.windows.get_diversity()
        return theSeq
    def get_consensus(self):
        return self.windows.get_consensus_string()
    def get_diversity(self):
        return self.windows.get_diversity()




class SequenceInput(object):
    def __init__(self,name=None,data=None,offset=0):
        self.name = name
        self.sequence = data
        self.start_offset = offset



if __name__ == '__main__':
    seq1 = "CGTACAGTTACAAAGTTTGTTACTTGTCGTTACGGGCCGATCG"
    seq2 = "GTACAGTTACAAAGTATGTTACTTGTCGTTACGGG"
    seq3 = "GGTACAGATACAAAGTTTGTAACTTGTCGTTACGCGCCG"
    seq4 = "TTGGTACAGATACAAAGTTTGTAACTTGTCGTTACGCGCTGATC"
    ref  = "TTCGAACAGTTACAAGGTTTGTTACTTGTGGTTACGGGCCGATCG"
    examples = [SequenceInput('seq1',seq1,2),SequenceInput('seq2',seq2,3),SequenceInput('seq3',seq3,2),SequenceInput('seq4',seq4,0)]

 
    testing = Alignment({'list':examples},10,ref)
    testing.run()

    print(testing.get_consensus())
    print(testing.get_diversity())
    

