from __future__ import annotations
import primer3
from Bio.Seq import Seq

class Primer:
    """a class for calculating and storing primer data
    """
    PLUS = "+"
    MINUS = "-"
    
    def __init__(self, seq:Seq, contig:str, start:int, length:int, strand:str) -> Primer:
        """creates a Primer object

        Args:
            seq (Seq): primer sequence
            contig (str): the contig name
            start (int): the start position on the contig's (+) strand
            length (int): the primer length
            strand (str): the strand of the contig; "+" or "-"

        Returns:
            Primer: the primer
        """
        # make sure valid strand was passed
        if strand not in (Primer.PLUS, Primer.MINUS):
            raise ValueError("invalid strand specified")
        
        # initialize the attributes
        self.seq:Seq = None
        self.start:int = start
        self.end:int = start + length - 1
        self.contig:str = contig
        self.strand:str = strand
        self.Tm:float = None
        self.gcPer:float = None
        self.__length:int = length
        
        # flip start and end if on the minus strand
        if self.strand == Primer.MINUS:
            self.start = self.end
            self.end = start
        
        # run import methods
        self.__importSeq(seq)
        self.__calcPerGc()
        self.__calculateTm()
    
    # overloads
    def __hash__(self) -> int:
        return self.seq.__hash__()
    
    def __len__(self) -> int:
        return self.__length
    
    def __str__(self) -> str:
        return self.seq.__str__()
    
    def __repr__(self) -> str:
        return str(self)
    
    def __eq__(self, other:Primer) -> bool:
        # compare two primers
        if type(other) is Primer:
            return self.seq == other.seq
        
        # handle comparisons to Seq and str
        elif type(other) in (Seq, str):
            return self.seq == other
        
        # handle comparisons to invalid objects
        else:
            raise ValueError(f"cannot compare Primer to type '{type(other)}'")
    
    def __ne__(self, other:Primer) -> bool:
        return not self.seq == other.seq
    
    def __format__(self, format_spec, /):
        return str(self.seq).__format__(format_spec)
    
    # private methods
    def __importSeq(self, seq:Seq) -> None:
        # make sure that the sequence is upper-case
        self.seq = Seq(seq.upper())
    
    def __calcPerGc(self) -> None:
        # constant
        GC = ("G", "C")
        
        # count number of GCs
        numGc = 0
        for base in self.seq:
            if base in GC:
                numGc += 1
        
        # calculate and store the percentage of GC
        self.gcPer = numGc / len(self) * 100

    def __calculateTm(self) -> None:
        # calculate and store the melting temp
        self.Tm = primer3.calc_tm(str(self.seq))
    
    # public methods
    def reverseComplement(self) -> Primer:
        """reverse complements the calling object

        Returns:
            Primer: the reverse complement of the calling object
        """
        # figure out the reverse strand
        if self.strand == Primer.PLUS:
            new = Primer(self.seq.reverse_complement(), self.contig, self.start, len(self), Primer.MINUS)
        else:
            new = Primer(self.seq.reverse_complement(), self.contig, self.end, len(self), Primer.PLUS)
            
        # make the new object
        return new

    def getMinimizer(self, lmerSize:int, strand:str) -> Seq:
        """finds the minimizer for the calling object

        Args:
            lmerSize (int): the length of the l-mer
            strand (str): the strand to evaluate; "+" or "-"

        Raises:
            ValueError: primer length must exceed l-mer length

        Returns:
            Seq: the minimizer sequence
        """
        # constants
        ERR_MSG_1 = "Window size should be less than or equal to the k-mer length."
        ERR_MSG_2 = "Invalid strand"
        
        # make sure the lmer is smaller than the primer
        if len(self.seq) < lmerSize:
            raise ValueError(ERR_MSG_1)
        
        # make sure the strand is valid
        if strand not in (Primer.PLUS, Primer.MINUS):
            raise ValueError(ERR_MSG_2)

        # get the appropriate sequence
        if self.strand == strand:
            seq = self.seq
        else:
            seq = self.seq.reverse_complement()

        # Initialize the minimizer and its position
        minimizer = seq[:lmerSize]

        # Iterate through the k-mer with the sliding window
        for i in range(1, len(seq) - lmerSize + 1):
            currentWindow = seq[i:i + lmerSize]

            # Update the minimizer if the current window is lexicographically smaller
            if currentWindow < minimizer:
                minimizer = currentWindow

        return minimizer
