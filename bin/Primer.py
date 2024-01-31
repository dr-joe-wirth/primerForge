from __future__ import annotations
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp

class Primer:
    """a class for calculating and storing primer data
    """
    def __init__(self, seq:Seq, contig:str, start:int, length:int) -> Primer:
        """creates a Primer object

        Args:
            seq (Seq): primer sequence
            contig (str): the contig name
            start (int): the start position on the contig
            length (int): the primer length

        Returns:
            Primer: the primer
        """
        # initialize the attributes
        self.seq:Seq = None
        self.start:int = start
        self.end:int = start + length - 1
        self.contig:str = contig
        self.Tm:float = None
        self.gcPer:float = None
        self.__length:int = length
        
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
        return self.seq == other.seq
    
    def __ne__(self, other:Primer) -> bool:
        return not self.seq == other.seq
    
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
        self.Tm = MeltingTemp.Tm_Wallace(self.seq)
    
    # public methods
    def reverseComplement(self) -> Primer:
        """reverse complements the calling object

        Returns:
            Primer: the reverse complement of the calling object
        """
        new = Primer(self.seq.reverse_complement(), self.contig, self.end, len(self))
        new.end = self.start
        
        return new

    def getMinimizer(self, lmerSize:int) -> Seq:
        """finds the minimizer for the calling object

        Args:
            lmerSize (int): the length of the l-mer

        Raises:
            ValueError: primer length must exceed l-mer length

        Returns:
            Seq: the minimizer sequence
        """
        # constants
        ERR_MSG = "Window size should be less than or equal to the k-mer length."
        
        # make sure the lmer is smaller than the primer
        if len(self.seq) < lmerSize:
            raise ValueError(ERR_MSG)

        # Initialize the minimizer and its position
        minimizer = self.seq[:lmerSize]

        # Iterate through the k-mer with the sliding window
        for i in range(1, len(self.seq) - lmerSize + 1):
            current_window = self.seq[i:i + lmerSize]

            # Update the minimizer if the current window is lexicographically smaller
            if current_window < minimizer:
                minimizer = current_window

        return minimizer
