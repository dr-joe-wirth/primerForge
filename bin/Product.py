from __future__ import annotations
from typing import Union

class Product:
    """class for storing PCR product data for a pair of primers
    """
    def __init__(self, contig:str, size:Union[int,str], fbin:int, rbin:int, dimerTm:float, fStart:int, fEnd:int, fStrand:str, rStart:int, rEnd:int, rStrand:str) -> Product:
        """constructor

        Args:
            contig (str): the contig where the product resides
            size (Union[int,str]): the contig size (will be str if multiple sizes)
            fbin (int): the bin the forward primer belongs to
            rbin (int): the bin the reverse primer belongs to
            dimerTm (float): the heterodimer Tm for the primer pair
            fStart (int): the start of the forward primer
            fEnd (int): the end on the forward primer
            fStrand (str): the strand of the forward primer
            rStart (int): the start of the reverse primer
            rEnd (int): the end on the reverse primer
            rStrand (str): the strand of the reverse primer            

        Returns:
            Product: a Product object
        """
        # set attributes
        self.contig:str = contig
        self.size:Union[int,str] = size
        self.dimerTm:float = dimerTm
        self.fbin:int = fbin
        self.rbin:int = rbin
        self.fStart:int = fStart
        self.fEnd:int = fEnd
        self.fStrand:str = fStrand
        self.rStart:int = rStart
        self.rEnd:int = rEnd
        self.rStrand:str = rStrand
    
    # overloads
    def __str__(self) -> str:
        return f'{self.contig}: {self.size} bp'
    
    def __repr__(self) -> str:
        return str(self)
