from __future__ import annotations
from bin.Primer import Primer


class Level():
    """ class for storing the dataset level
    """
    # private attributes
    __NUM_TO_STR = {0: 'candidate kmer',
                    1: 'unfiltered pair',
                    2: 'filtered pair',
                    3: 'final pair'}
    
    __STR_TO_NUM = {val:key for key,val in __NUM_TO_STR.items()}

    # overloads
    def __init__(self) -> Level:
        """constructor; takes no arguments

        Returns:
            Level: a Level object
        """
        # start levels at the lowest position
        self.__num = 0
      
    def __str__(self) -> str:
        """provides a string representation of the calling object

        Returns:
            str: the string representation
        """
        return Level.__NUM_TO_STR[self.__num]

    def __repr__(self) -> str:
        """provides a string representation of the calling object

        Returns:
            str: the string representation
        """
        return str(self)
    
    def __hash__(self) -> int:
        """hashes the calling object

        Returns:
            int: the hashed object
        """
        return str(self).__hash__()
    
    def __iadd__(self, other:int) -> Level:
        """increments the calling object

        Args:
            other (int): the integer to add

        Returns:
            Level: the modified calling object
        """
        self.__num += other
        return self
    
    def __isub__(self, other:int) -> Level:
        """decrements the calling object

        Args:
            other (int): the integer to subtract

        Returns:
            Level: the modified calling object
        """
        # use += operater on negative input
        self += -other
        return self
    
    def __add__(self, other:int) -> Level:
        """creates a new Level object by adding to it

        Args:
            other (int): the integer to add

        Raises:
            Exception: cannot create a valid object with the input

        Returns:
            Level: a new Level object
        """
        # messages
        ERR_MSG_A = "cannot add "
        ERR_MSG_B = " to "
        
        # create a new Level object at the specified position
        out = Level()
        out.__num = self.__num + other
        
        if out.__num > max(Level.__NUM_TO_STR.keys()) or out.__num < min(Level.__NUM_TO_STR.keys()):
            raise Exception(f"{ERR_MSG_A}{other}{ERR_MSG_B}{self}")
        
        return out
    
    def __sub__(self, other:int) -> Level:
        """creates a new Level object by subtracting from it

        Args:
            other (int): the integer to add

        Returns:
            Level: a new Level object
        """
        # use the + operator on the negative input
        return self + -other

    def __eq__(self, other) -> bool:
        """== operator

        Args:
            other (_type_): can be str or Level

        Raises:
            ValueError: other is not str or Level

        Returns:
            bool: indicates if the two objects are equivalent
        """
        # handle strings
        if type(other) is str:
            otherNum = Level.__STR_TO_NUM[other]
        
        # handle Levels
        elif type(other) is Level:
            otherNum = other.__num
        
        # raise if neither string nor Level
        else:
            raise ValueError(f"cannot compare '{type(other)}' objects to {Level.__name__}")
        
        # equivalence based on the number
        return self.__num == otherNum
    
    def __gt__(self, other) -> bool:
        """> operator

        Args:
            other (_type_): string or Level

        Raises:
            ValueError: not string or Level

        Returns:
            bool: self > other
        """
        # handle strings
        if type(other) is str:
            otherNum = Level.__STR_TO_NUM[other]
        
        # handle Levels
        elif type(other) is Level:
            otherNum = other.__num

        # raise error if neither string nor Level        
        else:
            raise ValueError(f"cannot compare '{type(other)}' objects to {Level.__name__}")
        
        return self.__num > otherNum
    
    def __lt__(self, other) -> bool:
        """< operator

        Args:
            other (_type_): string or Level

        Raises:
            ValueError: not string or Level

        Returns:
            bool: self < other
        """
        # handle strings
        if type(other) is str:
            otherNum = Level.__STR_TO_NUM[other]
        
        # handle Levels
        elif type(other) is Level:
            otherNum = other.__num
        
        # raise if not string nor Level
        else:
            raise ValueError(f"cannot compare '{type(other)}' objects to {Level.__name__}")
        
        return self.__num < otherNum
    
    def __ne__(self, other) -> bool:
        """!= operator

        Args:
            other (_type_): can be str or Level

        Returns:
            bool: indicates if the two objects are inequivalent
        """
        return not self == other
    
    def __ge__(self, other) -> bool:
        """>= operator

        Args:
            other (_type_): string or Level

        Returns:
            bool: self >= other
        """
        return self == other or self > other
    
    def __le__(self, other) -> bool:
        """<= operator

        Args:
            other (_type_): string or Level

        Returns:
            bool: self <= other
        """
        return self == other or self < other


class AnalysisData():
    """class for storing data for downstream analysis of kmer distribution
    """
    # public attributes
    LEVELS = (Level(),
              Level() + 1,
              Level() + 2,
              Level() + 3)
    
    # overloads
    def __init__(self, primer:Primer, index:int, name:str) -> AnalysisData:
        """constructor

        Args:
            primer (Primer): a Primer object
            index (int): the index of this object
            name (str): the name of the genome

        Returns:
            AnalysisData: an AnalysisData object
        """
        # public attributes
        self.primer:Primer = primer
        self.name:str = name
        
        # set private attributes
        self.__index:int = index
        self.__pair:set[int] = set()
        self.__level:Level = Level()
    
    def __str__(self) -> str:
        """provides a string repesentation of the calling object

        Returns:
            str: the string representation
        """
        return f"{str(self.__index) + ':':<8}{self.name} ({self.__level})"
    
    def __repr__(self) -> str:
        """provides a string repesentation of the calling object

        Returns:
            str: the string representation
        """
        return str(self)
    
    def __hash__(self) -> int:
        """hashes calling object

        Returns:
            int: the hash
        """
        return f"{str(self.primer)}{self.primer.contig}".__hash__()
    
    # public methods
    def updatePairs(self, pairNum:int) -> None:
        """adds another index to the set of pairs

        Args:
            pairNum (int): the other index
        """
        self.__pair.add(pairNum)
    
    def getPairs(self) -> list[int]:
        """returns a list of the pairs

        Returns:
            list[int]: a list of the pairs
        """
        return list(self.__pair)
    
    def incrementLevel(self) -> None:
        ERR_MSG = "cannot increment data set"
        
        if self.__level == max(AnalysisData.LEVELS):
            raise Exception(ERR_MSG)
        
        self.__level += 1
    
    def getLevel(self) -> str:
        return self.__level

    def getIndex(self) -> int:
        return self.__index
