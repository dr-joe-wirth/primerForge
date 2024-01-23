import sys

# class definition
class Clock:
    """ this class is used to easily track durations in a pretty format
    """
    import time as __time
    
    def __init__(self):
        """ accepts no inputs
            saves the start time and initializes member variables
        """
        self.__startTime:float = 0.0
        self.__duration:float = 0.0
        self.__start()
    
    def __start(self) -> None:
        """ accepts no inputs
            saves the start time
        """
        self.__startTime = Clock.__time.time()
    
    def __end(self) -> None:
        """ accepts no inputs
            saves the elapsed duration
        """
        self.__duration = Clock.__time.time() - self.__startTime
    
    def __parseDuration(self,digi:int) -> tuple[int,int,float]:
        """converts duration from seconds to hours, minutes, and seconds

        Args:
            digi (int): the number of decimal points for the seconds

        Returns:
            tuple[int,int,float]: hours, minutes, seconds
        """
        # constants
        TO_HRS = 3600
        TO_MIN = 60
        
        # parse the time; round seconds to specified digits
        hours = int(self.__duration // TO_HRS)
        minutes = int(self.__duration % TO_HRS // TO_MIN)
        seconds = round(self.__duration - minutes * TO_MIN - hours * TO_HRS, digi)
        
        return hours, minutes, seconds
    
    def __getDurationString(self, digi:int) -> str:
        """converts the duration to a formatted string

        Args:
            digi (int): the number of decimal points for the seconds
        
        Returns:
            str: the duration as a string in hh:mm:ss.ms format
        """
        # constant
        ZERO = "0"
        
        # parse the duration
        hours,minutes,seconds = self.__parseDuration(digi)
        
        # separate seconds from their decimals
        decimal = str(int(10**digi*seconds))
        seconds = str(int(seconds))
        
        # convert hours and minutes to strings
        hours = str(hours)
        minutes = str(minutes)
        
        # make sure hours are at least 2 chars long
        if len(hours) == 1:
            hours = ZERO + hours
        
        # make sure minutes are at least 2 chars long
        if len(minutes) == 1:
            minutes = ZERO + minutes
        
        # make sure seconds (int) are at least 2 chars long
        if len(seconds) == 1:
            seconds = ZERO + seconds
        
        # ensure that the decimal is exactly `digi` chars long
        if len(decimal) > digi:
            decimal = decimal[-digi:]
        elif len(decimal) < digi:
            decimal = ZERO * (digi - len(decimal)) + decimal
        
        # reconstruct the seconds to include decimals unless no decimals were requested
        if digi != 0:
            seconds = seconds + "." + decimal
        
        return ":".join([hours,minutes,seconds])
        
    def printTime(self, decimals:int=2) -> None:
        """prints the current time in hh:mm::ss format

        Args:
            decimals (int, optional): the number of decimal points for the seconds. Defaults to 2.
        """
        self.__end()
        print(self.__getDurationString(decimals))
    
    def getTime(self, decimals:int=2) -> str:
        self.__end()
        return self.__getDurationString(decimals)
    
    def restart(self) -> None:
        """ restarts the time on the clocker object
        """
        self.__start()

# functions
def _printStart(clock:Clock, msg:str, end:str=' ... ') -> None:
    """prints the start messages and restarts the clock

    Args:
        clock (Clock): a Clock object
        msg (str): the message to print
        end (str, optional): the end of the printed message. Defaults to ' ... '.
    """
    print(msg, end=end)
    sys.stdout.flush()
    clock.restart()


def _printDone(clock:Clock) -> None:
    """prints the end message and the duration

    Args:
        clock (Clock): a Clock object
    """
    print(f"done {clock.getTime()}" )
    sys.stdout.flush()
