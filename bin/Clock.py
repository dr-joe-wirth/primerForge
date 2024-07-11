from __future__ import annotations
import os, sys, time
from multiprocessing import Process, Event


class Wheel:
    """provides a wheel spinning
    """
    # public attribute
    PAUSE = 0.12
    
    # priave attribute
    __CHARS = "-\|/"
    __EVENT = Event()
        
    # overloads
    def __init__(self) -> Wheel:
        """constructor

        Returns:
            Wheel: a Wheel object
        """
        # initialize attributes
        self.__msg:str = ''
        self.__process:Process

    # private methods
    def __spin(self) -> None:
        """spinning function
        """
        # keep spinning until it is time to stop
        while not self.__EVENT.is_set():
            # print each character in the wheel then pause
            for char in Wheel.__CHARS:
                sys.stdout.write('\r' + self.__msg + char)
                sys.stdout.flush()
                time.sleep(Wheel.PAUSE)
    
    # public methods
    def start(self, msg:str) -> None:
        """starts the wheel spinning

        Args:
            msg (str): the message to preceed the wheel
        """
        # initialize attributes
        self.__msg = msg
        self.__process = Process(target=self.__spin)
        self.__EVENT.clear()
        
        # start spinning the wheel
        self.__process.start()
    
    def stop(self) -> None:
        """stops spinning the wheel
        """
        # this will raise AttributeError if the wheel isn't spinning
        try:
            # stop spinning the wheel
            if self.__process.is_alive():
                self.__EVENT.set()
                self.__process.join()
                
                # remove the wheel character
                sys.stdout.write('\r' + self.__msg)
                sys.stdout.flush()
        
        except AttributeError:
            pass


class Clock:
    """ this class is used to easily track durations in a pretty format
    """
    # private attributes; share a single wheel (easier to kill)
    __WHEEL = Wheel()
    
    # overload
    def __init__(self) -> Clock:
        """ constructor. accepts no inputs
            saves the start time and initializes member variables
        """
        # initialize attributes
        self.__startTime:float = 0.0
        self.__duration:float = 0.0
        self.__spin:bool = True
        
        # start the clock
        self.__start()
    
    # private methods
    def __start(self) -> None:
        """ accepts no inputs
            saves the start time
        """
        self.__startTime = time.time()
    
    def __end(self) -> None:
        """ accepts no inputs
            saves the elapsed duration
        """
        self.__duration = time.time() - self.__startTime
    
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

    # public methods
    def getTime(self, decimals:int=2) -> tuple[int,int,float]:
        """gets the elapsed time as a tuple

        Args:
            decimals (int, optional): the number of deimal points for the seconds. Defaults to 2.

        Returns:
            tuple[int,int,float]: hours, minutes, seconds
        """
        self.__end()
        return self.__parseDuration(decimals)

    def getTimeString(self, decimals:int=2) -> str:
        """gets the elapsed time as a string in hh:mm:ss.ms format

        Args:
            decimals (int, optional): the number of decimal points for the seconds. Defaults to 2.

        Returns:
            str: hh:mm:ss.ms
        """
        self.__end()
        return self.__getDurationString(decimals)
      
    def printTime(self, decimals:int=2) -> None:
        """prints the current time in hh:mm::ss.ms format

        Args:
            decimals (int, optional): the number of decimal points for the seconds. Defaults to 2.
        """
        self.__end()
        print(self.__getDurationString(decimals))

    def restart(self) -> None:
        """ restarts the clock object
        """
        self.__start()
    
    def printStart(self, msg:str, prefix:str='', end:str=' ... ', spin:bool=True,) -> None:
        """prints the start message and restarts the clock

        Args:
            msg (str): the message to print
            prefix (str, optional): the beginning of the printed message. Defaults to ''.
            end (str, optional): the end of the printed message. Defaults to ' ... '.
            spin (bool, optional): indicates if the wheel should spin. Defaults to True.
        """
        # save the spin status
        if not "CI" in os.environ:
            self.__spin = spin
        
        # do not spin if we are doing CI
        else:
            self.__spin = False
        
        # spin the wheel if requested; wheel handles printing
        if self.__spin:
            self.__WHEEL.start(prefix + msg + end)
        
        # otherwise print the message 
        else:
            print(prefix + msg, end=end)
            sys.stdout.flush()
        
        # restart the clock
        self.restart()
    
    def printDone(self) -> None:
        """prints the end message and the duration
        """
        # stop spinning the wheel if necessary
        Clock._killWheel()
        
        # print the done string
        print(f"done {self.getTimeString()}")
        sys.stdout.flush()
        
        Clock.__WHEEL.__msg = ''
    
    def _killWheel() -> None:
        """kills any spinning clocks
        """
        Clock.__WHEEL.stop()
