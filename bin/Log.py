from __future__ import annotations
import logging, os

class Log():
    """class for interacting with the logger 
    """
    # global constants
    __LOG_FN = "primerForge.log"
    
    def __init__(self, debugDir:str='', debug:bool=False, initialize=True) -> Log:
        """creates a Log object

        Args:
            debugDir (str, optional): the directory for debugging. Defaults to ''.
            debug (bool, optional): indicates if logging is debug level. Defaults to False.
            initialize (bool, optional): indicates if the logger should be initialize. Defaults to True

        Returns:
            Log: a Log object
        """
        # type hint attributes
        self.__logger:logging.Logger
        
        # set the debug directory
        if debugDir == '':
            self.debugDir:str = os.getcwd()
        else:
            self.debugDir = debugDir
        
        # make sure debug directory exists
        if not os.path.isdir(self.debugDir):
            os.mkdir(self.debugDir)
        
        # set the log file
        self.logFn:str = os.path.join(self.debugDir, Log.__LOG_FN)
        
        # point the log at the log file; level is DEBUG or INFO
        if debug and initialize:
            logging.basicConfig(filename=self.logFn, level=logging.DEBUG)
        elif initialize:
            logging.basicConfig(filename=self.logFn, level=logging.INFO)
    
    def rename(self, name:str):
        """renames the logger

        Args:
            name (str): the new name of the logger
        """
        self.__logger = logging.getLogger(name)
    
    def info(self, msg:str) -> None:
        """writes message as logger.info

        Args:
            msg (str): the message to write
        """
        self.__logger.info(msg)
    
    def debug(self, msg:str) -> None:
        """writes message as logger.debug

        Args:
            msg (str): the message to write
        """
        self.__logger.debug(msg)
    
    def error(self, msg:str) -> None:
        """writes message as logger.error

        Args:
            msg (str): the message to write
        """
        self.__logger.error(msg)
    
    def critical(self, msg:str) -> None:
        """writes message as logger.critical

        Args:
            msg (str): the message to write
        """
        self.__logger.critical(msg)
