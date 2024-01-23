from __future__ import annotations
import logging, os

class Log():
    """class to interact with the logger 
    """
    # global constants
    __DEBUG_DIR = "_debug"
    __LOG_FN = "primerForge.log"
    
    def __init__(self) -> Log:
        # type hint attributes
        self.__logger:logging.Logger
        self.debugDir:str = os.path.join(os.getcwd(), Log.__DEBUG_DIR)
        
        # make sure debug directory exists
        if not os.path.isdir(Log.__DEBUG_DIR):
            os.mkdir(Log.__DEBUG_DIR)
        
        # initialize the logging file
        logging.basicConfig(filename=os.path.join(self.debugDir, Log.__LOG_FN), level=logging.DEBUG)
    
    def setLogger(self, name:str) -> None:
        """sets the name of the logger

        Args:
            name (str): the name of the logger
        """
        self.__logger = logging.getLogger(name)
    
    def writeInfoMsg(self, msg:str) -> None:
        """writes message as logger.info

        Args:
            msg (str): the message to write
        """
        self.__logger.info(msg)
    
    def writeDebugMsg(self, msg:str) -> None:
        """writes message as logger.debug

        Args:
            msg (str): the message to write
        """
        self.__logger.debug(msg)
    
    def writeErrorMsg(self, msg:str) -> None:
        """writes message as logger.error

        Args:
            msg (str): the message to write
        """
        self.__logger.error(msg)
    
    def writeCriticalMsg(self, msg:str) -> None:
        """writes message as logger.critical

        Args:
            msg (str): the message to write
        """
        self.__logger.critical(msg)
