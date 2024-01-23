from __future__ import annotations
import logging, os

class Debugger():
    __DEBUG_DIR = "_debug"
    __LOG_FN = "primerForge.log"
    __DONE = "done"
    
    def __init__(self) -> Debugger:
        self.__logger:logging.Logger
        self.debugDir:str = os.path.join(os.getcwd(), Debugger.__DEBUG_DIR)
        
        if not os.path.isdir(Debugger.__DEBUG_DIR):
            os.mkdir(Debugger.__DEBUG_DIR)
        
        logging.basicConfig(filename=os.path.join(self.debugDir, Debugger.__LOG_FN), level=logging.DEBUG)
    
    def setLogger(self, name:str) -> None:
        self.__logger = logging.getLogger(name)
    
    def writeDebugMsg(self, msg:str) -> None:
        self.__logger.debug(msg)
    
    def writeErrorMsg(self, msg:str) -> None:
        self.__logger.error(msg)