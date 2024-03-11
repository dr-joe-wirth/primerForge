#!/usr/bin/env python3

__version__ = "0.5.3"
__author__ = "Joseph S. Wirth"

from bin.Parameters import Parameters
from bin.main import _main

if __name__ == "__main__":
    # parse command line arguments
    params = Parameters(__author__, __version__)
    
    if not params.helpRequested:
        try:
            # start up the logger
            params.log.rename(__name__)
            params.logRunDetails()
            
            # run primerForge
            _main(params)
        
        # catch all error messages
        except Exception as e:
            # save the error message if in debug mode
            params.log.rename(__name__)
            params.log.critical(e)

            # terminate        
            raise Exception(e)
