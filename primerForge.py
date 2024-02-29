#!/usr/bin/env python3

__version__ = "0.5.0"
__author__ = "Joseph S. Wirth"

from bin.Parameters import Parameters
from bin.main import _main

if __name__ == "__main__":
    # parse command line arguments
    params = Parameters(__author__, __version__)
    
    # save the run details if debugging
    if params.debug:
        params.log.initialize(__name__)
        params.logRunDetails()
    
    # catch all error messages
    try:
        _main(params)
        
    except Exception as e:
        # save the error message if in debug mode
        if params.debug:
            params.log.rename(__name__)
            params.log.critical(e)

        # terminate        
        raise Exception(e)
