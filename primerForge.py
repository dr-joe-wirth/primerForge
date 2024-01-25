#!/usr/bin/env python3
from bin.main import _main, Parameters, __version__, __author__

if __name__ == "__main__":
    # parse command line arguments
    params = Parameters()
    
    # save the run details if debugging
    if params.debug:
        params.log.setLogger(__name__)
        params.saveRunDetails()
    
    # catch all error messages
    try:
        _main(params)
        
    except Exception as e:
        # save the error message if in debug mode
        if params.debug:
            params.log.setLogger(__name__)
            params.log.writeCriticalMsg(e)

        # terminate        
        raise Exception(e)
