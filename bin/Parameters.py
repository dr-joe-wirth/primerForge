from __future__ import annotations
import getopt, glob, os, pickle, sys
from bin.Debugger import Debugger
from bin.main import __version__

class Parameters():
    def __init__(self) -> Parameters:
        # type hint attributes
        self.ingroupFns:list[str]
        self.outgroupFns:list[str]
        self.outFn:str
        self.format:str
        self.minLen:int
        self.maxLen:int
        self.minGc:float
        self.maxGc:float
        self.minTm:float
        self.maxTm:float
        self.minPcr:int
        self.maxPcr:int
        self.maxTmDiff:float
        self.disallowedLens:range
        self.numThreads:int
        self.debug:bool
        self.helpRequested:bool
        self.debugger:Debugger
        
        # parse the command line arguments
        self.__parseArgs()
        
        # initialize a log if debugging
        if self.debug:
            self.__debugSetup()
    
    def __parseArgs(self) -> None:
        """parses command line arguments

        Raises:
            ValueError: invalid ingroup file(s)
            ValueError: invalid outgroup file(s)
            ValueError: two bad pcr lengths no specified
            ValueError: bad pcr lengths are not integers
            ValueError: invalid file format
            ValueError: invalid primer length
            ValueError: primer lengths not integers
            ValueError: specify range of GC
            ValueError: GC values non-numeric
            ValueError: specify range of Tm
            ValueError: Tm values non-numeric
            ValueError: invalid PCR product length
            ValueError: PCR product lengths not integers
            ValueError: Tm difference is non-numeric
            ValueError: num threads not an integer
            ValueError: must specify an ingroup file
            ValueError: must specify an output file

        Returns:
            tuple[list[str],list[str],str,str,int,int,float,float,float,float,int,int,float,range,int,bool,bool]:
                ingroupFns,outgroupFns,outFN,format,minPrimerLen,maxPrimerLen,minGc,maxGc,minTm,maxTm,minPcrLen,maxPcrLen,maxTmDiff,disallowedLens,numThreads,helpRequested,debug
        """
        # constants
        ALLOWED_FORMATS = ('genbank', 'fasta')
        SEP = ","
        
        # flags
        INGROUP_FLAGS = ('-i', '--ingroup')
        OUT_FLAGS = ('-o', '--out')
        OUTGROUP_FLAGS = ('-u', '--outgroup')
        DISALLOW_FLAGS = ("-b", "--bad_outgroup_pcr_len")
        FMT_FLAGS = ('-f', '--format')
        PRIMER_LEN_FLAGS = ('-p', '--primer_len')
        GC_FLAGS = ('-g', '--gc_range')
        TM_FLAGS = ('-t', '--tm_range')
        THREADS_FLAGS = ('-n', '--num_threads')
        PCR_LEN_FLAGS = ('-r', '--pcr_prod_len')
        TM_DIFF_FLAGS = ('-d', '--tm_diff')
        HELP_FLAGS = ('-h', '--help')
        DEBUG_FLAGS = ("--debug",)
        SHORT_OPTS = INGROUP_FLAGS[0][-1] + ":" + \
                    OUT_FLAGS[0][-1] + ":" + \
                    OUTGROUP_FLAGS[0][-1] + ":" + \
                    DISALLOW_FLAGS[0][-1] + ":" + \
                    FMT_FLAGS[0][-1] + ":" + \
                    PRIMER_LEN_FLAGS[0][-1] + ":" + \
                    GC_FLAGS[0][-1] + ":" + \
                    TM_FLAGS[0][-1] + ":" + \
                    PCR_LEN_FLAGS[0][-1] + ":" + \
                    TM_DIFF_FLAGS[0][-1] + ":" + \
                    THREADS_FLAGS[0][-1] + ":" + \
                    HELP_FLAGS[0][-1]
        LONG_OPTS = (INGROUP_FLAGS[1][2:] + "=",
                    OUT_FLAGS[1][2:] + "=",
                    OUTGROUP_FLAGS[1][2:] + "=",
                    DISALLOW_FLAGS[1][2:] + "=",
                    FMT_FLAGS[1][2:] + "=",
                    PRIMER_LEN_FLAGS[1][2:] + "=",
                    GC_FLAGS[1][2:] + "=",
                    TM_FLAGS[1][2:] + "=",
                    PCR_LEN_FLAGS[1][2:] + "=",
                    TM_DIFF_FLAGS[1][2:] + "=",
                    THREADS_FLAGS[1][2:] + "=",
                    HELP_FLAGS[1][2:],
                    DEBUG_FLAGS[0][2:])

        # default values
        DEF_FRMT = ALLOWED_FORMATS[0]
        DEF_MIN_LEN = 16
        DEF_MAX_LEN = 20
        DEF_MIN_GC = 40.0
        DEF_MAX_GC = 60.0
        DEF_MIN_TM = 55.0
        DEF_MAX_TM = 68.0
        DEF_MIN_PCR = 120
        DEF_MAX_PCR = 2400
        DEF_MAX_TM_DIFF = 5.0
        DEF_NUM_THREADS = 1

        # messages
        IGNORE_MSG = 'ignoring unused argument: '
        ERR_MSG_1  = 'invalid or missing ingroup file(s)'
        ERR_MSG_2  = 'invalid or missing outgroup file(s)'
        ERR_MSG_3  = 'must specify exactly two integers for bad outgroup PCR product length'
        ERR_MSG_4  = 'bad outgroup PCR product lengths are not integers'
        ERR_MSG_5  = 'invalid format'
        ERR_MSG_6  = 'can only specify one primer length or a range (min,max)'
        ERR_MSG_7  = 'primer lengths are not integers'
        ERR_MSG_8  = 'must specify a range of GC values (min,max)'
        ERR_MSG_9  = 'gc values are not numeric'
        ERR_MSG_10 = 'must specify a range of Tm values (min, max)'
        ERR_MSG_11 = 'Tm values are not numeric'
        ERR_MSG_12 = 'can only specify one PCR product length or a range (min,max)'
        ERR_MSG_13 = 'PCR product lengths are not integers'
        ERR_MSG_14 = 'max Tm difference is not numeric'
        ERR_MSG_15 = 'num threads is not an integer'
        ERR_MSG_16 = 'must specify one or more ingroup files'
        ERR_MSG_17 = 'must specify an output file'

        def printHelp():
            GAP = " "*4
            EOL = "\n"
            SEP_1 = ", "
            SEP_2 = "|"
            SEP_3 = ","
            DEF_OPEN = ' (default: '
            CLOSE = ')'
            WIDTH = 30
            HELP_MSG = EOL + "Finds pairs of primers suitable for a group of input genomes" + EOL + \
                    GAP + "Joseph S. Wirth, 2023" + EOL*2 + \
                    "usage:" + EOL + \
                    GAP + "primerDesign.py [-iofpgtnrdh]" + EOL*2 + \
                    "required arguments:" + EOL + \
                    GAP + f'{INGROUP_FLAGS[0] + SEP_1 + INGROUP_FLAGS[1]:<{WIDTH}}[file] ingroup filename or a file pattern inside double-quotes (eg."*.gbff"){EOL}' + \
                    GAP + f"{OUT_FLAGS[0] + SEP_1 + OUT_FLAGS[1]:<{WIDTH}}[file] output filename{EOL*2}" + \
                    "optional arguments:" + EOL + \
                    GAP + f'{OUTGROUP_FLAGS[0] + SEP_1 + OUTGROUP_FLAGS[1]:<{WIDTH}}[file(s)] outgroup filename filename or a file pattern inside double-quotes (eg."*.gbff"){EOL}' + \
                    GAP + f"{DISALLOW_FLAGS[0] + SEP_1 + DISALLOW_FLAGS[1]:<{WIDTH}}[int,int] a range of PCR product lengths that the outgroup cannot produce{DEF_OPEN}same as '{PCR_LEN_FLAGS[1]}'{CLOSE}{EOL}" + \
                    GAP + f"{FMT_FLAGS[0] + SEP_1 + FMT_FLAGS[1]:<{WIDTH}}[str] file format of the ingroup and outgroup {ALLOWED_FORMATS[0]}{SEP_2}{ALLOWED_FORMATS[1]}{DEF_OPEN}{DEF_FRMT}{CLOSE}{EOL}" + \
                    GAP + f"{PRIMER_LEN_FLAGS[0] + SEP_1 + PRIMER_LEN_FLAGS[1]:<{WIDTH}}[int(s)] a single primer length or a range specified as 'min,max'{DEF_OPEN}{DEF_MIN_LEN}{SEP_3}{DEF_MAX_LEN}{CLOSE}{EOL}" + \
                    GAP + f"{GC_FLAGS[0] + SEP_1 + GC_FLAGS[1]:<{WIDTH}}[float,float] a min and max percent GC specified as a comma separated list{DEF_OPEN}{DEF_MIN_GC}{SEP_3}{DEF_MAX_GC}{CLOSE}{EOL}" + \
                    GAP + f"{TM_FLAGS[0] + SEP_1 + TM_FLAGS[1]:<{WIDTH}}[float,float] a min and max melting temp (Tm) specified as a comma separated list{DEF_OPEN}{DEF_MIN_TM}{SEP_3}{DEF_MAX_TM}{CLOSE}{EOL}" + \
                    GAP + f"{PCR_LEN_FLAGS[0] + SEP_1 + PCR_LEN_FLAGS[1]:<{WIDTH}}[int(s)] a single PCR product length or a range specified as 'min,max'{DEF_OPEN}{DEF_MIN_PCR}{SEP_3}{DEF_MAX_PCR}{CLOSE}{EOL}" + \
                    GAP + f"{TM_DIFF_FLAGS[0] + SEP_1 + TM_DIFF_FLAGS[1]:<{WIDTH}}[float] the maximum allowable Tm difference between a pair of primers{DEF_OPEN}{DEF_MAX_TM_DIFF}{CLOSE}{EOL}" + \
                    GAP + f"{THREADS_FLAGS[0] + SEP_1 + THREADS_FLAGS[1]:<{WIDTH}}[int] the number of threads for parallel processing{DEF_OPEN}{DEF_NUM_THREADS}{CLOSE}{EOL}" + \
                    GAP + f"{HELP_FLAGS[0] + SEP_1 + HELP_FLAGS[1]:<{WIDTH}}print this message{EOL*2}"
            
            print(HELP_MSG)
            
        # set default values
        self.ingroupFns = None
        self.outFn = None
        self.outgroupFns = list()
        self.disallowedLens = None # start as None; update at the end
        self.format = DEF_FRMT
        self.minLen = DEF_MIN_LEN
        self.maxLen = DEF_MAX_LEN
        self.minGc = DEF_MIN_GC
        self.maxGc = DEF_MAX_GC
        self.minTm = DEF_MIN_TM
        self.maxTm = DEF_MAX_TM
        self.minPcr = DEF_MIN_PCR
        self.maxPcr = DEF_MAX_PCR
        self.maxTmDiff = DEF_MAX_TM_DIFF
        self.numThreads = DEF_NUM_THREADS
        self.debug = False
        self.helpRequested = False
        
        # give help if requested
        if HELP_FLAGS[0] in sys.argv or HELP_FLAGS[1] in sys.argv or len(sys.argv) == 1:
            self.helpRequested = True
            printHelp()
        
        # parse command line arguments
        else:
            opts,args = getopt.getopt(sys.argv[1:], SHORT_OPTS, LONG_OPTS)
            for opt,arg in opts:
                # get the ingroup filenames
                if opt in INGROUP_FLAGS:
                    self.ingroupFns = glob.glob(arg)
                    for fn in self.ingroupFns:
                        if not os.path.isfile(fn):
                            raise ValueError(ERR_MSG_1)
                    if self.ingroupFns == []:
                        raise ValueError(ERR_MSG_1)
                
                # get output filehandle
                elif opt in OUT_FLAGS:
                    self.outFn = arg
                
                # get the outgroup filenames
                elif opt in OUTGROUP_FLAGS:
                    self.outgroupFns = glob.glob(arg)
                    for fn in self.outgroupFns:
                        if not os.path.isfile(fn):
                            raise ValueError(ERR_MSG_2)
                    if self.outgroupFns == []:
                        raise ValueError(ERR_MSG_2)
                
                # get the disallowed outgroup pcr product lengths
                elif opt in DISALLOW_FLAGS:
                    # split comma-separated list
                    disallowed = arg.split(SEP)
                    
                    # make sure exactly two values provided
                    if len(disallowed) != 2:
                        raise ValueError(ERR_MSG_3)
                    
                    # coerce lengths to ints
                    try:
                        disallowed = [int(x) for x in disallowed]
                    except:
                        raise ValueError(ERR_MSG_4)
                    
                    # save the values
                    self.disallowedLens = range(min(disallowed), max(disallowed)+1)
                
                # get the file format
                elif opt in FMT_FLAGS:
                    if arg not in ALLOWED_FORMATS:
                        raise ValueError(ERR_MSG_5)
                    self.format = arg
                
                # get the primer lengths
                elif opt in PRIMER_LEN_FLAGS:
                    # split comma-separated list
                    primerRange = arg.split(SEP)
                    
                    # make sure at one or two primers specified
                    if len(primerRange) not in {1,2}:
                        raise ValueError(ERR_MSG_6)
                    
                    # coerce to lengths to ints
                    try:
                        primerRange = [int(x) for x in primerRange]
                    except:
                        raise ValueError(ERR_MSG_7)
                    
                    # save values
                    self.minLen = min(primerRange)
                    self.maxLen = max(primerRange)
                
                # get the allowed GC range
                elif opt in GC_FLAGS:
                    # expecting two values separated by a comma
                    gcRange = arg.split(SEP)
                    if len(gcRange) != 2:
                        raise ValueError(ERR_MSG_8)
                    
                    # make sure the values are numeric
                    try:
                        gcRange = [float(x) for x in gcRange]
                    except:
                        raise ValueError(ERR_MSG_9)
                
                    # save values
                    self.minGc = min(gcRange)
                    self.maxGc = max(gcRange)
                
                # get the allowed Tm range
                elif opt in TM_FLAGS:
                    # expecting two values separated by a comma
                    tmRange = arg.split(SEP)
                    if len(tmRange) != 2:
                        raise ValueError(ERR_MSG_10)
                
                    # make sure the values are numeric
                    try:
                        tmRange = [float(x) for x in tmRange]
                    except:
                        raise ValueError(ERR_MSG_11)
                
                    # save values
                    self.minTm = min(tmRange)
                    self.maxTm = max(tmRange)
            
                # get the allowed PCR lengths
                elif opt in PCR_LEN_FLAGS:
                    # expecting one or two values
                    pcrRange = arg.split(SEP)
                    if len(pcrRange) not in {1,2}:
                        raise ValueError(ERR_MSG_12)
                    
                    # coerce to integers
                    try:
                        pcrRange = [int(x) for x in pcrRange]
                    except:
                        raise ValueError(ERR_MSG_13)
                
                    # save values
                    self.minPcr = min(pcrRange)
                    self.maxPcr = max(pcrRange)
                    
                    # see if the disallowed needs to be changed
                    if self.disallowedLens is None:
                        self.disallowedLens = range(self.minPcr,self.maxPcr+1)
                
                # get the allowed Tm difference between primer pairs
                elif opt in TM_DIFF_FLAGS:
                    # make sure input is numeric
                    try:
                        self.maxTmDiff = float(arg)
                    except:
                        raise ValueError(ERR_MSG_14)
                
                # get the number of threads to use
                elif opt in THREADS_FLAGS:
                    # make sure input is an integer
                    try:
                        self.numThreads = int(arg)
                    except:
                        raise ValueError(ERR_MSG_15)
                
                elif opt in DEBUG_FLAGS:
                    self.debug = True
                
                else:
                    print(IGNORE_MSG + opt + " " + arg)
            
            # update disallowed to match pcr parameters unless it was already specified
            if self.disallowedLens is None:
                self.disallowedLens = range(self.minPcr, self.maxPcr + 1)
            
            # make sure an input file was specified
            if self.ingroupFns is None:
                raise ValueError(ERR_MSG_16)
            
            # make sure an output file was specified
            if self.outFn is None:
                raise ValueError(ERR_MSG_17)

    def __debugSetup(self):
        self.debugger = Debugger()
    
    def logRunDetails(self, name:str, version:str) -> None:
        WIDTH = 32
        self.debugger.writeDebugMsg(name, f'{"version:":<{WIDTH}}{version}')
        self.debugger.writeDebugMsg(name, f'{"ingroup:":<{WIDTH}}{",".join(self.ingroupFns)}')
        self.debugger.writeDebugMsg(name, f'{"outgroup:":<{WIDTH}}{",".join(self.outgroupFns)}')
        self.debugger.writeDebugMsg(name, f'{"out filename:":{WIDTH}}{self.outFn}')
        self.debugger.writeDebugMsg(name, f'{"file format:":{WIDTH}}{self.format}')
        self.debugger.writeDebugMsg(name, f'{"min kmer len:":{WIDTH}}{self.minLen}')
        self.debugger.writeDebugMsg(name, f'{"max kmer len:":<{WIDTH}}{self.maxLen}')
        self.debugger.writeDebugMsg(name, f'{"min % G+C":<{WIDTH}}{self.minGc}')
        self.debugger.writeDebugMsg(name, f'{"max % G+C":<{WIDTH}}{self.maxGc}')
        self.debugger.writeDebugMsg(name, f'{"min Tm:":<{WIDTH}}{self.minTm}')
        self.debugger.writeDebugMsg(name, f'{"max Tm:":<{WIDTH}}{self.maxTm}')
        self.debugger.writeDebugMsg(name, f'{"max Tm difference:":<{WIDTH}}{self.maxTmDiff}')
        self.debugger.writeDebugMsg(name, f'{"min PCR len:":<{WIDTH}}{self.minPcr}')
        self.debugger.writeDebugMsg(name, f'{"max PCR len:":<{WIDTH}}{self.maxPcr}')
        self.debugger.writeDebugMsg(name, f'{"disallowed outgroup PCR lens:":<{WIDTH}}{",".join(map(str,[min(self.disallowedLens),max(self.disallowedLens)]))}')
        self.debugger.writeDebugMsg(name, f'{"num threads:":<{WIDTH}}{self.numThreads}')
    
    def dumpObj(self, obj, fn:str) -> None:
        with open(os.path.join(self.debugger.debugDir, fn)) as fh:
            pickle.dump(obj, fh)
    