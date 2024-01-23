from Bio import SeqIO
import getopt, glob, os, sys
from bin.Primer import Primer
from Bio.SeqRecord import SeqRecord
from bin.getPrimerPairs import _getPrimerPairs
from bin.Clock import Clock, _printStart, _printDone
from bin.getCandidateKmers import _getAllCandidateKmers
from bin.removeOutgroupPrimers import _removeOutgroupPrimers


def __parseArgs() -> tuple[list[str],list[str],str,str,int,int,float,float,float,float,int,int,float,range,int,bool]:
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
        tuple[list[str],list[str],str,str,int,int,float,float,float,float,int,int,float,range,int,bool]:
            ingroupFns,outgroupFns,outFN,format,minPrimerLen,maxPrimerLen,minGc,maxGc,minTm,maxTm,minPcrLen,maxPcrLen,maxTmDiff,disallowedLens,numThreads,helpRequested
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
                 HELP_FLAGS[1][2:])

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
        HELP_MSG = EOL + "Finds pairs of primers suitable for an input genome." + EOL + \
                   GAP + "Joseph S. Wirth, 2023" + EOL*2 + \
                   "usage:" + EOL + \
                   GAP + "primerDesign.py [-iofpgtnrdh]" + EOL*2 + \
                   "required arguments:" + EOL + \
                   GAP + f'{INGROUP_FLAGS[0] + SEP_1 + INGROUP_FLAGS[1]:<30}[file] ingroup filename or a file pattern inside double-quotes (eg."*.gbff"){EOL}' + \
                   GAP + f"{OUT_FLAGS[0] + SEP_1 + OUT_FLAGS[1]:<30}[file] output filename{EOL*2}" + \
                   "optional arguments:" + EOL + \
                   GAP + f'{OUTGROUP_FLAGS[0] + SEP_1 + OUTGROUP_FLAGS[1]:<30}[file(s)] outgroup filename filename or a file pattern inside double-quotes (eg."*.gbff"){EOL}' + \
                   GAP + f"{DISALLOW_FLAGS[0] + SEP_1 + DISALLOW_FLAGS[1]:<30}[int,int] a range of PCR product lengths that the outgroup cannot produce{DEF_OPEN}same as '{PCR_LEN_FLAGS[1]}'{CLOSE}{EOL}" + \
                   GAP + f"{FMT_FLAGS[0] + SEP_1 + FMT_FLAGS[1]:<30}[str] file format of the ingroup and outgroup {ALLOWED_FORMATS[0]}{SEP_2}{ALLOWED_FORMATS[1]}{DEF_OPEN}{DEF_FRMT}{CLOSE}{EOL}" + \
                   GAP + f"{PRIMER_LEN_FLAGS[0] + SEP_1 + PRIMER_LEN_FLAGS[1]:<30}[int(s)] a single primer length or a range specified as 'min,max'{DEF_OPEN}{DEF_MIN_LEN}{SEP_3}{DEF_MAX_LEN}{CLOSE}{EOL}" + \
                   GAP + f"{GC_FLAGS[0] + SEP_1 + GC_FLAGS[1]:<30}[float,float] a min and max percent GC specified as a comma separated list{DEF_OPEN}{DEF_MIN_GC}{SEP_3}{DEF_MAX_GC}{CLOSE}{EOL}" + \
                   GAP + f"{TM_FLAGS[0] + SEP_1 + TM_FLAGS[1]:<30}[float,float] a min and max melting temp (Tm) specified as a comma separated list{DEF_OPEN}{DEF_MIN_TM}{SEP_3}{DEF_MAX_TM}{CLOSE}{EOL}" + \
                   GAP + f"{PCR_LEN_FLAGS[0] + SEP_1 + PCR_LEN_FLAGS[1]:<30}[int(s)] a single PCR product length or a range specified as 'min,max'{DEF_OPEN}{DEF_MIN_PCR}{SEP_3}{DEF_MAX_PCR}{CLOSE}{EOL}" + \
                   GAP + f"{TM_DIFF_FLAGS[0] + SEP_1 + TM_DIFF_FLAGS[1]:<30}[float] the maximum allowable Tm difference between a pair of primers{DEF_OPEN}{DEF_MAX_TM_DIFF}{CLOSE}{EOL}" + \
                   GAP + f"{THREADS_FLAGS[0] + SEP_1 + THREADS_FLAGS[1]:<30}[int] the number of threads for parallel processing{DEF_OPEN}{DEF_NUM_THREADS}{CLOSE}{EOL}" + \
                   GAP + f"{HELP_FLAGS[0] + SEP_1 + HELP_FLAGS[1]:<30}print this message{EOL*2}"
        
        print(HELP_MSG)
        
    # set default values
    ingroupFns = None
    outFN = None
    outgroupFns = list()
    disallowedLens = None # start as None; update at the end
    frmt = DEF_FRMT
    minLen = DEF_MIN_LEN
    maxLen = DEF_MAX_LEN
    minGc = DEF_MIN_GC
    maxGc = DEF_MAX_GC
    minTm = DEF_MIN_TM
    maxTm = DEF_MAX_TM
    minPcr = DEF_MIN_PCR
    maxPcr = DEF_MAX_PCR
    maxTmDiff = DEF_MAX_TM_DIFF
    numThreads = DEF_NUM_THREADS
    helpRequested = False
    
    # give help if requested
    if HELP_FLAGS[0] in sys.argv or HELP_FLAGS[1] in sys.argv or len(sys.argv) == 1:
        helpRequested = True
        printHelp()
    
    # parse command line arguments
    else:
        opts,args = getopt.getopt(sys.argv[1:], SHORT_OPTS, LONG_OPTS)
        for opt,arg in opts:
            # get the ingroup filenames
            if opt in INGROUP_FLAGS:
                ingroupFns = glob.glob(arg)
                for fn in ingroupFns:
                    if not os.path.isfile(fn):
                        raise ValueError(ERR_MSG_1)
                if ingroupFns == []:
                    raise ValueError(ERR_MSG_1)
            
            # get output filehandle
            elif opt in OUT_FLAGS:
                outFN = arg
            
            # get the outgroup filenames
            elif opt in OUTGROUP_FLAGS:
                outgroupFns = glob.glob(arg)
                for fn in arg:
                    if not os.path.isfile(fn):
                        raise ValueError(ERR_MSG_2)
                if outgroupFns == []:
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
                disallowedLens = range(min(disallowed), max(disallowed)+1)
            
            # get the file format
            elif opt in FMT_FLAGS:
                if arg not in ALLOWED_FORMATS:
                    raise ValueError(ERR_MSG_5)
                frmt = arg
            
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
                minLen = min(primerRange)
                maxLen = max(primerRange)
            
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
                minGc = min(gcRange)
                maxGc = max(gcRange)
            
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
                minTm = min(tmRange)
                maxTm = max(tmRange)
        
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
                minPcr = min(pcrRange)
                maxPcr = max(pcrRange)
                
                # see if the disallowed needs to be changed
                if disallowedLens is None:
                    disallowedLens = range(minPcr,maxPcr+1)
            
            # get the allowed Tm difference between primer pairs
            elif opt in TM_DIFF_FLAGS:
                # make sure input is numeric
                try:
                    maxTmDiff = float(arg)
                except:
                    raise ValueError(ERR_MSG_14)
            
            # get the number of threads to use
            elif opt in THREADS_FLAGS:
                # make sure input is an integer
                try:
                    numThreads = int(arg)
                except:
                    raise ValueError(ERR_MSG_15)
            
            else:
                print(IGNORE_MSG + opt + " " + arg)
        
        # update disallowed to match pcr parameters unless it was already specified
        if disallowedLens is None:
            disallowedLens = range(minPcr, maxPcr + 1)
        
        # make sure an input file was specified
        if ingroupFns is None:
            raise ValueError(ERR_MSG_16)
        
        # make sure an output file was specified
        if outFN is None:
            raise ValueError(ERR_MSG_17)
    
    return ingroupFns,outgroupFns,outFN,frmt,minLen,maxLen,minGc,maxGc,minTm,maxTm,minPcr,maxPcr,maxTmDiff,disallowedLens,numThreads,helpRequested


def __readSequenceData(seqFiles:list[str], frmt:str) -> dict[str, list[SeqRecord]]:
    """reads sequence data into file

    Args:
        seqFiles (list[str]): a list of sequence files to read
        frmt (str): the format of the sequence files

    Returns:
        dict[str, list[SeqRecord]]: key=genome name; val=list of contigs as SeqRecords
    """
    # initialize output
    out = dict()
    
    # for each file in the list
    for fn in seqFiles:
        # get the genome name and use it as a key to store the list of parsed contigs
        name = os.path.basename(fn)
        out[name] = SeqIO.parse(fn, frmt)
    
    return out


def __writePrimerPairs(fn:str, pairs:dict[tuple[Primer,Primer],dict[str,int]]) -> None:
    """writes pairs of primers to file

    Args:
        fn (str): the filename to write
        pairs (dict[tuple[Primer,Primer],dict[str,int]]): key=primer pair; val=dict: key=genome name; val=pcr product length
    """
    # contants
    EOL = "\n"
    SEP = "\t"
    NUM_DEC = 1
    
    # helper function to create the headers
    def getHeaders(names) -> list[str]:
        # constants
        HEADERS = ('fwd_seq',
                   'fwd_Tm',
                   'fwd_GC',
                   'rev_seq',
                   'rev_Tm',
                   'rev_GC')
        CONTIG = "_contig"
        LENGTH = "_length"
        
        # each name will have a contig and a length
        headers = list(HEADERS)
        for name in names:
            headers.append(name + CONTIG)
            headers.append(name + LENGTH)
        
        return headers
    
    names = list(next(iter(pairs.values())).keys())
    headers = getHeaders(names)
    
    # open the file
    with open(fn, 'w') as fh:
        # write the headers
        fh.write(SEP.join(headers) + EOL)
        fh.flush()
        
        # for each primer pair
        for fwd,rev in pairs.keys():
            # save the primer pair data
            row = [fwd.seq,
                   round(fwd.Tm, NUM_DEC),
                   round(fwd.gcPer, NUM_DEC),
                   rev.seq,
                   round(rev.Tm, NUM_DEC),
                   round(rev.gcPer, NUM_DEC)]
            
            # then save the contig name and PCR product length for each genome
            for name in names:
                row.extend(pairs[(fwd,rev)][name])
            
            fh.write(SEP.join(map(str, row)) + EOL)
            fh.flush()


def _main() -> None:
    """main runner function:
        * reads ingroup and outgroup sequences into memory
        * gets candidate kmers to use to search for primer pairs
        * gets primer pairs that are present in all ingroup and excluded from all outgroup
    """
    # messages
    MSG_1  = "identifying kmers suitable for use as primers"
    MSG_2  = "identifying primer pairs suitable for use in PCR"
    MSG_3A = "    identified "
    MSG_3B = " candidate primer pairs shared in all ingroup sequences"
    MSG_4  = "removing primer pairs present in the outgroup sequences"
    MSG_5A = "writing "
    MSG_5B = " primer pairs to file"
    MSG_6  = '\ntotal runtime: '

    # parse command line arguments
    ingroupFiles,outgroupFiles,outFn,frmt,minPrimerLen,maxPrimerLen,minGc,maxGc,minTm,maxTm,minPcrLen,maxPcrLen,maxTmDiff,disallowedLens,numThreads,helpRequested = __parseArgs()
    
    # start the timers
    totalClock = Clock()
    clock = Clock()
    
    # only do work if help was not requested
    if not helpRequested:
        # read the ingroup and outgroup sequences into memory
        ingroupSeqs = __readSequenceData(ingroupFiles, frmt)

        # get the candidate kmers for the ingroup
        _printStart(clock, MSG_1, '\n')
        candidateKmers = _getAllCandidateKmers(ingroupSeqs, minPrimerLen, maxPrimerLen, minGc, maxGc, minTm, maxTm, numThreads)
        _printDone(clock)
        del ingroupSeqs
        
        # get the suitable primer pairs for the ingroup
        _printStart(clock, MSG_2)
        pairs = _getPrimerPairs(candidateKmers, minPrimerLen, minPcrLen, maxPcrLen, maxTmDiff, numThreads)
        _printDone(clock)
        del candidateKmers
        
        # print the number of candidate primer pairs
        print(f"{MSG_3A}{len(pairs)}{MSG_3B}")
        
        # remove primers that make products in the outgroup
        _printStart(clock, MSG_4)
        outgroupSeqs = __readSequenceData(outgroupFiles, frmt)
        _removeOutgroupPrimers(outgroupSeqs, pairs, disallowedLens, numThreads)   ######## what is that final parameter??
        _printDone(clock)
        del outgroupSeqs
        
        # write results to file
        _printStart(clock, f"{MSG_5A}{len(pairs)}{MSG_5B}")
        __writePrimerPairs(outFn, pairs)
        _printDone(clock)
        
        # print the total runtime
        print(MSG_6, end='', flush=True)
        totalClock.printTime()
